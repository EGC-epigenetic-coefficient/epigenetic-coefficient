if (!exists("wasserstein_1_transport", mode = "function")) {
  stop("source egc_reference_R.R before egc_inference_R.R")
}

egc_validate_inference_inputs <- function(burden, groups, group_ranks, weights = NULL) {
  burden <- as.numeric(burden)
  groups <- as.character(groups)
  if (is.null(weights)) weights <- rep(1, length(burden))
  weights <- as.numeric(weights)
  if (length(burden) == 0L || length(groups) != length(burden) || length(weights) != length(burden)) {
    stop("invalid input length")
  }
  if (any(!is.finite(burden)) || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("burden and strictly positive weights must be finite")
  }
  observed <- unique(groups)
  if (!all(observed %in% names(group_ranks))) stop("missing group rank")
  ranks <- as.numeric(group_ranks[observed])
  names(ranks) <- observed
  if (any(!is.finite(ranks)) || any(ranks < 0) || any(ranks > 1) || length(unique(ranks)) < 2L) {
    stop("invalid group ranks")
  }
  list(burden = burden, groups = groups, weights = weights, ranks = ranks)
}

kish_effective_n <- function(weights) {
  sum(weights)^2 / sum(weights^2)
}

probability_superiority_fast <- function(x, y, x_weights = NULL, y_weights = NULL) {
  if (is.null(x_weights)) x_weights <- rep(1, length(x))
  if (is.null(y_weights)) y_weights <- rep(1, length(y))
  x <- as.numeric(x)
  y <- as.numeric(y)
  x_weights <- as.numeric(x_weights) / sum(x_weights)
  y_weights <- as.numeric(y_weights) / sum(y_weights)
  oy <- order(y, method = "radix")
  sy <- y[oy]
  sw <- y_weights[oy]
  cw <- cumsum(sw)
  n_le <- findInterval(x, sy)
  n_lt <- findInterval(x, sy, left.open = TRUE)
  p_le <- ifelse(n_le > 0L, cw[pmax(n_le, 1L)], 0)
  p_lt <- ifelse(n_lt > 0L, cw[pmax(n_lt, 1L)], 0)
  p_gt <- 1 - p_le
  as.numeric(sum(x_weights * (p_lt - p_gt)))
}

estimate_egc_inference_statistics <- function(burden, groups, group_ranks, weights = NULL) {
  checked <- egc_validate_inference_inputs(burden, groups, group_ranks, weights)
  burden <- checked$burden
  groups <- checked$groups
  weights <- checked$weights
  ranks <- checked$ranks
  ordered_groups <- names(sort(ranks))
  total_weight <- sum(weights)
  probabilities <- setNames(numeric(length(ordered_groups)), ordered_groups)
  values <- vector("list", length(ordered_groups)); names(values) <- ordered_groups
  group_weights <- vector("list", length(ordered_groups)); names(group_weights) <- ordered_groups
  group_n <- setNames(integer(length(ordered_groups)), ordered_groups)
  group_neff <- setNames(numeric(length(ordered_groups)), ordered_groups)
  for (label in ordered_groups) {
    mask <- groups == label
    values[[label]] <- burden[mask]
    group_weights[[label]] <- weights[mask]
    probabilities[[label]] <- sum(weights[mask]) / total_weight
    group_n[[label]] <- sum(mask)
    group_neff[[label]] <- kish_effective_n(weights[mask])
  }

  dm_u <- 0
  dm_r_num <- 0
  dm_r_den <- 0
  A <- 0
  D <- 0
  pairwise <- list()
  pair_index <- 1L
  for (lower_index in seq_len(length(ordered_groups) - 1L)) {
    lower <- ordered_groups[[lower_index]]
    for (higher_index in seq.int(lower_index + 1L, length(ordered_groups))) {
      higher <- ordered_groups[[higher_index]]
      pair_weight <- probabilities[[lower]] * probabilities[[higher]]
      rank_distance <- ranks[[higher]] - ranks[[lower]]
      w1 <- wasserstein_1_transport(
        values[[lower]], values[[higher]], group_weights[[lower]], group_weights[[higher]]
      )
      delta <- probability_superiority_fast(
        values[[lower]], values[[higher]], group_weights[[lower]], group_weights[[higher]]
      )
      orientation_weight <- pair_weight * rank_distance
      dm_u <- dm_u + 2 * pair_weight * w1
      dm_r_num <- dm_r_num + pair_weight * rank_distance * w1
      dm_r_den <- dm_r_den + pair_weight * rank_distance^2
      A <- A + orientation_weight * delta
      D <- D + orientation_weight * abs(delta)
      pairwise[[pair_index]] <- data.frame(
        lower_group = lower,
        higher_group = higher,
        pair_weight = pair_weight,
        rank_distance = rank_distance,
        w1 = w1,
        delta = delta,
        orientation_weight = orientation_weight,
        stringsAsFactors = FALSE
      )
      pair_index <- pair_index + 1L
    }
  }
  dd <- if (D <= 1e-14) NA_real_ else A / D
  coherence <- if (D <= 1e-14) NA_real_ else abs(A) / D
  list(
    statistics = c(
      egc_dm_u = dm_u,
      egc_dm_r = dm_r_num / dm_r_den,
      orientation_A = A,
      orientation_D = D,
      egc_dd = dd,
      directional_coherence = coherence
    ),
    group_n = group_n,
    group_neff = group_neff,
    pairwise = do.call(rbind, pairwise)
  )
}

egc_permutation_test <- function(burden, groups, group_ranks, weights = NULL, B = 499L) {
  checked <- egc_validate_inference_inputs(burden, groups, group_ranks, weights)
  observed <- estimate_egc_inference_statistics(
    checked$burden, checked$groups, checked$ranks, checked$weights
  )$statistics
  permuted <- matrix(NA_real_, nrow = B, ncol = 4L)
  colnames(permuted) <- c("egc_dm_u", "orientation_A", "orientation_D", "egc_dd")
  for (b in seq_len(B)) {
    perm_groups <- sample(checked$groups, length(checked$groups), replace = FALSE)
    stat <- estimate_egc_inference_statistics(
      checked$burden, perm_groups, checked$ranks, checked$weights
    )$statistics
    permuted[b, ] <- stat[c("egc_dm_u", "orientation_A", "orientation_D", "egc_dd")]
  }
  p_dm <- (1 + sum(permuted[, "egc_dm_u"] >= observed[["egc_dm_u"]])) / (B + 1)
  p_D <- (1 + sum(permuted[, "orientation_D"] >= observed[["orientation_D"]])) / (B + 1)
  list(
    observed = observed,
    permuted = permuted,
    p_dm = p_dm,
    p_D = p_D,
    null_mean_dm = mean(permuted[, "egc_dm_u"]),
    null_median_dm = median(permuted[, "egc_dm_u"]),
    null_sd_dm = sd(permuted[, "egc_dm_u"]),
    null_quantiles_dm = quantile(permuted[, "egc_dm_u"], c(0.5, 0.9, 0.95, 0.99), names = FALSE),
    null_quantile95_D = unname(quantile(permuted[, "orientation_D"], 0.95))
  )
}

egc_stratified_resample_indices <- function(groups, fraction = 1, replace = TRUE, minimum_per_group = 1L) {
  groups <- as.character(groups)
  labels <- unique(groups)
  result <- integer()
  for (label in labels) {
    idx <- which(groups == label)
    target <- if (fraction == 1) length(idx) else floor(length(idx) * fraction)
    target <- max(minimum_per_group, target)
    if (!replace) target <- min(target, length(idx))
    if (target <= 0L || (!replace && target > length(idx))) stop("invalid resample size")
    result <- c(result, sample(idx, target, replace = replace))
  }
  result
}

egc_stratified_bootstrap <- function(
    burden, groups, group_ranks, weights = NULL, B = 499L,
    fraction = 1, replace = TRUE, minimum_per_group = 1L) {
  checked <- egc_validate_inference_inputs(burden, groups, group_ranks, weights)
  observed <- estimate_egc_inference_statistics(
    checked$burden, checked$groups, checked$ranks, checked$weights
  )$statistics
  boot <- matrix(NA_real_, nrow = B, ncol = length(observed))
  colnames(boot) <- names(observed)
  failures <- character()
  for (b in seq_len(B)) {
    idx <- tryCatch(
      egc_stratified_resample_indices(
        checked$groups, fraction = fraction, replace = replace,
        minimum_per_group = minimum_per_group
      ),
      error = function(e) e
    )
    if (inherits(idx, "error")) {
      failures <- c(failures, conditionMessage(idx))
      next
    }
    stat <- tryCatch(
      estimate_egc_inference_statistics(
        checked$burden[idx], checked$groups[idx], checked$ranks, checked$weights[idx]
      )$statistics,
      error = function(e) e
    )
    if (inherits(stat, "error")) {
      failures <- c(failures, conditionMessage(stat))
      next
    }
    boot[b, ] <- stat
  }
  valid <- is.finite(boot[, "egc_dm_u"])
  boot_valid <- boot[valid, , drop = FALSE]
  percentile <- apply(boot_valid, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  if (is.null(dim(percentile))) percentile <- matrix(percentile, nrow = 2)
  colnames(percentile) <- colnames(boot_valid)
  sign_A <- sign(observed[["orientation_A"]])
  sign_stability <- if (sign_A == 0) NA_real_ else mean(sign(boot_valid[, "orientation_A"]) == sign_A)
  list(
    observed = observed,
    bootstrap = boot,
    percentile_ci = percentile,
    valid_replicates = sum(valid),
    failure_fraction = 1 - mean(valid),
    failure_messages = sort(table(failures), decreasing = TRUE),
    sign_stability_A = sign_stability,
    fraction_D_below_tolerance = mean(boot_valid[, "orientation_D"] <= 1e-12),
    bootstrap_bias_corrected_dm = max(0, 2 * observed[["egc_dm_u"]] - mean(boot_valid[, "egc_dm_u"]))
  )
}

egc_m_out_of_n_interval <- function(
    burden, groups, group_ranks, weights = NULL, B = 499L, gamma = 0.7,
    replace = TRUE, minimum_per_group = 10L) {
  n <- length(burden)
  m <- floor(n^gamma)
  fraction <- m / n
  group_counts <- as.numeric(table(groups))
  realized_group_m <- pmax(minimum_per_group, floor(group_counts * fraction))
  if (!replace) realized_group_m <- pmin(realized_group_m, group_counts)
  realized_m <- sum(realized_group_m)
  res <- egc_stratified_bootstrap(
    burden, groups, group_ranks, weights, B = B, fraction = fraction,
    replace = replace, minimum_per_group = minimum_per_group
  )
  valid <- is.finite(res$bootstrap[, "egc_dm_u"])
  roots <- sqrt(realized_m) * (
    res$bootstrap[valid, "egc_dm_u"] - res$observed[["egc_dm_u"]]
  )
  q <- quantile(roots, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  ci <- c(
    lower = max(0, res$observed[["egc_dm_u"]] - q[[2L]] / sqrt(n)),
    upper = max(0, res$observed[["egc_dm_u"]] - q[[1L]] / sqrt(n))
  )
  list(
    observed = res$observed,
    gamma = gamma,
    nominal_m = m,
    realized_m = realized_m,
    replace = replace,
    ci_dm_u = ci,
    failure_fraction = res$failure_fraction
  )
}

egc_bias_candidates <- function(permutation_result, bootstrap_result) {
  raw <- permutation_result$observed[["egc_dm_u"]]
  perm <- permutation_result$permuted[, "egc_dm_u"]
  boot <- bootstrap_result$bootstrap[, "egc_dm_u"]
  boot <- boot[is.finite(boot)]
  c(
    raw = raw,
    permutation_mean_subtraction = max(0, raw - mean(perm)),
    permutation_median_subtraction = max(0, raw - median(perm)),
    within_stratum_bootstrap_bias_correction = max(0, 2 * raw - mean(boot)),
    quadrature_null_correction = sqrt(max(0, raw^2 - mean(perm^2)))
  )
}

egc_subsampling_extrapolation <- function(
    burden, groups, group_ranks, weights = NULL,
    fractions = c(0.5, 0.7, 0.9), K = 199L) {
  checked <- egc_validate_inference_inputs(burden, groups, group_ranks, weights)
  rows <- list()
  pos <- 1L
  for (fraction in fractions) {
    values <- numeric(K)
    realized_m <- integer(K)
    for (k in seq_len(K)) {
      idx <- egc_stratified_resample_indices(
        checked$groups, fraction = fraction, replace = FALSE, minimum_per_group = 1L
      )
      values[[k]] <- estimate_egc_inference_statistics(
        checked$burden[idx], checked$groups[idx], checked$ranks, checked$weights[idx]
      )$statistics[["egc_dm_u"]]
      realized_m[[k]] <- length(idx)
    }
    rows[[pos]] <- data.frame(
      fraction = fraction,
      m = round(mean(realized_m)),
      mean_dm_u = mean(values),
      sd_dm_u = sd(values)
    )
    pos <- pos + 1L
  }
  profile <- do.call(rbind, rows)
  fit <- lm(mean_dm_u ~ I(1 / sqrt(m)), data = profile)
  theta <- unname(coef(fit)[[1L]])
  c_term <- unname(coef(fit)[[2L]])
  r_squared <- summary(fit)$r.squared
  admissible <- is.finite(theta) && is.finite(c_term) && c_term >= 0 && r_squared >= 0.8
  list(
    estimate = if (admissible) max(0, theta) else NA_real_,
    theta_untruncated = theta,
    c = c_term,
    r_squared = r_squared,
    admissible = admissible,
    profile = profile
  )
}
