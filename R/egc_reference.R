weighted_mean <- function(values, weights) {
  sum(weights * values) / sum(weights)
}

weighted_covariance <- function(x, y, weights) {
  x_bar <- weighted_mean(x, weights)
  y_bar <- weighted_mean(y, weights)
  sum(weights * (x - x_bar) * (y - y_bar)) / sum(weights)
}

weighted_correlation <- function(x, y, weights, tolerance = 1e-14) {
  var_x <- weighted_covariance(x, x, weights)
  var_y <- weighted_covariance(y, y, weights)
  if (var_x <= tolerance || var_y <= tolerance) return(NA_real_)
  weighted_covariance(x, y, weights) / sqrt(var_x * var_y)
}

wasserstein_1_transport <- function(x, y, x_weights = NULL, y_weights = NULL) {
  if (is.null(x_weights)) x_weights <- rep(1, length(x))
  if (is.null(y_weights)) y_weights <- rep(1, length(y))
  stopifnot(length(x) > 0, length(y) > 0, length(x_weights) == length(x), length(y_weights) == length(y))
  if (any(x_weights <= 0) || any(y_weights <= 0)) stop("W1 weights must be strictly positive")
  order_x <- order(x, method = "radix")
  order_y <- order(y, method = "radix")
  sorted_x <- x[order_x]
  sorted_y <- y[order_y]
  probability_x <- x_weights[order_x] / sum(x_weights)
  probability_y <- y_weights[order_y] / sum(y_weights)
  i <- 1L
  j <- 1L
  remaining_x <- probability_x[[1L]]
  remaining_y <- probability_y[[1L]]
  distance <- 0
  tolerance <- 1e-15
  while (i <= length(sorted_x) && j <= length(sorted_y)) {
    mass <- min(remaining_x, remaining_y)
    distance <- distance + mass * abs(sorted_x[[i]] - sorted_y[[j]])
    remaining_x <- remaining_x - mass
    remaining_y <- remaining_y - mass
    if (remaining_x <= tolerance) {
      i <- i + 1L
      if (i <= length(sorted_x)) remaining_x <- probability_x[[i]]
    }
    if (remaining_y <= tolerance) {
      j <- j + 1L
      if (j <= length(sorted_y)) remaining_y <- probability_y[[j]]
    }
  }
  as.numeric(distance)
}

probability_superiority_difference <- function(x, y, x_weights = NULL, y_weights = NULL) {
  if (is.null(x_weights)) x_weights <- rep(1, length(x))
  if (is.null(y_weights)) y_weights <- rep(1, length(y))
  probability_x <- x_weights / sum(x_weights)
  probability_y <- y_weights / sum(y_weights)
  comparison <- outer(x, y, FUN = function(a, b) sign(a - b))
  pair_weights <- outer(probability_x, probability_y)
  as.numeric(sum(pair_weights * comparison))
}

estimate_egc_components <- function(burden, groups, group_ranks, weights = NULL) {
  burden <- as.numeric(burden)
  groups <- as.character(groups)
  if (is.null(weights)) weights <- rep(1, length(burden))
  weights <- as.numeric(weights)
  if (length(burden) == 0 || length(groups) != length(burden) || length(weights) != length(burden)) stop("invalid input length")
  if (any(!is.finite(burden)) || any(!is.finite(weights)) || any(weights <= 0)) stop("burden and positive weights must be finite")
  observed <- unique(groups)
  if (!all(observed %in% names(group_ranks))) stop("missing group rank")
  ranks <- as.numeric(group_ranks[observed])
  names(ranks) <- observed
  if (any(!is.finite(ranks)) || any(ranks < 0) || any(ranks > 1) || length(unique(ranks)) < 2) stop("invalid ranks")

  individual_ranks <- unname(ranks[groups])
  rank_variance <- weighted_covariance(individual_ranks, individual_ranks, weights)
  if (rank_variance <= 1e-14) stop("rank variance must be positive")
  beta <- weighted_covariance(individual_ranks, burden, weights) / rank_variance
  egc_s <- -beta

  ordered_groups <- names(sort(ranks))
  total_weight <- sum(weights)
  probabilities <- setNames(numeric(length(ordered_groups)), ordered_groups)
  means <- probabilities
  values <- vector("list", length(ordered_groups)); names(values) <- ordered_groups
  group_weights <- vector("list", length(ordered_groups)); names(group_weights) <- ordered_groups
  for (label in ordered_groups) {
    mask <- groups == label
    values[[label]] <- burden[mask]
    group_weights[[label]] <- weights[mask]
    probabilities[[label]] <- sum(weights[mask]) / total_weight
    means[[label]] <- weighted_mean(burden[mask], weights[mask])
  }

  mean_burden <- weighted_mean(burden, weights)
  egc_lm <- sqrt(sum(probabilities * (means - mean_burden)^2))
  egc_ld <- -weighted_correlation(unname(ranks[ordered_groups]), unname(means), unname(probabilities))

  dm_u <- 0
  dm_r_numerator <- 0
  dm_r_denominator <- 0
  dd_numerator <- 0
  dd_denominator <- 0
  pairwise <- list()
  pair_index <- 1L
  for (lower_index in seq_len(length(ordered_groups) - 1L)) {
    lower <- ordered_groups[[lower_index]]
    for (higher_index in seq.int(lower_index + 1L, length(ordered_groups))) {
      higher <- ordered_groups[[higher_index]]
      pair_weight <- probabilities[[lower]] * probabilities[[higher]]
      rank_distance <- ranks[[higher]] - ranks[[lower]]
      w1 <- wasserstein_1_transport(values[[lower]], values[[higher]], group_weights[[lower]], group_weights[[higher]])
      delta <- probability_superiority_difference(values[[lower]], values[[higher]], group_weights[[lower]], group_weights[[higher]])
      dm_u <- dm_u + 2 * pair_weight * w1
      dm_r_numerator <- dm_r_numerator + pair_weight * rank_distance * w1
      dm_r_denominator <- dm_r_denominator + pair_weight * rank_distance^2
      orientation_weight <- pair_weight * rank_distance
      dd_numerator <- dd_numerator + orientation_weight * delta
      dd_denominator <- dd_denominator + orientation_weight * abs(delta)
      pairwise[[pair_index]] <- data.frame(
        lower_group = lower,
        higher_group = higher,
        p_lower = probabilities[[lower]],
        p_higher = probabilities[[higher]],
        rank_distance = rank_distance,
        w1 = w1,
        delta = delta,
        stringsAsFactors = FALSE
      )
      pair_index <- pair_index + 1L
    }
  }
  if (dm_r_denominator <= 1e-14) stop("DM-R denominator must be positive")
  egc_dd <- if (dd_denominator <= 1e-14) NA_real_ else dd_numerator / dd_denominator
  list(
    components = c(
      egc_s = egc_s,
      egc_lm = egc_lm,
      egc_ld = egc_ld,
      egc_dm_u = dm_u,
      egc_dm_r = dm_r_numerator / dm_r_denominator,
      egc_dd = egc_dd,
      mean_burden = mean_burden
    ),
    pairwise = do.call(rbind, pairwise)
  )
}
