if (!exists("egc_generate_screening_dataset", mode = "function")) {
  stop("source egc_screening_scenarios_R.R before egc_screening_engine_R.R")
}
if (!exists("estimate_egc_inference_statistics", mode = "function")) {
  stop("source egc_reference_R.R and egc_inference_R.R before egc_screening_engine_R.R")
}

egc_set_stream_seed <- function(seed) {
  assign(".Random.seed", as.integer(seed), envir = .GlobalEnv)
  invisible(NULL)
}

egc_generate_job_streams <- function(number, master_seed) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(master_seed))
  seed <- .Random.seed
  streams <- vector("list", number)
  if (number == 0L) return(streams)
  for (index in seq_len(number)) {
    streams[[index]] <- seed
    seed <- parallel::nextRNGStream(seed)
  }
  streams
}

egc_generate_selected_job_streams <- function(stream_indices, master_seed) {
  stream_indices <- as.integer(stream_indices)
  if (length(stream_indices) == 0L || any(stream_indices < 1L) || anyDuplicated(stream_indices)) {
    stop("stream indices must be unique positive integers")
  }
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(master_seed))
  seed <- .Random.seed
  wanted <- setNames(seq_along(stream_indices), as.character(stream_indices))
  streams <- vector("list", length(stream_indices))
  for (index in seq_len(max(stream_indices))) {
    key <- as.character(index)
    if (key %in% names(wanted)) streams[[wanted[[key]]]] <- seed
    seed <- parallel::nextRNGStream(seed)
  }
  streams
}

egc_module_streams <- function(job_seed, number) {
  streams <- vector("list", number)
  seed <- job_seed
  for (index in seq_len(number)) {
    streams[[index]] <- seed
    seed <- parallel::nextRNGSubStream(seed)
  }
  streams
}

egc_make_jobs <- function(registry, condition_ids, sample_sizes, replications) {
  selected <- registry[registry$condition_id %in% condition_ids, , drop = FALSE]
  missing <- setdiff(condition_ids, selected$condition_id)
  if (length(missing) > 0L) stop(paste("unknown condition ids", paste(missing, collapse = ",")))
  selected <- selected[order(selected$condition_id, method = "radix"), , drop = FALSE]
  sample_sizes <- sort(unique(as.integer(sample_sizes)))
  rows <- list()
  position <- 1L
  for (condition_index in seq_len(nrow(selected))) {
    condition <- selected[condition_index, , drop = FALSE]
    for (n in sample_sizes) {
      for (replication in seq_len(as.integer(replications))) {
        row <- condition
        row$N <- as.integer(n)
        row$replication <- as.integer(replication)
        row$job_id <- sprintf("%s__N%05d__R%06d", row$condition_id, n, replication)
        rows[[position]] <- row
        position <- position + 1L
      }
    }
  }
  jobs <- do.call(rbind, rows)
  jobs <- jobs[order(jobs$condition_id, jobs$N, jobs$replication, method = "radix"), , drop = FALSE]
  jobs$stream_index <- seq_len(nrow(jobs))
  rownames(jobs) <- NULL
  jobs
}

egc_secondary_point_statistics <- function(data, group_ranks) {
  individual_ranks <- unname(group_ranks[as.character(data$group)])
  beta <- weighted_covariance(individual_ranks, data$burden, data$weight) /
    weighted_covariance(individual_ranks, individual_ranks, data$weight)
  labels <- names(sort(group_ranks))
  probabilities <- setNames(numeric(length(labels)), labels)
  means <- probabilities
  for (label in labels) {
    mask <- data$group == label
    probabilities[[label]] <- sum(data$weight[mask]) / sum(data$weight)
    means[[label]] <- weighted_mean(data$burden[mask], data$weight[mask])
  }
  overall <- weighted_mean(data$burden, data$weight)
  c(
    egc_s = -beta,
    egc_lm = sqrt(sum(probabilities * (means - overall)^2)),
    egc_ld = -weighted_correlation(unname(group_ranks[labels]), unname(means), unname(probabilities)),
    mean_burden = overall
  )
}

egc_rescaled_interval_multi <- function(
    burden, groups, group_ranks, weights, B, gamma, replace, minimum_per_group = 10L) {
  group_counts <- as.numeric(table(groups))
  if (any(group_counts < minimum_per_group)) {
    return(list(
      applicable = FALSE,
      ci_dm_u = c(lower = NA_real_, upper = NA_real_),
      ci_A = c(lower = NA_real_, upper = NA_real_),
      failure_fraction = 1,
      sign_stability_A = NA_real_,
      realized_m = NA_integer_
    ))
  }
  n <- length(burden)
  nominal_m <- floor(n^gamma)
  fraction <- nominal_m / n
  realized_group_m <- pmax(minimum_per_group, floor(group_counts * fraction))
  if (!replace) realized_group_m <- pmin(realized_group_m, group_counts)
  realized_m <- sum(realized_group_m)
  result <- egc_stratified_bootstrap(
    burden, groups, group_ranks, weights, B = B, fraction = fraction,
    replace = replace, minimum_per_group = minimum_per_group
  )
  valid <- is.finite(result$bootstrap[, "egc_dm_u"]) &
    is.finite(result$bootstrap[, "orientation_A"])
  if (sum(valid) < 2L) {
    return(list(
      applicable = TRUE,
      ci_dm_u = c(lower = NA_real_, upper = NA_real_),
      ci_A = c(lower = NA_real_, upper = NA_real_),
      failure_fraction = 1 - mean(valid),
      sign_stability_A = NA_real_,
      realized_m = realized_m
    ))
  }
  interval_for <- function(name, clamp_zero = FALSE) {
    observed <- result$observed[[name]]
    roots <- sqrt(realized_m) * (result$bootstrap[valid, name] - observed)
    q <- unname(quantile(roots, c(0.025, 0.975), na.rm = TRUE))
    interval <- c(
      lower = observed - q[[2L]] / sqrt(n),
      upper = observed - q[[1L]] / sqrt(n)
    )
    if (clamp_zero) interval <- pmax(0, interval)
    interval
  }
  sign_A <- sign(result$observed[["orientation_A"]])
  list(
    applicable = TRUE,
    ci_dm_u = interval_for("egc_dm_u", clamp_zero = TRUE),
    ci_A = interval_for("orientation_A", clamp_zero = FALSE),
    failure_fraction = 1 - mean(valid),
    sign_stability_A = if (sign_A == 0) NA_real_ else mean(sign(result$bootstrap[valid, "orientation_A"]) == sign_A),
    realized_m = realized_m
  )
}

egc_method_row <- function(
    job, method, estimand, estimate = NA_real_, ci = c(NA_real_, NA_real_),
    p_value = NA_real_, B = NA_integer_, failure_fraction = 0,
    admissible = TRUE, notes = "") {
  data.frame(
    job_id = as.character(job$job_id),
    condition_id = as.character(job$condition_id),
    scenario = as.character(job$scenario),
    N = as.integer(job$N),
    replication = as.integer(job$replication),
    method = as.character(method),
    estimand = as.character(estimand),
    estimate = as.numeric(estimate),
    ci_lower = as.numeric(ci[[1L]]),
    ci_upper = as.numeric(ci[[2L]]),
    p_value = as.numeric(p_value),
    B = as.integer(B),
    failure_fraction = as.numeric(failure_fraction),
    admissible = isTRUE(admissible),
    notes = as.character(notes),
    stringsAsFactors = FALSE
  )
}

egc_run_screening_job <- function(
    job, job_seed, B_permutation, B_bootstrap, gammas, K_extrapolation,
    alpha = 0.05, minimum_per_group = 10L, dry_run_excluded = TRUE) {
  module_count <- 4L + 2L * length(gammas)
  module_seeds <- egc_module_streams(job_seed, module_count)
  seed_position <- 1L

  egc_set_stream_seed(module_seeds[[seed_position]])
  seed_position <- seed_position + 1L
  data <- egc_generate_screening_dataset(job, as.integer(job$N))
  ranks <- tapply(data$rank, data$group, unique)
  point <- estimate_egc_inference_statistics(data$burden, data$group, ranks, data$weight)
  secondary <- egc_secondary_point_statistics(data, ranks)

  egc_set_stream_seed(module_seeds[[seed_position]])
  seed_position <- seed_position + 1L
  permutation <- egc_permutation_test(
    data$burden, data$group, ranks, data$weight, B = as.integer(B_permutation)
  )

  egc_set_stream_seed(module_seeds[[seed_position]])
  seed_position <- seed_position + 1L
  bootstrap_n <- egc_stratified_bootstrap(
    data$burden, data$group, ranks, data$weight, B = as.integer(B_bootstrap)
  )
  bias <- egc_bias_candidates(permutation, bootstrap_n)

  sensitivity <- list()
  sensitivity_rows <- list()
  sensitivity_position <- 1L
  for (gamma in gammas) {
    for (replace in c(TRUE, FALSE)) {
      egc_set_stream_seed(module_seeds[[seed_position]])
      seed_position <- seed_position + 1L
      result <- egc_rescaled_interval_multi(
        data$burden, data$group, ranks, data$weight,
        B = as.integer(B_bootstrap), gamma = gamma, replace = replace,
        minimum_per_group = minimum_per_group
      )
      label <- sprintf("CI-%s-gamma%.2f", if (replace) "M" else "S", gamma)
      sensitivity[[label]] <- result
      sensitivity_rows[[sensitivity_position]] <- egc_method_row(
        job, label, "egc_dm_u", point$statistics[["egc_dm_u"]], result$ci_dm_u,
        B = B_bootstrap, failure_fraction = result$failure_fraction,
        admissible = result$applicable,
        notes = if (result$applicable) sprintf("realized_m=%s", result$realized_m) else "minimum group size not reached"
      )
      sensitivity_position <- sensitivity_position + 1L
      sensitivity_rows[[sensitivity_position]] <- egc_method_row(
        job, paste0(label, "-A"), "orientation_A", point$statistics[["orientation_A"]], result$ci_A,
        B = B_bootstrap, failure_fraction = result$failure_fraction,
        admissible = result$applicable,
        notes = if (result$applicable) sprintf("realized_m=%s", result$realized_m) else "minimum group size not reached"
      )
      sensitivity_position <- sensitivity_position + 1L
    }
  }

  egc_set_stream_seed(module_seeds[[seed_position]])
  extrapolation <- egc_subsampling_extrapolation(
    data$burden, data$group, ranks, data$weight,
    fractions = c(0.5, 0.7, 0.9), K = as.integer(K_extrapolation)
  )

  observed_A <- point$statistics[["orientation_A"]]
  ci_A_n <- bootstrap_n$percentile_ci[, "orientation_A"]
  ci_A_excludes_zero <- all(is.finite(ci_A_n)) && (ci_A_n[[1L]] > 0 || ci_A_n[[2L]] < 0)
  sensitivity_midpoint_signs <- vapply(
    sensitivity,
    function(item) {
      if (!item$applicable || any(!is.finite(item$ci_A))) return(NA_real_)
      sign(mean(item$ci_A))
    },
    numeric(1)
  )
  sensitivity_sign_ok <- all(is.na(sensitivity_midpoint_signs) | sensitivity_midpoint_signs == sign(observed_A))
  group_min_n <- min(point$group_n)
  group_min_neff <- min(point$group_neff)
  support_critical <- group_min_n < 20 || group_min_neff < 20
  orientation_gate <- isTRUE(
    permutation$p_dm <= alpha &&
      permutation$p_D <= alpha &&
      ci_A_excludes_zero &&
      is.finite(bootstrap_n$sign_stability_A) && bootstrap_n$sign_stability_A >= 0.95 &&
      is.finite(point$statistics[["directional_coherence"]]) && point$statistics[["directional_coherence"]] >= 0.50 &&
      sensitivity_sign_ok &&
      !support_critical &&
      bootstrap_n$failure_fraction <= 0.01
  )

  replication_row <- data.frame(
    job_id = as.character(job$job_id),
    condition_id = as.character(job$condition_id),
    scenario = as.character(job$scenario),
    block = as.character(job$block),
    profile = as.character(job$profile),
    variant = as.character(job$variant),
    effect_level = as.character(job$effect_level),
    N = as.integer(job$N),
    replication = as.integer(job$replication),
    stream_index = as.integer(job$stream_index),
    stream_signature = paste(job_seed[1:7], collapse = ":"),
    egc_s = secondary[["egc_s"]],
    egc_lm = secondary[["egc_lm"]],
    egc_ld = secondary[["egc_ld"]],
    egc_dm_u = point$statistics[["egc_dm_u"]],
    egc_dm_r = point$statistics[["egc_dm_r"]],
    orientation_A = observed_A,
    orientation_D = point$statistics[["orientation_D"]],
    egc_dd = point$statistics[["egc_dd"]],
    directional_coherence = point$statistics[["directional_coherence"]],
    mean_burden = secondary[["mean_burden"]],
    p_dm = permutation$p_dm,
    p_D = permutation$p_D,
    ci_n_dm_lower = bootstrap_n$percentile_ci[["2.5%", "egc_dm_u"]],
    ci_n_dm_upper = bootstrap_n$percentile_ci[["97.5%", "egc_dm_u"]],
    ci_n_A_lower = ci_A_n[[1L]],
    ci_n_A_upper = ci_A_n[[2L]],
    sign_stability_A = bootstrap_n$sign_stability_A,
    bootstrap_failure_fraction = bootstrap_n$failure_fraction,
    group_min_n = group_min_n,
    group_min_neff = group_min_neff,
    support_critical = support_critical,
    orientation_gate = orientation_gate,
    bias_raw = bias[["raw"]],
    bias_perm_mean = bias[["permutation_mean_subtraction"]],
    bias_perm_median = bias[["permutation_median_subtraction"]],
    bias_bootstrap = bias[["within_stratum_bootstrap_bias_correction"]],
    bias_quadrature = bias[["quadrature_null_correction"]],
    bias_extrapolation = extrapolation$estimate,
    extrapolation_r_squared = extrapolation$r_squared,
    extrapolation_admissible = extrapolation$admissible,
    dry_run_excluded = isTRUE(dry_run_excluded),
    stringsAsFactors = FALSE
  )

  methods <- list(
    egc_method_row(job, "POINT", "egc_dm_u", point$statistics[["egc_dm_u"]]),
    egc_method_row(job, "POINT", "egc_dm_r", point$statistics[["egc_dm_r"]]),
    egc_method_row(job, "POINT", "orientation_A", observed_A),
    egc_method_row(job, "POINT", "orientation_D", point$statistics[["orientation_D"]]),
    egc_method_row(job, "POINT", "egc_dd", point$statistics[["egc_dd"]]),
    egc_method_row(job, "PERM-DM", "egc_dm_u", point$statistics[["egc_dm_u"]], p_value = permutation$p_dm, B = B_permutation),
    egc_method_row(job, "PERM-D", "orientation_D", point$statistics[["orientation_D"]], p_value = permutation$p_D, B = B_permutation),
    egc_method_row(job, "CI-N", "egc_dm_u", point$statistics[["egc_dm_u"]], bootstrap_n$percentile_ci[, "egc_dm_u"], B = B_bootstrap, failure_fraction = bootstrap_n$failure_fraction),
    egc_method_row(job, "CI-N-A", "orientation_A", observed_A, ci_A_n, B = B_bootstrap, failure_fraction = bootstrap_n$failure_fraction),
    egc_method_row(job, "BC-PM", "egc_dm_u", bias[["permutation_mean_subtraction"]], B = B_permutation, notes = "negative control"),
    egc_method_row(job, "BC-PD", "egc_dm_u", bias[["permutation_median_subtraction"]], B = B_permutation, notes = "negative control"),
    egc_method_row(job, "BC-B", "egc_dm_u", bias[["within_stratum_bootstrap_bias_correction"]], B = B_bootstrap),
    egc_method_row(job, "BC-Q", "egc_dm_u", bias[["quadrature_null_correction"]], B = B_permutation),
    egc_method_row(job, "BC-E", "egc_dm_u", extrapolation$estimate, B = K_extrapolation, admissible = extrapolation$admissible, notes = sprintf("R2=%.12g", extrapolation$r_squared))
  )
  methods <- c(methods, sensitivity_rows)
  list(replication = replication_row, methods = do.call(rbind, methods), error = NULL)
}

egc_run_screening_jobs <- function(
    jobs, streams, B_permutation, B_bootstrap, gammas, K_extrapolation,
    alpha = 0.05, cores = 1L, dry_run_excluded = TRUE) {
  if (length(streams) != nrow(jobs)) stop("stream count differs from job count")
  worker <- function(index) {
    tryCatch(
      egc_run_screening_job(
        jobs[index, , drop = FALSE], streams[[index]],
        B_permutation = B_permutation, B_bootstrap = B_bootstrap,
        gammas = gammas, K_extrapolation = K_extrapolation, alpha = alpha,
        dry_run_excluded = dry_run_excluded
      ),
      error = function(error) list(
        replication = NULL,
        methods = NULL,
        error = data.frame(
          job_id = as.character(jobs$job_id[[index]]),
          condition_id = as.character(jobs$condition_id[[index]]),
          scenario = as.character(jobs$scenario[[index]]),
          N = as.integer(jobs$N[[index]]),
          replication = as.integer(jobs$replication[[index]]),
          error_type = class(error)[[1L]],
          message = conditionMessage(error),
          stringsAsFactors = FALSE
        )
      )
    )
  }
  indices <- seq_len(nrow(jobs))
  if (as.integer(cores) > 1L && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(indices, worker, mc.cores = as.integer(cores), mc.set.seed = FALSE, mc.preschedule = FALSE)
  } else {
    results <- lapply(indices, worker)
  }
  replication_rows <- Filter(Negate(is.null), lapply(results, `[[`, "replication"))
  method_rows <- Filter(Negate(is.null), lapply(results, `[[`, "methods"))
  error_rows <- Filter(Negate(is.null), lapply(results, `[[`, "error"))
  list(
    replication = if (length(replication_rows)) do.call(rbind, replication_rows) else data.frame(),
    methods = if (length(method_rows)) do.call(rbind, method_rows) else data.frame(),
    errors = if (length(error_rows)) do.call(rbind, error_rows) else data.frame()
  )
}
