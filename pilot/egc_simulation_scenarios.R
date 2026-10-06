egc_screening_profiles <- list(
  P1 = c(0.20, 0.20, 0.20, 0.20, 0.20),
  P2 = c(0.10, 0.15, 0.20, 0.25, 0.30),
  P3 = c(0.05, 0.10, 0.15, 0.25, 0.45)
)

egc_profile_ridits <- function(probabilities) {
  cumulative <- c(0, cumsum(probabilities)[-length(probabilities)])
  cumulative + probabilities / 2
}

egc_profile <- function(profile_id) {
  if (!profile_id %in% names(egc_screening_profiles)) stop("unknown profile")
  probabilities <- egc_screening_profiles[[profile_id]]
  labels <- paste0("G", seq_along(probabilities))
  ranks <- egc_profile_ridits(probabilities)
  names(probabilities) <- labels
  names(ranks) <- labels
  list(probabilities = probabilities, ranks = ranks)
}

egc_allocate_counts <- function(n, probabilities) {
  expected <- n * probabilities
  counts <- floor(expected)
  remaining <- as.integer(n - sum(counts))
  if (remaining > 0L) {
    fractional <- expected - counts
    order_fractional <- order(-fractional, seq_along(fractional), method = "radix")
    counts[order_fractional[seq_len(remaining)]] <- counts[order_fractional[seq_len(remaining)]] + 1L
  }
  counts <- as.integer(counts)
  names(counts) <- names(probabilities)
  if (sum(counts) != n || any(counts <= 0L)) stop("invalid deterministic group allocation")
  counts
}

egc_weighted_center_and_range <- function(raw, probabilities, amplitude) {
  centered <- raw - sum(probabilities * raw)
  observed_range <- max(centered) - min(centered)
  if (!is.finite(observed_range) || observed_range <= 0) stop("shape has zero range")
  centered * amplitude / observed_range
}

egc_condition_value <- function(condition, name) {
  value <- condition[[name]]
  if (length(value) != 1L || is.na(value)) stop(paste("missing condition field", name))
  value
}

egc_scenario_signature <- function(condition) {
  scenario <- as.character(egc_condition_value(condition, "scenario"))
  profile_id <- as.character(egc_condition_value(condition, "profile"))
  profile <- egc_profile(profile_id)
  probabilities <- profile$probabilities
  ranks <- profile$ranks
  delta <- as.numeric(egc_condition_value(condition, "delta"))
  gamma <- as.numeric(egc_condition_value(condition, "gamma"))
  means <- setNames(rep(0, length(ranks)), names(ranks))
  sigmas <- setNames(rep(1, length(ranks)), names(ranks))
  tail_probabilities <- setNames(rep(0, length(ranks)), names(ranks))
  cross_strength <- setNames(rep(0, length(ranks)), names(ranks))
  distribution <- "normal"

  if (scenario == "C0") {
    # Complete equality.
  } else if (scenario == "C1") {
    means <- delta * (0.5 - ranks)
  } else if (scenario == "C2") {
    means <- -delta * (0.5 - ranks)
  } else if (scenario == "C3") {
    threshold_q <- as.numeric(egc_condition_value(condition, "threshold_q"))
    means <- egc_weighted_center_and_range(as.numeric(ranks <= threshold_q), probabilities, delta)
    names(means) <- names(ranks)
  } else if (scenario == "C4") {
    means <- egc_weighted_center_and_range((1 - ranks)^2, probabilities, delta)
    names(means) <- names(ranks)
  } else if (scenario == "C5") {
    shape <- as.character(egc_condition_value(condition, "shape"))
    raw <- (ranks - 0.5)^2
    if (shape == "INVERTED_U") raw <- -raw
    if (!shape %in% c("U", "INVERTED_U")) stop("invalid C5 shape")
    means <- egc_weighted_center_and_range(raw, probabilities, delta)
    names(means) <- names(ranks)
  } else if (scenario == "C6") {
    sigmas <- exp(gamma * (0.5 - ranks))
  } else if (scenario == "C7") {
    distribution <- "adverse_tail_mixture"
    tail_max <- as.numeric(egc_condition_value(condition, "tail_max"))
    tail_probabilities <- tail_max * (1 - ranks) / max(1 - ranks)
  } else if (scenario == "C8") {
    distribution <- "symmetric_crossing_mixture"
    cross_strength <- abs(ranks - 0.5) / max(abs(ranks - 0.5))
  } else if (scenario == "C9") {
    means <- delta * (0.5 - ranks)
    sigmas <- exp(gamma * (0.5 - ranks))
  } else {
    stop(paste("unsupported screening scenario", scenario))
  }

  list(
    scenario = scenario,
    profile = profile_id,
    probabilities = probabilities,
    ranks = ranks,
    means = means,
    sigmas = sigmas,
    distribution = distribution,
    tail_probabilities = tail_probabilities,
    tail_shift = if (scenario == "C7") 3 else 0,
    cross_strength = cross_strength,
    cross_a = if (scenario == "C8") as.numeric(egc_condition_value(condition, "cross_a")) else 0
  )
}

egc_generate_screening_dataset <- function(condition, n) {
  n <- as.integer(n)
  if (!is.finite(n) || n < 10L) stop("n must be at least 10")
  signature <- egc_scenario_signature(condition)
  counts <- egc_allocate_counts(n, signature$probabilities)
  labels <- names(counts)
  burden <- numeric(n)
  groups <- character(n)
  ranks <- numeric(n)
  cursor <- 1L

  for (label in labels) {
    size <- counts[[label]]
    idx <- seq.int(cursor, cursor + size - 1L)
    if (signature$distribution == "normal") {
      burden[idx] <- rnorm(size, mean = signature$means[[label]], sd = signature$sigmas[[label]])
    } else if (signature$distribution == "adverse_tail_mixture") {
      probability <- signature$tail_probabilities[[label]]
      adverse <- rbinom(size, size = 1L, prob = probability) == 1L
      main_mean <- -probability * signature$tail_shift / (1 - probability)
      observation_mean <- ifelse(adverse, signature$tail_shift, main_mean)
      burden[idx] <- rnorm(size, mean = observation_mean, sd = 1)
    } else if (signature$distribution == "symmetric_crossing_mixture") {
      strength <- signature$cross_strength[[label]]
      use_mixture <- rbinom(size, size = 1L, prob = strength) == 1L
      component_sign <- sample(c(-1, 1), size, replace = TRUE)
      residual_sd <- sqrt(max(0, 1 - signature$cross_a^2))
      mixture_value <- component_sign * signature$cross_a + rnorm(size, 0, residual_sd)
      burden[idx] <- ifelse(use_mixture, mixture_value, rnorm(size, 0, 1))
    } else {
      stop("unknown distribution")
    }
    groups[idx] <- label
    ranks[idx] <- signature$ranks[[label]]
    cursor <- cursor + size
  }

  data.frame(
    observation_id = seq_len(n),
    group = groups,
    rank = ranks,
    burden = burden,
    weight = rep(1, n),
    stringsAsFactors = FALSE
  )
}
