script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)))) else normalizePath(getwd())
})
repo_root <- normalizePath(file.path(script_dir, ".."))
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else script_dir
dataset_dir <- file.path(output_dir, "datasets")
dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

options(digits = 17, scipen = 999)
RNGkind("L'Ecuyer-CMRG")
master_seed <- 20260730L
set.seed(master_seed)

streams <- vector("list", 20L)
streams[[1L]] <- .Random.seed
for (index in 2:20) {
  streams[[index]] <- parallel::nextRNGStream(streams[[index - 1L]])
}

with_stream <- function(index, expression) {
  assign(".Random.seed", streams[[index]], envir = .GlobalEnv)
  force(expression)
}

rank_map <- function(group_count) {
  setNames((seq_len(group_count) - 0.5) / group_count, paste0("G", seq_len(group_count)))
}

make_frame <- function(dataset_id, burden, group, ranks, weight = NULL) {
  if (is.null(weight)) weight <- rep(1, length(burden))
  stopifnot(length(burden) == length(group), length(weight) == length(burden))
  observed_rank <- unname(ranks[as.character(group)])
  stopifnot(all(is.finite(burden)), all(is.finite(weight)), all(weight > 0), all(is.finite(observed_rank)))
  data.frame(
    dataset_id = dataset_id,
    row_id = sprintf("%s_%05d", dataset_id, seq_along(burden)),
    burden = as.numeric(burden),
    group = as.character(group),
    rank = as.numeric(observed_rank),
    weight = as.numeric(weight),
    stringsAsFactors = FALSE
  )
}

five_groups <- function(counts = rep(60L, 5L)) {
  rep(paste0("G", 1:5), times = counts)
}

five_ranks <- rank_map(5L)
datasets <- vector("list", 20L)
descriptions <- character(20L)
scenarios <- character(20L)

# DS01: exact two-group manual example.
datasets[[1L]] <- with_stream(1L, {
  groups <- rep(c("G1", "G2"), each = 8L)
  make_frame("DS01", c(rep(1, 8L), rep(-1, 8L)), groups, c(G1 = 0.25, G2 = 0.75))
})
descriptions[[1L]] <- "two degenerate strata; exact manual example"
scenarios[[1L]] <- "manual_exact"

# DS02: identical empirical distributions in every stratum.
datasets[[2L]] <- with_stream(2L, {
  base <- seq(-2, 2, length.out = 21L)
  make_frame("DS02", rep(base, times = 5L), rep(paste0("G", 1:5), each = length(base)), five_ranks)
})
descriptions[[2L]] <- "five-stratum exact null"
scenarios[[2L]] <- "exact_null"

# DS03: deterministic linear location shift with common shape.
datasets[[3L]] <- with_stream(3L, {
  base <- seq(-1.5, 1.5, length.out = 31L)
  means <- 0.4 * (0.5 - unname(five_ranks))
  make_frame(
    "DS03",
    rep(base, times = 5L) + rep(means, each = length(base)),
    rep(paste0("G", 1:5), each = length(base)),
    five_ranks
  )
})
descriptions[[3L]] <- "deterministic linear location gradient"
scenarios[[3L]] <- "linear_exact"

# DS04: stochastic linear normal gradient.
datasets[[4L]] <- with_stream(4L, {
  groups <- five_groups()
  means <- 0.4 * (0.5 - unname(five_ranks))
  burden <- unlist(lapply(seq_len(5L), function(g) rnorm(60L, means[[g]], 1)), use.names = FALSE)
  make_frame("DS04", burden, groups, five_ranks)
})
descriptions[[4L]] <- "balanced stochastic linear normal gradient"
scenarios[[4L]] <- "linear_normal"

# DS05: inverse stochastic linear gradient.
datasets[[5L]] <- with_stream(5L, {
  groups <- five_groups()
  means <- -0.4 * (0.5 - unname(five_ranks))
  burden <- unlist(lapply(seq_len(5L), function(g) rnorm(60L, means[[g]], 1)), use.names = FALSE)
  make_frame("DS05", burden, groups, five_ranks)
})
descriptions[[5L]] <- "balanced inverse linear normal gradient"
scenarios[[5L]] <- "inverse_linear"

# DS06: symmetric U-shaped location pattern.
datasets[[6L]] <- with_stream(6L, {
  groups <- five_groups()
  raw <- (unname(five_ranks) - 0.5)^2
  means <- (raw - mean(raw)) * 0.4 / diff(range(raw))
  burden <- unlist(lapply(seq_len(5L), function(g) rnorm(60L, means[[g]], 1)), use.names = FALSE)
  make_frame("DS06", burden, groups, five_ranks)
})
descriptions[[6L]] <- "symmetric U-shaped location pattern"
scenarios[[6L]] <- "u_shape"

# DS07: inverted U-shaped location pattern.
datasets[[7L]] <- with_stream(7L, {
  groups <- five_groups()
  raw <- -(unname(five_ranks) - 0.5)^2
  means <- (raw - mean(raw)) * 0.4 / diff(range(raw))
  burden <- unlist(lapply(seq_len(5L), function(g) rnorm(60L, means[[g]], 1)), use.names = FALSE)
  make_frame("DS07", burden, groups, five_ranks)
})
descriptions[[7L]] <- "inverted U-shaped location pattern"
scenarios[[7L]] <- "inverted_u"

# DS08: scale gradient with exactly symmetric, zero-mean empirical distributions.
datasets[[8L]] <- with_stream(8L, {
  sigmas <- exp(0.5 * (0.5 - unname(five_ranks)))
  values <- unlist(lapply(seq_len(5L), function(g) {
    half <- abs(rnorm(30L)) * sigmas[[g]]
    c(-half, half)
  }), use.names = FALSE)
  make_frame("DS08", values, five_groups(), five_ranks)
})
descriptions[[8L]] <- "ordered scale gradient with exact zero means"
scenarios[[8L]] <- "scale_exact_mean_zero"

# DS09: skewed tails with sample means centered exactly at zero.
datasets[[9L]] <- with_stream(9L, {
  values <- unlist(lapply(seq_len(5L), function(g) {
    rate <- 0.7 + 0.25 * g
    x <- rexp(70L, rate = rate)
    if (g >= 4L) x <- -x
    x - mean(x)
  }), use.names = FALSE)
  make_frame("DS09", values, five_groups(rep(70L, 5L)), five_ranks)
})
descriptions[[9L]] <- "different skewed tails with exact zero sample means"
scenarios[[9L]] <- "tail_shape_mean_zero"

# DS10: symmetric crossing distributions with equal means and distinct spreads.
datasets[[10L]] <- with_stream(10L, {
  amplitudes <- c(2.0, 0.8, 1.4, 0.8, 2.0)
  values <- unlist(lapply(amplitudes, function(a) rep(c(-a, a), each = 30L)), use.names = FALSE)
  make_frame("DS10", values, five_groups(), five_ranks)
})
descriptions[[10L]] <- "symmetric crossing distributions"
scenarios[[10L]] <- "crossing_symmetric"

# DS11 and DS12 are deterministic transformations of DS03.
datasets[[11L]] <- with_stream(11L, {
  source <- datasets[[3L]]
  make_frame("DS11", source$burden - 5, source$group, setNames(source$rank[match(unique(source$group), source$group)], unique(source$group)), source$weight)
})
descriptions[[11L]] <- "negative translation of DS03 by five units"
scenarios[[11L]] <- "translation_check"

datasets[[12L]] <- with_stream(12L, {
  source <- datasets[[3L]]
  make_frame("DS12", source$burden * 0.08, source$group, setNames(source$rank[match(unique(source$group), source$group)], unique(source$group)), source$weight)
})
descriptions[[12L]] <- "positive scale transformation of DS03 by 0.08"
scenarios[[12L]] <- "scale_check"

# DS13: balanced strata with unequal individual weights.
datasets[[13L]] <- with_stream(13L, {
  groups <- five_groups()
  means <- 0.4 * (0.5 - unname(five_ranks))
  burden <- unlist(lapply(seq_len(5L), function(g) rnorm(60L, means[[g]], 1)), use.names = FALSE)
  weights <- exp(rnorm(length(groups), 0, 0.55))
  make_frame("DS13", burden, groups, five_ranks, weights)
})
descriptions[[13L]] <- "balanced strata with unequal positive weights"
scenarios[[13L]] <- "unequal_weights"

# DS14: moderately unbalanced strata (P2).
datasets[[14L]] <- with_stream(14L, {
  counts <- c(30L, 45L, 60L, 75L, 90L)
  groups <- five_groups(counts)
  means <- 0.4 * (0.5 - unname(five_ranks))
  burden <- unlist(lapply(seq_len(5L), function(g) rnorm(counts[[g]], means[[g]], 1)), use.names = FALSE)
  make_frame("DS14", burden, groups, five_ranks)
})
descriptions[[14L]] <- "moderately unbalanced five-stratum sample"
scenarios[[14L]] <- "unbalanced_p2"

# DS15: severely unbalanced strata (P3-like).
datasets[[15L]] <- with_stream(15L, {
  counts <- c(10L, 20L, 30L, 50L, 90L)
  groups <- five_groups(counts)
  means <- 0.4 * (0.5 - unname(five_ranks))
  burden <- unlist(lapply(seq_len(5L), function(g) rnorm(counts[[g]], means[[g]], 1)), use.names = FALSE)
  make_frame("DS15", burden, groups, five_ranks)
})
descriptions[[15L]] <- "severely unbalanced five-stratum sample"
scenarios[[15L]] <- "unbalanced_p3"

# DS16: discrete outcome with many ties.
datasets[[16L]] <- with_stream(16L, {
  probabilities <- list(
    c(0.05, 0.10, 0.20, 0.30, 0.35),
    c(0.10, 0.15, 0.25, 0.30, 0.20),
    c(0.15, 0.20, 0.30, 0.20, 0.15),
    c(0.20, 0.30, 0.25, 0.15, 0.10),
    c(0.35, 0.30, 0.20, 0.10, 0.05)
  )
  burden <- unlist(lapply(seq_len(5L), function(g) sample(-2:2, 80L, replace = TRUE, prob = probabilities[[g]])), use.names = FALSE)
  make_frame("DS16", burden, five_groups(rep(80L, 5L)), five_ranks)
})
descriptions[[16L]] <- "discrete burden with many ties"
scenarios[[16L]] <- "discrete_ties"

# DS17: heavy t3 residuals.
datasets[[17L]] <- with_stream(17L, {
  groups <- five_groups()
  means <- 0.4 * (0.5 - unname(five_ranks))
  burden <- unlist(lapply(seq_len(5L), function(g) means[[g]] + rt(60L, df = 3) / sqrt(3)), use.names = FALSE)
  make_frame("DS17", burden, groups, five_ranks)
})
descriptions[[17L]] <- "linear gradient with t3 residuals"
scenarios[[17L]] <- "heavy_tails"

# DS18: normal data with symmetric extreme contamination.
datasets[[18L]] <- with_stream(18L, {
  groups <- five_groups()
  means <- 0.4 * (0.5 - unname(five_ranks))
  blocks <- lapply(seq_len(5L), function(g) {
    x <- rnorm(60L, means[[g]], 1)
    x[1L] <- x[1L] + 7
    x[2L] <- x[2L] - 7
    x
  })
  make_frame("DS18", unlist(blocks, use.names = FALSE), groups, five_ranks)
})
descriptions[[18L]] <- "linear normal gradient with symmetric outliers"
scenarios[[18L]] <- "contamination"

# DS19: ten ordered strata.
datasets[[19L]] <- with_stream(19L, {
  ranks <- rank_map(10L)
  groups <- rep(paste0("G", 1:10), each = 25L)
  means <- 0.4 * (0.5 - unname(ranks))
  burden <- unlist(lapply(seq_len(10L), function(g) rnorm(25L, means[[g]], 1)), use.names = FALSE)
  make_frame("DS19", burden, groups, ranks)
})
descriptions[[19L]] <- "ten ordered strata with linear gradient"
scenarios[[19L]] <- "ten_strata"

# DS20: three strata with unequal sizes, weights, means, and scales.
datasets[[20L]] <- with_stream(20L, {
  counts <- c(40L, 70L, 110L)
  labels <- paste0("G", 1:3)
  probabilities <- counts / sum(counts)
  ranks <- setNames(c(probabilities[[1L]] / 2, probabilities[[1L]] + probabilities[[2L]] / 2, sum(probabilities[1:2]) + probabilities[[3L]] / 2), labels)
  groups <- rep(labels, times = counts)
  means <- c(0.25, 0.05, -0.10)
  sigmas <- c(1.25, 0.75, 1.05)
  burden <- unlist(lapply(seq_len(3L), function(g) rnorm(counts[[g]], means[[g]], sigmas[[g]])), use.names = FALSE)
  weights <- exp(rnorm(length(groups), 0, 0.35))
  make_frame("DS20", burden, groups, ranks, weights)
})
descriptions[[20L]] <- "three unequal strata with weights, means, and scales"
scenarios[[20L]] <- "three_strata_mixed"

manifest_rows <- vector("list", 20L)
for (index in seq_len(20L)) {
  dataset <- datasets[[index]]
  dataset_id <- sprintf("DS%02d", index)
  file_path <- file.path(dataset_dir, paste0(dataset_id, ".csv"))
  write.table(dataset, file_path, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE, na = "", qmethod = "double")
  manifest_rows[[index]] <- data.frame(
    dataset_id = dataset_id,
    description = descriptions[[index]],
    scenario = scenarios[[index]],
    row_count = nrow(dataset),
    group_count = length(unique(dataset$group)),
    total_weight = sum(dataset$weight),
    min_group_n = min(table(dataset$group)),
    max_group_n = max(table(dataset$group)),
    rng_seed = paste(streams[[index]], collapse = ";"),
    stringsAsFactors = FALSE
  )
}

combined <- do.call(rbind, datasets)
manifest <- do.call(rbind, manifest_rows)
write.table(combined, file.path(output_dir, "fixed_datasets_combined.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE, na = "", qmethod = "double")
write.table(manifest, file.path(output_dir, "dataset_manifest.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE, na = "", qmethod = "double")

environment_lines <- capture.output({
  cat("Master seed:", master_seed, "\n")
  cat("RNGkind:", paste(RNGkind(), collapse = " | "), "\n")
  cat("Dataset count:", length(datasets), "\n")
  cat("Total rows:", nrow(combined), "\n")
  print(sessionInfo())
})
writeLines(environment_lines, file.path(output_dir, "r_generation_environment.txt"), useBytes = TRUE)
cat(sprintf("Generated %d datasets and %d rows in %s\n", length(datasets), nrow(combined), output_dir))
