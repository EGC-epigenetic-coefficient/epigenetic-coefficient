script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)))) else normalizePath(getwd())
})
repo_root <- normalizePath(file.path(script_dir, ".."))
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else file.path(script_dir, "fixed_datasets_combined.csv")
output_dir <- if (length(args) >= 2) args[[2]] else file.path(script_dir, "results")
source(file.path(repo_root, "R/egc_reference.R"), local = TRUE)
options(digits = 17, scipen = 999)
RNGkind("L'Ecuyer-CMRG")
set.seed(20260730L)

data <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("dataset_id", "row_id", "burden", "group", "rank", "weight")
if (!all(required %in% names(data))) stop("missing required dataset columns")
dataset_ids <- unique(data$dataset_id)
component_rows <- list()
pairwise_rows <- list()
for (index in seq_along(dataset_ids)) {
  dataset_id <- dataset_ids[[index]]
  subset <- data[data$dataset_id == dataset_id, , drop = FALSE]
  rank_check <- aggregate(rank ~ group, subset, function(x) length(unique(x)))
  if (any(rank_check$rank != 1L)) stop(paste("inconsistent rank in", dataset_id))
  rank_values <- aggregate(rank ~ group, subset, function(x) unique(x)[[1L]])
  ranks <- setNames(rank_values$rank, rank_values$group)
  result <- estimate_egc_components(subset$burden, subset$group, ranks, subset$weight)
  component_rows[[index]] <- data.frame(
    dataset_id = dataset_id,
    egc_s = unname(result$components[["egc_s"]]),
    egc_lm = unname(result$components[["egc_lm"]]),
    egc_ld = unname(result$components[["egc_ld"]]),
    egc_dm_u = unname(result$components[["egc_dm_u"]]),
    egc_dm_r = unname(result$components[["egc_dm_r"]]),
    egc_dd = unname(result$components[["egc_dd"]]),
    mean_burden = unname(result$components[["mean_burden"]]),
    stringsAsFactors = FALSE
  )
  pairwise <- result$pairwise
  pairwise$dataset_id <- dataset_id
  pairwise <- pairwise[, c("dataset_id", "lower_group", "higher_group", "p_lower", "p_higher", "rank_distance", "w1", "delta")]
  pairwise_rows[[index]] <- pairwise
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(do.call(rbind, component_rows), file.path(output_dir, "r_components.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE, na = "", qmethod = "double")
write.table(do.call(rbind, pairwise_rows), file.path(output_dir, "r_pairwise.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE, na = "", qmethod = "double")
environment_lines <- capture.output({
  cat("Input:", input_file, "\n")
  cat("Datasets:", length(dataset_ids), "\n")
  cat("RNGkind:", paste(RNGkind(), collapse = " | "), "\n")
  print(sessionInfo())
})
writeLines(environment_lines, file.path(output_dir, "r_evaluation_environment.txt"), useBytes = TRUE)
cat(sprintf("Evaluated %d datasets in R\n", length(dataset_ids)))
