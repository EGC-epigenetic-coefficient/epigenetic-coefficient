command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
if (length(file_arg) == 1L) {
  script_path <- sub("^--file=", "", file_arg)
  # Rscript encodes spaces as "~+~" in --file on some platforms.
  script_path <- gsub("~\\+~", " ", script_path)
  root <- dirname(normalizePath(script_path))
} else {
  root <- normalizePath(getwd())
}

source(file.path(root, "..", "..", "R", "egc_reference.R"))

run_example <- function(example, input_file) {
  data <- read.csv(file.path(root, input_file), stringsAsFactors = FALSE)
  ranks <- tapply(data$ses_rank, data$ses_group, unique)
  result <- estimate_egc_components(data$z_burden, data$ses_group, ranks, data$weight)$components
  data.frame(
    example = example,
    metric = names(result),
    value = as.numeric(result),
    implementation = "R",
    stringsAsFactors = FALSE
  )
}

output <- rbind(
  run_example("A", "minimal_example_input.csv"),
  run_example("B", "minimal_example_B_input.csv")
)
write.csv(output, file.path(root, "minimal_example_results_R.csv"), row.names = FALSE, quote = TRUE)
print(output, row.names = FALSE)
