args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(values) {
  result <- list()
  for (value in values) {
    if (!startsWith(value, "--") || !grepl("=", value, fixed = TRUE)) {
      stop(paste("arguments must use --name=value:", value))
    }
    parts <- strsplit(sub("^--", "", value), "=", fixed = TRUE)[[1L]]
    result[[parts[[1L]]]] <- paste(parts[-1L], collapse = "=")
  }
  result
}

options(stringsAsFactors = FALSE, digits = 17)
parsed <- parse_args(args)
`%||%` <- function(x, y) if (is.null(x)) y else x
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)))) else normalizePath(getwd())
})
repo_root <- normalizePath(file.path(script_dir, ".."))
registry_path <- parsed$registry %||% stop("--registry is required")
output_dir <- parsed$`output-dir` %||% stop("--output-dir is required")
condition_ids <- strsplit(parsed$`condition-ids` %||% stop("--condition-ids is required"), ",", fixed = TRUE)[[1L]]
sample_sizes <- as.integer(strsplit(parsed$`sample-sizes` %||% "250,1000", ",", fixed = TRUE)[[1L]])
replications <- as.integer(parsed$replications %||% "1")
B_permutation <- as.integer(parsed$`B-permutation` %||% "499")
B_bootstrap <- as.integer(parsed$`B-bootstrap` %||% "499")
gammas <- as.numeric(strsplit(parsed$gammas %||% "0.6,0.7,0.8", ",", fixed = TRUE)[[1L]])
K_extrapolation <- as.integer(parsed$`K-extrapolation` %||% "199")
master_seed <- as.integer(parsed$`master-seed` %||% "20260731")
cores <- as.integer(parsed$cores %||% "1")

source(file.path(repo_root, "R/egc_reference.R"))
source(file.path(repo_root, "R/egc_inference.R"))
source(file.path(script_dir, "egc_simulation_scenarios.R"))
source(file.path(script_dir, "egc_simulation_engine.R"))

if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("append-only contract: profiled output directory is not empty")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  line <- system2("shasum", c("-a", "256", shQuote(normalizePath(path))), stdout = TRUE)
  strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

write_deterministic_csv <- function(data, path) {
  write.table(
    data, path, sep = ",", row.names = FALSE, col.names = TRUE,
    quote = TRUE, na = "", qmethod = "double", fileEncoding = "UTF-8"
  )
}

registry <- read.csv(registry_path, check.names = FALSE, stringsAsFactors = FALSE)
jobs <- egc_make_jobs(registry, condition_ids, sample_sizes, replications)
streams <- egc_generate_selected_job_streams(jobs$stream_index, master_seed)

.benchmark_depth <- 0L
.benchmark_job <- NULL
.benchmark_timings <- list()
.benchmark_timing_position <- 1L

record_timed_call <- function(module, original, arguments) {
  if (.benchmark_depth > 0L) return(do.call(original, arguments))
  .benchmark_depth <<- .benchmark_depth + 1L
  before <- proc.time()
  succeeded <- FALSE
  message_text <- ""
  on.exit({
    elapsed <- proc.time() - before
    .benchmark_timings[[.benchmark_timing_position]] <<- data.frame(
      job_id = as.character(.benchmark_job$job_id),
      condition_id = as.character(.benchmark_job$condition_id),
      scenario = as.character(.benchmark_job$scenario),
      profile = as.character(.benchmark_job$profile),
      N = as.integer(.benchmark_job$N),
      replication = as.integer(.benchmark_job$replication),
      module = as.character(module),
      user_seconds = unname(elapsed[["user.self"]]),
      system_seconds = unname(elapsed[["sys.self"]]),
      elapsed_seconds = unname(elapsed[["elapsed"]]),
      status = if (succeeded) "PASS" else "ERROR",
      message = message_text,
      benchmark_excluded = TRUE,
      stringsAsFactors = FALSE
    )
    .benchmark_timing_position <<- .benchmark_timing_position + 1L
    .benchmark_depth <<- .benchmark_depth - 1L
  }, add = TRUE)
  value <- tryCatch(
    do.call(original, arguments),
    error = function(error) {
      message_text <<- conditionMessage(error)
      stop(error)
    }
  )
  succeeded <- TRUE
  value
}

wrap_simple <- function(name, module) {
  original <- get(name, envir = .GlobalEnv)
  assign(name, function(...) record_timed_call(module, original, list(...)), envir = .GlobalEnv)
}

wrap_simple("egc_generate_screening_dataset", "DATA_GENERATION")
wrap_simple("estimate_egc_inference_statistics", "POINT_PRIMARY")
wrap_simple("egc_secondary_point_statistics", "POINT_SECONDARY")
wrap_simple("egc_permutation_test", "PERMUTATION")
wrap_simple("egc_stratified_bootstrap", "BOOTSTRAP_N")
wrap_simple("egc_subsampling_extrapolation", "BC_E_EXTRAPOLATION")

original_rescaled <- egc_rescaled_interval_multi
egc_rescaled_interval_multi <- function(...) {
  arguments <- list(...)
  gamma <- if (!is.null(arguments$gamma)) arguments$gamma else arguments[[6L]]
  replace <- if (!is.null(arguments$replace)) arguments$replace else arguments[[7L]]
  module <- sprintf("CI_%s_GAMMA_%0.2f", if (isTRUE(replace)) "M" else "S", as.numeric(gamma))
  record_timed_call(module, original_rescaled, arguments)
}

started <- proc.time()[["elapsed"]]

worker <- function(index) {
  .benchmark_depth <<- 0L
  .benchmark_timings <<- list()
  .benchmark_timing_position <<- 1L
  .benchmark_job <- jobs[index, , drop = FALSE]
  assign(".benchmark_job", .benchmark_job, envir = .GlobalEnv)
  invisible(gc(reset = TRUE))
  result <- tryCatch(
    egc_run_screening_job(
      .benchmark_job, streams[[index]],
      B_permutation = B_permutation,
      B_bootstrap = B_bootstrap,
      gammas = gammas,
      K_extrapolation = K_extrapolation,
      alpha = 0.05,
      dry_run_excluded = TRUE
    ),
    error = function(error) list(
      replication = NULL,
      methods = NULL,
      error = data.frame(
        job_id = as.character(.benchmark_job$job_id),
        condition_id = as.character(.benchmark_job$condition_id),
        scenario = as.character(.benchmark_job$scenario),
        N = as.integer(.benchmark_job$N),
        replication = as.integer(.benchmark_job$replication),
        error_type = class(error)[[1L]],
        message = conditionMessage(error),
        stringsAsFactors = FALSE
      )
    )
  )
  memory <- gc()
  memory_row <- data.frame(
    job_id = as.character(.benchmark_job$job_id),
    condition_id = as.character(.benchmark_job$condition_id),
    scenario = as.character(.benchmark_job$scenario),
    profile = as.character(.benchmark_job$profile),
    N = as.integer(.benchmark_job$N),
    replication = as.integer(.benchmark_job$replication),
    Ncells_max_mb = as.numeric(memory["Ncells", 7L]),
    Vcells_max_mb = as.numeric(memory["Vcells", 7L]),
    R_heap_max_mb_sum = as.numeric(memory["Ncells", 7L] + memory["Vcells", 7L]),
    benchmark_excluded = TRUE,
    stringsAsFactors = FALSE
  )
  list(result = result, timings = .benchmark_timings, memory = memory_row)
}

indices <- seq_len(nrow(jobs))
if (cores > 1L && .Platform$OS.type != "windows") {
  worker_results <- parallel::mclapply(
    indices, worker, mc.cores = cores, mc.set.seed = FALSE, mc.preschedule = FALSE
  )
} else {
  worker_results <- lapply(indices, worker)
}

replication_rows <- Filter(Negate(is.null), lapply(worker_results, function(item) item$result$replication))
method_rows <- Filter(Negate(is.null), lapply(worker_results, function(item) item$result$methods))
error_rows <- Filter(Negate(is.null), lapply(worker_results, function(item) item$result$error))
timing_groups <- lapply(worker_results, function(item) item$timings)
memory_rows <- lapply(worker_results, function(item) item$memory)

replication <- if (length(replication_rows)) do.call(rbind, replication_rows) else data.frame()
methods <- if (length(method_rows)) do.call(rbind, method_rows) else data.frame()
errors <- if (length(error_rows)) do.call(rbind, error_rows) else data.frame(
  job_id = character(), condition_id = character(), scenario = character(),
  N = integer(), replication = integer(), error_type = character(), message = character(),
  stringsAsFactors = FALSE
)

if (nrow(replication)) replication <- replication[order(replication$job_id, method = "radix"), , drop = FALSE]
if (nrow(methods)) methods <- methods[order(methods$job_id, methods$method, methods$estimand, method = "radix"), , drop = FALSE]
if (nrow(errors)) errors <- errors[order(errors$job_id, method = "radix"), , drop = FALSE]

replication_path <- file.path(output_dir, "replication_results.csv")
method_path <- file.path(output_dir, "method_results.csv")
error_path <- file.path(output_dir, "errors.csv")
parameter_path <- file.path(output_dir, "run_parameters.csv")
write_deterministic_csv(replication, replication_path)
write_deterministic_csv(methods, method_path)
write_deterministic_csv(errors, error_path)

parameters <- data.frame(
  key = c(
    "mode", "scientific_results", "master_seed", "rng_kind", "condition_ids",
    "sample_sizes", "replications", "B_permutation", "B_bootstrap", "gammas",
    "K_extrapolation", "jobs_total", "job_start", "job_end", "jobs_expected", "dry_run_excluded"
  ),
  value = c(
    "dry_run", "false", master_seed, "L'Ecuyer-CMRG", paste(sort(condition_ids), collapse = ";"),
    paste(sort(sample_sizes), collapse = ";"), replications, B_permutation, B_bootstrap,
    paste(gammas, collapse = ";"), K_extrapolation, nrow(jobs), 1L, nrow(jobs),
    nrow(jobs), "true"
  ),
  stringsAsFactors = FALSE
)
write_deterministic_csv(parameters, parameter_path)

manifest <- data.frame(
  file = c("replication_results.csv", "method_results.csv", "errors.csv", "run_parameters.csv"),
  sha256 = vapply(c(replication_path, method_path, error_path, parameter_path), sha256, character(1)),
  stringsAsFactors = FALSE
)
manifest$jobs_expected <- nrow(jobs)
manifest$jobs_completed <- nrow(replication)
manifest$errors <- nrow(errors)
manifest_path <- file.path(output_dir, "deterministic_manifest.csv")
write_deterministic_csv(manifest, manifest_path)

timing_rows <- unlist(timing_groups, recursive = FALSE)
timings <- if (length(timing_rows)) do.call(rbind, timing_rows) else data.frame()
if (nrow(timings)) timings <- timings[order(timings$job_id, timings$module, method = "radix"), , drop = FALSE]
write.csv(timings, file.path(output_dir, "module_timings.csv"), row.names = FALSE, na = "")

memory_table <- do.call(rbind, memory_rows)
memory_table <- memory_table[order(memory_table$job_id, method = "radix"), , drop = FALSE]
write.csv(memory_table, file.path(output_dir, "job_memory.csv"), row.names = FALSE, na = "")

elapsed <- proc.time()[["elapsed"]] - started
runtime <- data.frame(
  key = c("cores", "elapsed_seconds", "platform", "R_version", "note"),
  value = c(
    cores, format(elapsed, digits = 12), R.version$platform, R.version.string,
    "profiled excluded benchmark; timing and memory files excluded from deterministic comparison"
  ),
  stringsAsFactors = FALSE
)
write.csv(runtime, file.path(output_dir, "runtime_metadata.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

cat(sprintf(
  "COMPLETE benchmark_excluded jobs=%d completed=%d errors=%d timing_rows=%d memory_rows=%d cores=%d elapsed=%.3f seconds\n",
  nrow(jobs), nrow(replication), nrow(errors), nrow(timings), nrow(memory_table), cores, elapsed
))
