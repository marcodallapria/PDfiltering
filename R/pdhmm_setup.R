# Project-wide setup helpers for the reproducibility archive.
#
# This file provides:
# - package checks;
# - helper paths built with `here`;
# - utilities shared by the split manuscript scripts.

required_packages <- function() {
  c(
    "here", "DescTools", "partitions", "collections", "LaplacesDemon",
    "extraDistr", "scales", "vioplot", "pbapply", "ggplot2",
    "RColorBrewer", "patchwork", "RcppAlgos", "microbenchmark",
    "dplyr", "pbmcapply", "R.utils", "geomtextpath", "igraph",
    "lubridate", "knitr", "rmarkdown"
  )
}

check_required_packages <- function(pkgs = required_packages()) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      ". Please see REQUIREMENTS.md.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

project_path <- function(...) {
  here::here(...)
}

raw_data_path <- function(...) {
  project_path("data", "raw", ...)
}

cached_results_path <- function(...) {
  project_path("results", "cached", ...)
}

figure_path <- function(...) {
  project_path("results", "figures", ...)
}

table_path <- function(...) {
  project_path("results", "tables", ...)
}

log_path <- function(...) {
  project_path("results", "logs", ...)
}

ensure_output_directories <- function() {
  dirs <- c(
    project_path("results", "figures"),
    project_path("results", "tables"),
    project_path("results", "logs")
  )
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  invisible(dirs)
}

save_reference_session <- function() {
  ensure_output_directories()
  writeLines(capture.output(sessionInfo()), log_path("sessionInfo_runtime.txt"))
}

safe_notify <- function(message_text = "Computation finished") {
  if (Sys.info()[["sysname"]] == "Darwin" && nzchar(Sys.which("say"))) {
    system2("say", message_text)
  }
  invisible(NULL)
}
