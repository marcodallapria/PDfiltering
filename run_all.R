# Fuller reproduction driver.
#
# This script renders the section-level notebooks. The notebooks are
# self-contained and therefore serve as the main full-workflow documents.

here::i_am("run_all.R")
source(here::here("R", "pdhmm_setup.R"))
check_required_packages(c("here", "knitr", "rmarkdown"))
ensure_output_directories()

if (!rmarkdown::pandoc_available()) {
  stop("Pandoc is required to render the notebooks. Please use RStudio or install Pandoc.", call. = FALSE)
}

rmarkdown::render(here::here("notebooks", "infectious_analysis.Rmd"))
rmarkdown::render(here::here("notebooks", "simulation_study.Rmd"))
rmarkdown::render(here::here("notebooks", "benchmarks.Rmd"))

save_reference_session()
message("Full reproduction workflow completed")
