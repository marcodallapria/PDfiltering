# Quick reproduction driver.
#
# This script is intended for referee-style reproduction runs. It uses the
# split analysis scripts and prefers cached computational outputs where the
# manuscript workflow is especially expensive.

source(here::here("analysis", "00_setup.R"))

message("Running quick reproduction workflow")
source(here::here("analysis", "01_infectious_analysis.R"), local = FALSE)
source(here::here("analysis", "02_simulation_study.R"), local = FALSE)
source(here::here("analysis", "03_benchmarks_and_supplement.R"), local = FALSE)
save_reference_session()
message("Quick reproduction workflow completed")
