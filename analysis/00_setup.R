# Shared setup for manuscript-specific scripts.
#
# All split analysis scripts should source this file first. It ensures that
# relative paths are resolved with `here`, checks package availability, and
# loads the reusable PD-HMM methods.

source(here::here("R", "pdhmm_setup.R"))
here::i_am("analysis/00_setup.R")
check_required_packages()
ensure_output_directories()
source(here::here("R", "pdhmm_reusable.R"))
