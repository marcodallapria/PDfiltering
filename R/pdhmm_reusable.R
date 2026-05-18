# Reusable software entry point.
#
# The documented public API lives in `R/pdhmm_api.R`. The file
# `PDfiltering_library.R` contains the lower-level implementation used by
# that interface.

source(here::here("R", "pdhmm_setup.R"))
source(here::here("R", "pdhmm_api.R"))
