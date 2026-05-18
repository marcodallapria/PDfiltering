# Requirements

This file records the software environment used for the reference run of the reproducibility archive.

## Reference environment

- R version: `4.5.2`
- OS: `macOS Sonoma 14.8.1`
- Platform: `aarch64-apple-darwin20`

## Required R packages

The following package versions were present in the reference environment.

| Package | Version |
| --- | --- |
| here | 1.0.1 |
| DescTools | 0.99.60 |
| partitions | 1.10-9 |
| collections | 0.3.9 |
| LaplacesDemon | 16.1.6 |
| extraDistr | 1.10.0 |
| scales | 1.4.0 |
| vioplot | 0.5.1 |
| pbapply | 1.7-4 |
| ggplot2 | 3.5.2 |
| RColorBrewer | 1.1-3 |
| patchwork | 1.3.2 |
| RcppAlgos | 2.9.3 |
| microbenchmark | 1.5.0 |
| dplyr | 1.1.4 |
| pbmcapply | 1.5.1 |
| R.utils | 2.13.0 |
| geomtextpath | 0.2.0 |
| igraph | 2.1.4 |
| lubridate | 1.9.4 |
| knitr | 1.50 |
| rmarkdown | 2.29 |

The package `renv` was not installed in the local reference environment used for this revision, so no `renv.lock` file is included here. Instead, the exact package versions are recorded in:

- `results/logs/package_versions.csv`
- `results/logs/sessionInfo.txt`

## Installation

If the packages above are not already installed, they can be installed with:

```r
install.packages(c(
  "here", "DescTools", "partitions", "collections", "LaplacesDemon",
  "extraDistr", "scales", "vioplot", "pbapply", "ggplot2",
  "RColorBrewer", "patchwork", "RcppAlgos", "microbenchmark",
  "dplyr", "pbmcapply", "R.utils", "geomtextpath", "igraph",
  "lubridate", "knitr", "rmarkdown"
))
```
