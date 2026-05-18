# Reproducibility materials for *Dalla Pria, Ruggiero, Spanò, “Exact inference via quasi-conjugacy in two-parameter Poisson–Dirichlet hidden Markov models”, JASA (forthcoming).*

[Preprint (arXiv)](https://arxiv.org/abs/2512.22098) · [JASA DOI](https://doi.org/10.1080/01621459.2026.2676717)

This directory contains the code, cached computational outputs, and supporting data needed to reproduce the numerical results in the manuscript and supplementary material.

The materials are organised into:

- reusable software for Poisson--Dirichlet hidden Markov model inference;
- a documented public API for the reusable software layer;
- manuscript-specific analysis scripts;
- cached outputs used to regenerate figures quickly;
- narrative notebooks corresponding to the main empirical, simulation, and benchmarking sections.

## Directory structure

- `R/`: reusable helper code and project setup utilities.
- `R/pdhmm_api.R`: documented public API for reusable filtering, smoothing, and simulation methods.
- `analysis/`: manuscript-specific scripts, split by task.
- `notebooks/`: shorter R Markdown notebooks corresponding to different sections of the paper.
- `data/raw/`: raw external data bundled for reproduction of the Infectious analysis.
- `results/cached/`: cached intermediate outputs used to regenerate expensive figures quickly.
- `results/figures/`: figures written by the driver scripts.
- `results/tables/`: tables written by the driver scripts.
- `results/logs/`: session and package version logs.
- `PDfiltering_library.R`: low-level implementation of the methods used by the public API.

## Software requirements

The reference run for this archive used:

- R `4.5.2`
- macOS Sonoma `14.8.1`
- platform `aarch64-apple-darwin20`

The main package versions used in the reference environment were:

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

The complete software-environment records are listed in:

- [REQUIREMENTS.md](REQUIREMENTS.md)
- [results/logs/package_versions.csv](results/logs/package_versions.csv)
- [results/logs/sessionInfo.txt](results/logs/sessionInfo.txt)

The scripts use the `here` package for relative paths. No `setwd()` call is required.

Users who want to reuse the software outside the manuscript workflows should start from `R/pdhmm_api.R`.

## Quick reproduction

For a quick run that regenerates the manuscript figures from cached intermediate outputs whenever available, start an R session in the repository root directory and run:

```r
source("run_quick.R")
```

This driver script:

- checks that the required packages are installed;
- creates output directories if needed;
- sources the split analysis scripts;
- prefers cached `.csv` files in `results/cached/` for the most expensive benchmark outputs;
- saves session information to `results/logs/`.

This is the recommended entry point for an external user who wants to verify the computational materials with minimal effort.

**Note (large cached BPF output).** The cached CSV output for the bootstrap particle filter (BPF) in the Infectious analysis is not hosted in this GitHub repository due to file-size constraints. To reproduce the BPF panel/outputs, please re-run the corresponding chunk/script locally (expected runtime ≈ 1 hour), as documented in the notebooks/scripts.

## Full recomputation from raw data

For a fuller run that renders the section notebooks as standalone documents, start an R session in the `Code` directory and run:

```r
source("run_all.R")
```

This renders the notebooks in `notebooks/`. Since notebook rendering uses `rmarkdown`, a working `pandoc` installation is required; in practice, using RStudio is the simplest option. The most computationally demanding steps remain stochastic and may exhibit small run-to-run variation unless the same random seeds and Monte Carlo budgets are used.

## How an external user should use the archive

An external user should proceed as follows.

1. Unzip the archive and open an R session with working directory equal to `Code`.
2. Install the packages listed in the software-requirements section if they are not already available.
3. Run `source("run_quick.R")` to reproduce the main outputs using the bundled cached files.
4. Inspect the notebooks in `notebooks/` for section-by-section documentation.
5. Optionally run `source("run_all.R")` to render the notebooks as standalone HTML documents.

Users who prefer to run only one part of the workflow can also:

- source one of the scripts in `analysis/`;
- knit one notebook from `notebooks/`.

## Runtime expectations

The timings below are approximate and hardware-dependent.

- `analysis/01_infectious_analysis.R`: about 10--20 minutes if the dual filter and smoother are recomputed; about 1 hour if the bootstrap particle filter is also recomputed rather than read from cache.
- `analysis/02_simulation_study.R`: roughly 20--50 minutes depending on whether the marginal-likelihood grid search is rerun.
- `analysis/03_benchmarks_and_supplement.R`: a few minutes when reading cached supplementary outputs; substantially longer if all benchmark code is rerun from scratch.
- `run_quick.R`: intended for a referee-style run using cached intermediates wherever possible.
- `run_all.R`: intended to render the section notebooks as standalone documents.

## Data provenance and licences

The empirical analysis uses the Infectious SocioPatterns data:

- dataset page: <https://sociopatterns.org/datasets/infectious-exhibition-dynamic-contact-networks/>
- bundled file used in the manuscript: `data/raw/listcontacts_2009_06_28.txt`

**Data license.** The SocioPatterns data are distributed under **CC BY-NC-SA 3.0**. Please acknowledge the SocioPatterns collaboration and cite:

Isella, L., Stehlé, J., Barrat, A., Cattuto, C., Pinton, J.-F., and Van den Broeck, W. (2011). *What’s in a crowd? Analysis of face-to-face behavioral networks.* **Journal of Theoretical Biology**, 271(1), 166–180.

**Code license.** The code in this repository is released under the **MIT License** (see `LICENSE`).

## Figure and output map

The high-level mapping from scripts and cached files to manuscript outputs is documented in:

- [FIGURE_MAP.md](FIGURE_MAP.md)
- [MANIFEST.md](MANIFEST.md)

## Notes on stochastic variation

Some reported quantities are obtained by Monte Carlo procedures. In particular:

- grid-search likelihoods based on the dual particle filter;
- filtering and smoothing summaries based on Monte Carlo draws;
- bootstrap particle-filter outputs.

Accordingly, full recomputation can produce small numerical differences relative to the cached outputs bundled in this archive.

## Main entry points

The recommended entry points for reproduction are:

- `run_quick.R`
- `run_all.R`
- the scripts in `analysis/`
- the shorter notebooks in `notebooks/`
