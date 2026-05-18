# Figure and output map

This file gives a high-level map from manuscript outputs to the scripts in this directory.

## Empirical Infectious analysis

- Contact-network reconstruction and partition extraction:
  `analysis/01_infectious_analysis.R`
- Maximum-likelihood estimation for the empirical dataset:
  `analysis/01_infectious_analysis.R`
- Dual filtering and smoothing plots for heterozygosity:
  `analysis/01_infectious_analysis.R`
- Bootstrap particle-filter comparison:
  `analysis/01_infectious_analysis.R`
  with cached inputs
  `results/cached/my_data_filter_bootstrap_graphs.csv`
  and
  `results/cached/my_quantiles_filter_bootstrap_graphs.csv`

## Simulation study

- Simulated signal and partition observations:
  `analysis/02_simulation_study.R`
- Marginal-likelihood grid search:
  `analysis/02_simulation_study.R`
- Dual filter and smoother on synthetic data:
  `analysis/02_simulation_study.R`
- Synthetic bootstrap particle-filter comparison:
  `analysis/02_simulation_study.R`
- Independent-prior baseline on synthetic data:
  `analysis/02_simulation_study.R`

## Supplementary coordinate plots and benchmark outputs

- First three coordinates, filter and smoother:
  `analysis/03_benchmarks_and_supplement.R`
  using cached files
  `results/cached/my_data_filter_coord*.csv`,
  `results/cached/my_data_smoothing_coord*.csv`,
  `results/cached/my_quantiles_filter_coord*.csv`,
  `results/cached/my_quantiles_smoothing_coord*.csv`,
  `results/cached/my_true_top_coord*.csv`,
  and
  `results/cached/my_true_top_and_Y_coord*.csv`
- Timing, scores, and pruning-related summaries:
  `analysis/03_benchmarks_and_supplement.R`
  using the cached benchmark `.csv` files listed in `MANIFEST.md`

## Notebooks

The notebooks in `notebooks/` mirror the same division:

- `infectious_analysis.Rmd`
- `simulation_study.Rmd`
- `benchmarks.Rmd`

These notebooks are intended for readable section-level inspection, while the scripts in `analysis/` are intended as the main execution units.
