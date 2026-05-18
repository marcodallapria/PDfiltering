# Benchmark and supplementary-output script.
#
# Purpose:
# - read the cached outputs used in the supplement;
# - regenerate the first-three-coordinate plots;
# - centralise the benchmark-oriented cached files in one place.
#
# This script focuses on the supplementary material that relies on cached
# `.csv` files bundled with the submission.

source(here::here("R", "pdhmm_setup.R"))
check_required_packages(c("here", "ggplot2"))
ensure_output_directories()
library(ggplot2)

labs_x_title <- c("st", "nd", "rd", "th", "th")
show_top <- 3

for (j in seq_len(show_top)) {
  my_data_filter_coord <- read.csv(cached_results_path(sprintf("my_data_filter_coord%d.csv", j)))
  my_data_filter_coord$time <- as.factor(my_data_filter_coord$time)

  my_data_smoothing_coord <- read.csv(cached_results_path(sprintf("my_data_smoothing_coord%d.csv", j)))
  my_data_smoothing_coord$time <- as.factor(my_data_smoothing_coord$time)

  my_true_top_coord <- read.csv(cached_results_path(sprintf("my_true_top_coord%d.csv", j)))
  my_true_top_coord$time <- as.factor(my_true_top_coord$time)
  my_true_top_coord$type <- as.factor(my_true_top_coord$type)

  my_true_top_and_Y_coord <- read.csv(cached_results_path(sprintf("my_true_top_and_Y_coord%d.csv", j)))
  my_true_top_and_Y_coord$time <- as.factor(my_true_top_and_Y_coord$time)
  my_true_top_and_Y_coord$type <- as.factor(my_true_top_and_Y_coord$type)

  my_quantiles_filter_coord <- read.csv(cached_results_path(sprintf("my_quantiles_filter_coord%d.csv", j)))
  my_quantiles_filter_coord$time <- as.factor(my_quantiles_filter_coord$time)

  my_quantiles_smoothing_coord <- read.csv(cached_results_path(sprintf("my_quantiles_smoothing_coord%d.csv", j)))
  my_quantiles_smoothing_coord$time <- as.factor(my_quantiles_smoothing_coord$time)

  plot_filter_coord <- ggplot() +
    geom_violin(data = my_data_filter_coord, aes(y = data, x = time), bw = 0.01, fill = "orchid1", colour = "orchid1", alpha = 1, scale = "count") +
    geom_line(data = my_quantiles_filter_coord, aes(y = upper, x = time, group = 1), linetype = "dashed", colour = "blue", linewidth = 1) +
    geom_line(data = my_quantiles_filter_coord, aes(y = lower, x = time, group = 1), linetype = "dashed", color = "blue", linewidth = 1) +
    geom_line(data = my_true_top_coord, aes(y = data, x = time, group = 1), linewidth = 1) +
    geom_point(data = my_true_top_and_Y_coord, aes(y = data, x = time, shape = type, fill = type), size = 2.5, show.legend = FALSE) +
    scale_fill_manual(values = c("white", "black")) +
    scale_shape_manual(values = c(23, 21)) +
    theme_minimal() +
    theme(text = element_text(size = 15)) +
    xlab("") +
    ylab("") +
    ylim(0, 1 / j)

  plot_smoothing_coord <- ggplot() +
    geom_violin(data = my_data_smoothing_coord, aes(y = data, x = time), bw = 0.01, fill = "orange", colour = "orange", alpha = 1, scale = "count") +
    geom_line(data = my_quantiles_smoothing_coord, aes(y = upper, x = time, group = 1), linetype = "dashed", colour = "red", linewidth = 1) +
    geom_line(data = my_quantiles_smoothing_coord, aes(y = lower, x = time, group = 1), linetype = "dashed", colour = "red", linewidth = 1) +
    geom_line(data = my_true_top_coord, aes(y = data, x = time, group = 1), linewidth = 1) +
    geom_point(data = my_true_top_and_Y_coord, aes(y = data, x = time, shape = type, fill = type), size = 2.5, color = "black", show.legend = FALSE) +
    scale_fill_manual(values = c("white", "black")) +
    scale_shape_manual(values = c(23, 21)) +
    theme_minimal() +
    theme(text = element_text(size = 15)) +
    xlab("") +
    ylab("") +
    ylim(0, 1 / j)

  ggsave(figure_path(sprintf("supplement_coord%d_filter.png", j)), plot_filter_coord, width = 8, height = 4.5, dpi = 300, bg = 'white')
  ggsave(figure_path(sprintf("supplement_coord%d_smoothing.png", j)), plot_smoothing_coord, width = 8, height = 4.5, dpi = 300, bg = 'white')
}

# Read benchmark-oriented cached outputs so that the script documents their role
# and fails early if required files are missing.
benchmark_files <- c(
  "my_MLEs_filter3.csv", "my_score_filter3.csv", "my_score_smooth3.csv",
  "my_score_boot3.csv", "my_score_ind.csv", "my_times.csv", "my_scores.csv",
  "size_top_propagation.csv", "size_top_update.csv",
  "top_propagation_relative.csv", "top_update_relative.csv"
)

benchmark_cache <- lapply(benchmark_files, function(file) read.csv(cached_results_path(file)))
names(benchmark_cache) <- benchmark_files

safe_notify("Supplementary benchmark plots finished")
