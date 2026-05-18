# Simulation-study workflow for the manuscript.
#
# Purpose:
# - simulate a PD diffusion trajectory and the associated partition-valued data;
# - estimate (alpha, theta) by marginal likelihood on a grid or reuse the
#   manuscript values;
# - reproduce the synthetic filtering, smoothing, bootstrap, and
#   independent-prior comparisons.
#
# This script intentionally retains the same mathematical workflow as the
# original `PDfiltering.Rmd`, but isolates it from the empirical analysis.

source(here::here("analysis", "00_setup.R"))

set.seed(223)

my_alpha <- 0.1
my_theta <- 1.5
my_steps <- 19
my_T <- my_steps + 1
my_t <- 0.025
my_epsilon <- 0.005
obs <- 50
show_top <- 3
my_q <- 2

signal <- Petrov(my_alpha, my_theta, my_steps, my_t, my_epsilon, N = obs, show_top = Inf)
true_entropy <- sapply(signal, function(el) qentropy(el$x / sum(el$x), q = my_q))
true_top <- sapply(signal, function(el) el$x[1:show_top])
Y <- lapply(signal, function(x) x$gamma[[1]])

run_simulation_mle <- FALSE
if (run_simulation_mle) {
  min_alpha <- 0
  max_alpha <- 0.2
  min_theta <- 0.75
  max_theta <- 2.25
  grid_alpha <- seq(min_alpha, max_alpha, by = 0.05)
  grid_theta <- seq(min_theta, max_theta, by = 0.25)
  my_grid <- expand.grid(alpha = grid_alpha, theta = grid_theta)
  my_grid <- my_grid[-1, ]
  temp <- my_grid
  temp$psf <- pbapply(
    my_grid, 1,
    function(el) {
      sum(log(Particle_filter_likelihood(
        Y,
        alpha = el[1],
        theta = el[2],
        t = my_t,
        pruning_update = 0,
        pruning_propagation = 10,
        N = 0,
        nParticles = 10^4,
        store_mixture = FALSE
      )))
    }
  )
  my_MLEs <- temp[which.max(temp$psf), ]
} else {
  my_MLEs <- c(0, 1.5)
}

MCrep <- 10^4
MCparticles <- 10^4
my_pruning_propagation <- 10

Y_ts <- sapply(Y, function(el) qentropy(el / sum(el), q = 2))
Y_ts <- data.frame(data = Y_ts, type = as.factor("Observed entropy"), time = as.factor(seq(from = 0, length.out = my_steps + 1, by = my_t)))
true_entropy_df <- data.frame(data = as.numeric(true_entropy), type = as.factor("Signal entropy"), time = as.factor(seq(from = 0, length.out = my_steps + 1, by = my_t)))
true_entropy_overlay <- rbind(true_entropy_df, Y_ts)

my_filter <- Particle_filter_qentropy(
  Y,
  alpha = as.numeric(my_MLEs[1]),
  theta = as.numeric(my_MLEs[2]),
  t = my_t,
  epsilon = my_epsilon,
  pruning_update = 0,
  pruning_propagation = my_pruning_propagation,
  N = MCrep,
  nParticles = MCparticles,
  show_top = Inf,
  verbose = FALSE,
  q = 2
)

quantiles_filter <- sapply(seq_len(my_T), function(t) quantile(my_filter[[t]], probs = c(0.025, 0.975), na.rm = TRUE))
my_data_filter <- lapply(seq_len(my_T), function(t) data.frame(data = my_filter[[t]], type = rep("filter", MCrep), time = (t - 1) * my_t))
my_data_filter <- do.call("rbind", my_data_filter)
my_data_filter$time <- as.factor(my_data_filter$time)
my_quantiles_filter <- data.frame(lower = quantiles_filter[1, ], upper = quantiles_filter[2, ], time = as.factor(seq(from = 0, length.out = my_steps + 1, by = my_t)))

my_smooth <- Smoothing_qentropy(
  Y,
  alpha = as.numeric(my_MLEs[1]),
  theta = as.numeric(my_MLEs[2]),
  t = my_t,
  epsilon = my_epsilon,
  pruning_update = 10,
  pruning_propagation = my_pruning_propagation,
  N = MCrep,
  nParticles = MCparticles,
  show_top = Inf,
  verbose = FALSE,
  q = 2
)[[2]]

quantiles_smooth <- sapply(seq_len(my_T), function(t) quantile(my_smooth[[t]], probs = c(0.025, 0.975), na.rm = TRUE))
my_data_smooth <- lapply(seq_len(my_T), function(t) data.frame(data = my_smooth[[t]], type = rep("smooth", MCrep), time = (t - 1) * my_t))
my_data_smooth <- do.call("rbind", my_data_smooth)
my_data_smooth$time <- as.factor(my_data_smooth$time)
my_quantiles_smooth <- data.frame(lower = quantiles_smooth[1, ], upper = quantiles_smooth[2, ], time = as.factor(seq(from = 0, length.out = my_steps + 1, by = my_t)))

plot_filter <- ggplot() +
  geom_violin(data = my_data_filter, aes(y = data, x = time), bw = 0.01, fill = "orchid1", colour = "orchid1", alpha = 0.66, scale = "count") +
  geom_line(data = my_quantiles_filter, aes(y = upper, x = time, group = 1), linetype = "dashed", colour = "blue", linewidth = 1) +
  geom_line(data = my_quantiles_filter, aes(y = lower, x = time, group = 1), linetype = "dashed", colour = "blue", linewidth = 1) +
  geom_line(data = true_entropy_df, aes(y = data, x = time, group = 1), linewidth = 1) +
  geom_point(data = true_entropy_overlay, aes(y = data, x = time, shape = type, fill = type), size = 3, color = "black", show.legend = FALSE) +
  scale_shape_manual(values = c(21, 23)) +
  scale_fill_manual(values = c("black", "white")) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("") +
  ylab("")

plot_smooth <- ggplot() +
  geom_violin(data = my_data_smooth, aes(y = data, x = time), bw = 0.01, fill = "orange", colour = "orange", alpha = 0.66, scale = "count") +
  geom_line(data = my_quantiles_smooth, aes(y = upper, x = time, group = 1), linetype = "dashed", colour = "red", linewidth = 1) +
  geom_line(data = my_quantiles_smooth, aes(y = lower, x = time, group = 1), linetype = "dashed", colour = "red", linewidth = 1) +
  geom_line(data = true_entropy_df, aes(y = data, x = time, group = 1), linewidth = 1) +
  geom_point(data = true_entropy_overlay, aes(y = data, x = time, shape = type, fill = type), size = 3, color = "black", show.legend = FALSE) +
  scale_shape_manual(values = c(21, 23)) +
  scale_fill_manual(values = c("black", "white")) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("") +
  ylab("")

ggsave(figure_path("simulation_filtering.png"), plot_filter, width = 11, height = 5, dpi = 300, bg = 'white')
ggsave(figure_path("simulation_smoothing.png"), plot_smooth, width = 11, height = 5, dpi = 300, bg = 'white')

safe_notify("Simulation study finished")
