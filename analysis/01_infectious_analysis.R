# Empirical analysis of the Infectious SocioPatterns data.
#
# Purpose:
# - reconstruct the partition-valued observations from the bundled raw data;
# - estimate the global parameters (alpha, theta) or reuse the fixed values
#   employed in the manuscript;
# - reproduce the main filtering, smoothing, and benchmark plots.
#
# Inputs:
# - data/raw/listcontacts_2009_06_28.txt
# - results/cached/my_data_filter_bootstrap_graphs.csv
# - results/cached/my_quantiles_filter_bootstrap_graphs.csv
#
# Outputs:
# - figure files in results/figures/
#
# Expected runtime:
# - about 10--20 minutes for the dual-filter and smoother computations;
# - about 1 hour if the bootstrap particle filter is recomputed from scratch.

if (file.exists(here::here("analysis", "00_setup.R"))) {
  source(here::here("analysis", "00_setup.R"))   
} else {
  source("00_setup.R")                           
}

day <- 28
contact_file <- raw_data_path(sprintf("listcontacts_2009_06_%02d.txt", day))
listcontacts <- read.delim(contact_file, header = FALSE)

# Sort edge endpoints so that each undirected contact has a canonical order.
listcontacts_sorted <- lapply(
  seq_len(nrow(listcontacts)),
  function(row) sort(as.numeric(listcontacts[row, c(2, 3)]), decreasing = FALSE)
)
listcontacts <- cbind(listcontacts[, 1], do.call(rbind, listcontacts_sorted))
colnames(listcontacts) <- c("time", "from", "to")
listcontacts <- as.data.frame(listcontacts)
mode(listcontacts$from) <- "character"
mode(listcontacts$to) <- "character"

listcontacts$hour <- vapply(listcontacts$time, function(el) hour(as_datetime(el)), numeric(1))
listcontacts$minute <- vapply(listcontacts$time, function(el) minute(as_datetime(el)), numeric(1))
listcontacts <- listcontacts[, -1]

by_hour <- lapply(unique(listcontacts$hour), function(t) subset(listcontacts, hour == t))
by_window <- lapply(by_hour, function(block) {
  if (block[1, 4] < 30 && block[nrow(block), 4] >= 30) {
    list(subset(block, minute < 30), subset(block, minute >= 30))
  } else {
    list(block)
  }
})
by_window <- unlist(by_window, recursive = FALSE)
labels <- vapply(
  by_window,
  function(block) paste0(block[1, 3], "h ", if (block[1, 4] < 30) "\n<30m" else "\n>=30m"),
  character(1)
)
edge_lists <- lapply(by_window, function(block) unique(block[, c(1, 2)]))
names(edge_lists) <- labels

for (i in seq_along(edge_lists)) {
  edge_lists[[i]]$from <- as.character(edge_lists[[i]]$from)
  edge_lists[[i]]$to <- as.character(edge_lists[[i]]$to)
}

graphs <- lapply(edge_lists, function(el) {
  if (nrow(el) > 1) {
    graph_from_edgelist(as.matrix(el[, c(1, 2)]), directed = FALSE)
  } else {
    NULL
  }
})
graphs <- graphs[!vapply(graphs, is.null, logical(1))]
Y_graphs <- lapply(graphs, function(g) sort(components(g)$csize, decreasing = TRUE))
graph_labels <- names(Y_graphs)

my_steps <- length(Y_graphs) - 1
my_T <- my_steps + 1
my_t <- 0.1
my_epsilon <- 0.005

# The manuscript uses fixed values unless the expensive grid-search MLE is rerun.
run_empirical_mle <- FALSE
if (run_empirical_mle) {
  min_alpha <- 0
  max_alpha <- 0.05
  min_theta <- 0.7
  max_theta <- 0.9
  grid_alpha <- seq(min_alpha, max_alpha, by = 0.01)
  grid_theta <- seq(min_theta, max_theta, by = 0.05)
  my_grid <- expand.grid(alpha = grid_alpha, theta = grid_theta)
  my_grid <- my_grid[-1, ]
  temp_graphs <- my_grid
  temp_graphs$psf <- pbapply(
    my_grid, 1,
    function(el) {
      sum(log(Particle_filter_likelihood(
        Y_graphs,
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
  my_MLEs_graphs <- temp_graphs[which.max(temp_graphs$psf), ]
} else {
  my_MLEs_graphs <- c(0, 0.75)
}

MCrep <- 10^4
MCparticles <- 10^5
my_q <- 2

entropy_graphs <- sapply(Y_graphs, function(el) qentropy(el / sum(el), q = my_q))
entropy_graphs <- data.frame(
  data = entropy_graphs,
  type = as.factor("Observed entropy"),
  time = as.factor(seq(from = 0, length.out = my_steps + 1, by = my_t))
)

my_filter_graphs <- Particle_filter_qentropy(
  Y_graphs,
  alpha = as.numeric(my_MLEs_graphs[1]),
  theta = as.numeric(my_MLEs_graphs[2]),
  t = my_t,
  epsilon = my_epsilon,
  pruning_update = 0,
  pruning_propagation = 10,
  N = MCrep,
  nParticles = MCparticles,
  show_top = Inf,
  verbose = FALSE,
  q = my_q
)

quantiles_filter <- sapply(
  seq_len(my_T),
  function(time) quantile(my_filter_graphs[[time]], probs = c(0.025, 0.975), na.rm = TRUE)
)

my_data_filter_graphs <- lapply(
  seq_len(my_T),
  function(time) {
    data.frame(data = my_filter_graphs[[time]], type = rep("filter", MCrep), time = (time - 1) * my_t)
  }
)
my_data_filter_graphs <- do.call("rbind", my_data_filter_graphs)
my_data_filter_graphs$time <- as.factor(my_data_filter_graphs$time)

my_quantiles_filter_graphs <- data.frame(
  lower = quantiles_filter[1, ],
  upper = quantiles_filter[2, ],
  time = as.factor(seq(from = 0, length.out = my_steps + 1, by = my_t))
)

my_smooth_graphs <- Smoothing_qentropy(
  Y_graphs,
  alpha = as.numeric(my_MLEs_graphs[1]),
  theta = as.numeric(my_MLEs_graphs[2]),
  t = my_t,
  epsilon = my_epsilon,
  pruning_update = 10,
  pruning_propagation = 10,
  N = MCrep,
  nParticles = MCparticles,
  show_top = Inf,
  verbose = FALSE,
  q = my_q
)[[2]]

quantiles_smooth <- sapply(
  seq_len(my_T),
  function(time) quantile(my_smooth_graphs[[time]], probs = c(0.025, 0.975), na.rm = TRUE)
)

my_data_smooth_graphs <- lapply(
  seq_len(my_T),
  function(time) {
    data.frame(data = my_smooth_graphs[[time]], type = rep("smooth", MCrep), time = (time - 1) * my_t)
  }
)
my_data_smooth_graphs <- do.call("rbind", my_data_smooth_graphs)
my_data_smooth_graphs$time <- as.factor(my_data_smooth_graphs$time)

my_quantiles_smooth_graphs <- data.frame(
  lower = quantiles_smooth[1, ],
  upper = quantiles_smooth[2, ],
  time = as.factor(seq(from = 0, length.out = my_steps + 1, by = my_t))
)

# The bootstrap benchmark is expensive, so this script reads the cached
# output by default and only recomputes it if explicitly requested.
run_empirical_bpf <- FALSE
if (run_empirical_bpf) {
  my_filter_bootstrap_graphs <- Bootstrap_filter_King_qentropy(
    Y_graphs,
    alpha = as.numeric(my_MLEs_graphs[1]),
    theta = as.numeric(my_MLEs_graphs[2]),
    delta = my_t,
    epsilon = my_epsilon,
    nParticles = 10^4,
    q = my_q
  )
  my_data_filter_bootstrap_graphs <- lapply(
    seq_len(my_T),
    function(time) {
      data.frame(
        data = my_filter_bootstrap_graphs$Filter[[time]],
        importance = my_filter_bootstrap_graphs$Importance_weights[[time]],
        type = rep("bootstrap", length(my_filter_bootstrap_graphs$Filter[[time]])),
        time = (time - 1) * my_t
      )
    }
  )
  my_data_filter_bootstrap_graphs <- do.call("rbind", my_data_filter_bootstrap_graphs)
  my_data_filter_bootstrap_graphs$time <- as.factor(my_data_filter_bootstrap_graphs$time)
} else {
  my_data_filter_bootstrap_graphs <- read.csv(cached_results_path("my_data_filter_bootstrap_graphs.csv"))
  my_data_filter_bootstrap_graphs$time <- as.factor(my_data_filter_bootstrap_graphs$time)
  my_quantiles_filter_bootstrap_graphs <- read.csv(cached_results_path("my_quantiles_filter_bootstrap_graphs.csv"))
  my_quantiles_filter_bootstrap_graphs$time <- as.factor(my_quantiles_filter_bootstrap_graphs$time)
}

plot_filter <- ggplot() +
  geom_violin(data = my_data_filter_graphs, aes(y = data, x = time), bw = 0.01, fill = "orchid1", colour = "orchid1", alpha = 0.66, scale = "width") +
  geom_line(data = entropy_graphs, aes(y = data, x = time, group = 1), linewidth = 1, color = "blue") +
  geom_point(data = entropy_graphs, aes(y = data, x = time), size = 3, fill = "white", color = "black", shape = 23) +
  geom_point(data = my_quantiles_filter_graphs, aes(y = upper, x = time), shape = 25, fill = "blue", color = "white", size = 3) +
  geom_point(data = my_quantiles_filter_graphs, aes(y = lower, x = time), shape = 24, fill = "blue", color = "white", size = 3) +
  ylim(c(0, 1)) +
  scale_x_discrete(labels = graph_labels) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("") +
  ylab("")

plot_smooth <- ggplot() +
  geom_violin(data = my_data_smooth_graphs, aes(y = data, x = time), bw = 0.01, fill = "orange", colour = "orange", alpha = 0.66, scale = "width") +
  geom_line(data = entropy_graphs, aes(y = data, x = time, group = 1), linewidth = 1, color = "red") +
  geom_point(data = entropy_graphs, aes(y = data, x = time), size = 3, fill = "white", color = "black", shape = 23) +
  geom_point(data = my_quantiles_smooth_graphs, aes(y = upper, x = time), shape = 25, fill = "red", color = "white", size = 3) +
  geom_point(data = my_quantiles_smooth_graphs, aes(y = lower, x = time), shape = 24, fill = "red", color = "white", size = 3) +
  ylim(c(0, 1)) +
  scale_x_discrete(labels = graph_labels) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("") +
  ylab("")

plot_bootstrap <- ggplot() +
  geom_violin(data = my_data_filter_bootstrap_graphs, aes(y = data, x = time, weight = importance), bw = 0.01, fill = "hotpink", colour = "hotpink", alpha = 0.5, scale = "width") +
  geom_line(data = entropy_graphs, aes(y = data, x = time, group = 1), linewidth = 1, color = "purple") +
  geom_point(data = entropy_graphs, aes(y = data, x = time), size = 3, fill = "white", color = "black", shape = 23) +
  geom_point(data = my_quantiles_filter_bootstrap_graphs, aes(y = upper, x = time), shape = 25, fill = "purple", color = "white", size = 3) +
  geom_point(data = my_quantiles_filter_bootstrap_graphs, aes(y = lower, x = time), shape = 24, fill = "purple", color = "white", size = 3) +
  ylim(c(0, 1)) +
  scale_x_discrete(labels = graph_labels) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  xlab("") +
  ylab("")

ggsave(figure_path("infectious_filtering.png"), plot_filter, width = 11, height = 5, dpi = 300, bg = 'white')
ggsave(figure_path("infectious_smoothing.png"), plot_smooth, width = 11, height = 5, dpi = 300, bg = 'white')
ggsave(figure_path("infectious_bootstrap.png"), plot_bootstrap, width = 11, height = 5, dpi = 300, bg = 'white')

safe_notify("Empirical analysis finished")
