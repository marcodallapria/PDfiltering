# Public reusable API for the PD-HMM software layer.
#
# This file provides a documented interface for users who want to reuse the
# methods outside the manuscript scripts. The low-level research
# implementation remains in `PDfiltering_library.R`; the functions below give
# stable entry points with explicit names and arguments.

source(here::here("R", "pdhmm_setup.R"))
check_required_packages()
source(here::here("PDfiltering_library.R"))

# Sample an epsilon-truncated two-parameter Poisson--Dirichlet law.
sample_pd <- function(alpha, theta, epsilon = 0.005, show_top = Inf, min_length = 0) {
  pd(alpha = alpha, theta = theta, epsilon = epsilon, show_top = show_top, min_length = min_length)
}

# Sample from the conditional PD posterior given an observed partition.
sample_pd_posterior <- function(alpha, theta, partition, epsilon = 0.005, show_top = Inf, min_length = 0) {
  pd_posterior(
    alpha = alpha,
    theta = theta,
    partition = partition,
    epsilon = epsilon,
    show_top = show_top,
    min_length = min_length
  )
}

# Simulate a discrete-time PD diffusion path together with sampled partitions.
simulate_petrov_path <- function(alpha, theta, steps, delta, epsilon = 0.005, sample_sizes = c(10, 20, 30), show_top = Inf, start = NA, verbose = FALSE) {
  Petrov(
    alpha = alpha,
    theta = theta,
    steps = steps,
    delta = delta,
    epsilon = epsilon,
    N = sample_sizes,
    show_top = show_top,
    start = start,
    verbose = verbose
  )
}

# Run the dual particle filter and return draws of a q-entropy summary.
dual_filter_qentropy <- function(observation_list, alpha, theta, t, epsilon = 0.005, pruning_update = 0, pruning_propagation = 0, n_draws = 10^4, n_particles = 10^4, q = 2, verbose = FALSE) {
  Particle_filter_qentropy(
    observation_list = observation_list,
    alpha = alpha,
    theta = theta,
    t = t,
    epsilon = epsilon,
    pruning_update = pruning_update,
    pruning_propagation = pruning_propagation,
    N = n_draws,
    nParticles = n_particles,
    show_top = Inf,
    verbose = verbose,
    q = q
  )
}

# Run the dual smoother and return draws of a q-entropy summary.
dual_smoother_qentropy <- function(observation_list, alpha, theta, t, epsilon = 0.005, pruning_update = 0, pruning_propagation = 0, n_draws = 10^4, n_particles = 10^4, q = 2, verbose = FALSE) {
  Smoothing_qentropy(
    observation_list = observation_list,
    alpha = alpha,
    theta = theta,
    t = t,
    epsilon = epsilon,
    pruning_update = pruning_update,
    pruning_propagation = pruning_propagation,
    N = n_draws,
    nParticles = n_particles,
    show_top = Inf,
    verbose = verbose,
    q = q
  )[[2]]
}

# Run the bootstrap particle filter on the latent state and return q-entropy summaries.
bootstrap_filter_qentropy <- function(observation_list, alpha, theta, delta, epsilon = 0.005, n_particles = 10^4, q = 2) {
  Bootstrap_filter_King_qentropy(
    observation_list = observation_list,
    alpha = alpha,
    theta = theta,
    delta = delta,
    epsilon = epsilon,
    nParticles = n_particles,
    show_top = Inf,
    q = q
  )
}
