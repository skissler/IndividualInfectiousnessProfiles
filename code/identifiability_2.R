# ==============================================================================
# Identifiability of the punctuation parameter psi
# ==============================================================================
#
# All information about psi comes from within-cluster correlations between
# sibling serial intervals. This script:
#   1. Implements convolution densities for Gamma+Gamma and Gamma-Gamma
#   2. Computes the cluster likelihood via 1D quadrature over shared component
#   3. Runs a simulation study to characterize posterior precision vs data size
#   4. Produces diagnostic and summary figures
#
# Depends on: parslist (from parameters.R), save_fig (from utils.R)
# ==============================================================================

#' Simulate transmission clusters from the full observation model
#'
#' @param K Number of index cases (clusters)
#' @param psi_true True punctuation parameter
#' @param R0 Basic reproduction number
#' @param alpha,beta Generation interval Gamma parameters
#' @param a_obs,b_obs Detection delay Gamma parameters
#' @param p_asc Ascertainment probability
#' @return List of numeric vectors, one per informative cluster (m >= 2),
#'   containing observed serial intervals
simulate_clusters <- function(K, psi_true, R0, alpha, beta, a_obs, b_obs, p_asc) {

	stopifnot(K > 0)
	stopifnot(psi_true >= 0, psi_true <= 1)
	stopifnot(R0 > 0)
	stopifnot(alpha >= 0)
	stopifnot(beta >= 0)
	stopifnot(a_obs >= 0)
	stopifnot(b_obs >= 0)
	stopifnot(p_asc >= 0, p_asc <= 1)

	clusters <- list()
	idx <- 0

	for (i in seq_len(K)) {
		# Number of offspring
		N <- rpois(1, R0)
		if (N == 0) next

		# Ascertainment
		observed <- rbinom(1, N, p_asc)
		if (observed < 2) next  # need >= 2 for information

		# Shared latent component
		l <- rgamma(1, shape = (1 - psi_true) * alpha, rate = beta)

		# Individual jitters
		eps <- rgamma(observed, shape = psi_true * alpha, rate = beta)

		# Generation intervals
		tau <- l + eps

		# Detection delays
		d0 <- rgamma(1, shape = a_obs, rate = b_obs)  # index case delay
		dj <- rgamma(observed, shape = a_obs, rate = b_obs)  # offspring delays

		# Serial intervals
		s <- tau + dj - d0

		idx <- idx + 1
		clusters[[idx]] <- s
	}

	clusters
}

# ==============================================================================
# Stan-based inference for psi
# ==============================================================================

#' Prepare cluster data for Stan
#'
#' Converts the list-of-vectors format from simulate_clusters() into the named
#' list that the Stan model expects.
#'
#' @param clusters List of numeric vectors (each cluster's serial intervals)
#' @param alpha Generation interval Gamma shape
#' @param beta Generation interval Gamma rate
#' @param a_obs Detection delay Gamma shape
#' @param b_obs Detection delay Gamma rate
#' @return Named list suitable for passing to Stan's sampling()
prepare_stan_data <- function(clusters, alpha, beta, a_obs, b_obs) {
	K <- length(clusters)
	m <- vapply(clusters, length, integer(1))
	N <- sum(m)
	s <- unlist(clusters)

	# Build cluster membership index
	cluster_idx <- rep(seq_len(K), times = m)

	list(
		K        = K,
		m        = m,
		N        = N,
		s        = s,
		cluster  = cluster_idx,
		alpha_gi = alpha,
		beta_gi  = beta,
		a_obs    = a_obs,
		b_obs    = b_obs
	)
}

#' Fit the psi inference Stan model
#'
#' Compiles the Stan model (cached), prepares data, runs MCMC sampling,
#' and returns the stanfit object.
#'
#' @param clusters List of numeric vectors (serial intervals per cluster)
#' @param alpha Generation interval Gamma shape
#' @param beta Generation interval Gamma rate
#' @param a_obs Detection delay Gamma shape
#' @param b_obs Detection delay Gamma rate
#' @param ... Additional arguments passed to rstan::sampling() (e.g. chains, iter, warmup)
#' @return A stanfit object
fit_psi <- function(clusters, alpha, beta, a_obs, b_obs, ...) {
	if (!requireNamespace("rstan", quietly = TRUE))
		stop("Package 'rstan' is required for fit_psi()")

	stan_file <- file.path("code", "psiinference.stan")
	if (!file.exists(stan_file))
		stop("Stan model not found at: ", stan_file)

	stan_data <- prepare_stan_data(clusters, alpha, beta, a_obs, b_obs)

	# Generate initial values guaranteed to satisfy d[n] > 0 for all n.
	# Since d[n] = s[n] - l[cluster[n]] - eps[n] + d0[cluster[n]],
	# we set l and eps small, and d0 per-cluster large enough that
	# d[n] >= 1.0 for the worst-case observation in each cluster.
	init_fn <- function() {
		K <- stan_data$K
		N <- stan_data$N
		s_vec <- stan_data$s
		cl <- stan_data$cluster

		l_init   <- rep(0.1, K)
		eps_init <- rep(0.1, N)
		d0_init  <- numeric(K)

		for (k in seq_len(K)) {
			idx <- which(cl == k)
			min_s_k <- min(s_vec[idx])
			# Need d0[k] > l[k] + eps_max - min(s) = 0.1 + 0.1 - min_s_k
			d0_init[k] <- max(a_obs / b_obs, 0.2 - min_s_k + 1.0)
		}

		list(
			psi = runif(1, 0.3, 0.7),
			l   = l_init,
			d0  = d0_init,
			eps = eps_init
		)
	}

	model <- rstan::stan_model(stan_file)
	rstan::sampling(model, data = stan_data, init = init_fn, ...)
}

# ==============================================================================
# Verification: simulate at known psi, check posterior recovery
# ==============================================================================

#' Quick verification test for the Stan psi inference model
#'
#' Simulates clusters at a known psi_true, fits the model, and checks that the
#' posterior concentrates around the true value.
#'
#' @param psi_true True punctuation parameter to simulate from
#' @param K Number of index cases to simulate
#' @param R0 Basic reproduction number
#' @param alpha Generation interval Gamma shape
#' @param beta Generation interval Gamma rate
#' @param a_obs Detection delay Gamma shape
#' @param b_obs Detection delay Gamma rate
#' @param p_asc Ascertainment probability
#' @param seed Random seed for reproducibility
#' @param ... Additional arguments passed to fit_psi() / rstan::sampling()
#' @return Invisible list with fit object and summary
verify_psi_inference <- function(psi_true = 0.5, K = 500, R0 = 2,
                                  alpha = 2.32, beta = 0.725,
                                  a_obs = 4, b_obs = 1, p_asc = 1.0,
                                  seed = 42, ...) {
	cat(sprintf("=== Verification: psi_true = %.2f, K = %d ===\n", psi_true, K))

	set.seed(seed)
	clusters <- simulate_clusters(K, psi_true, R0, alpha, beta, a_obs, b_obs, p_asc)
	cat(sprintf("  Simulated %d informative clusters (of %d index cases)\n",
	    length(clusters), K))

	if (length(clusters) < 5) {
		cat("  Too few clusters for meaningful test. Increase K.\n")
		return(invisible(NULL))
	}

	fit <- fit_psi(clusters, alpha, beta, a_obs, b_obs,
	               chains = 4, iter = 2000, warmup = 1000, seed = seed, ...)

	psi_draws <- rstan::extract(fit, "psi")$psi
	post_mean   <- mean(psi_draws)
	post_median <- median(psi_draws)
	ci <- quantile(psi_draws, c(0.025, 0.975))
	covers <- (psi_true >= ci[1]) && (psi_true <= ci[2])

	cat(sprintf("  Posterior mean:   %.3f\n", post_mean))
	cat(sprintf("  Posterior median: %.3f\n", post_median))
	cat(sprintf("  95%% CI:          (%.3f, %.3f)\n", ci[1], ci[2]))
	cat(sprintf("  CI width:         %.3f\n", ci[2] - ci[1]))
	cat(sprintf("  Covers true psi:  %s\n", ifelse(covers, "YES", "NO")))

	invisible(list(fit = fit, psi_true = psi_true,
	               mean = post_mean, median = post_median,
	               ci_lo = ci[1], ci_hi = ci[2], covers = covers))
}