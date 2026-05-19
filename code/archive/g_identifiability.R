# ==============================================================================
# GI parameter identifiability under punctuated infectiousness
# ==============================================================================
#
# Demonstrates that the standard iid serial-interval likelihood produces
# overconfident posteriors for GI shape (alpha) and rate (beta) when
# infectiousness is punctuated (psi < 1), because it ignores within-cluster
# correlation among sibling serial intervals.
#
# Approach:
#   1. Compute marginal serial interval density f_S(s; alpha, beta, a_obs, b_obs)
#      via numerical convolution (same for all psi)
#   2. Evaluate the iid posterior on a 2D grid over (alpha, beta)
#   3. Repeat across simulated datasets at different psi values
#   4. Show coverage of 95% CIs drops as psi -> 0
#
# Depends on: parslist (from parameters.R), save_fig, simulate_clusters (from utils.R)
# ==============================================================================

cat("=== GI parameter identifiability ===\n")

# ==============================================================================
# 1. Convolution density helpers
# ==============================================================================

#' Density of X - Y where X ~ Gamma(shape1, rate1), Y ~ Gamma(shape2, rate2)
dgamma_diff_gi <- function(d, shape1, rate1, shape2, rate2) {
	sapply(d, function(di) {
		integrand <- function(y) {
			dgamma(di + y, shape1, rate1) * dgamma(y, shape2, rate2)
		}
		lower <- max(0, -di)
		tryCatch(
			integrate(integrand, lower = lower, upper = Inf,
			          rel.tol = 1e-8, abs.tol = 1e-12)$value,
			error = function(e) 0
		)
	})
}

# ==============================================================================
# 2. Optimized serial interval density via matrix convolution
# ==============================================================================
#
# The serial interval is S = tau + (d_j - d_0), where:
#   tau ~ Gamma(alpha, beta)     (generation interval)
#   d_j, d_0 ~ Gamma(a_obs, b_obs)  (detection delays)
#
# The density is: f_S(s) = integral f_tau(tau) * f_{d-d0}(s - tau) dtau
#
# Strategy: pre-tabulate f_{d-d0} on a fine grid (fixed, known params),
# then for each (alpha, beta), form g(tau_k) = dgamma(tau_k; alpha, beta)
# and compute f_S(s_j) = dtau * sum_k D[j,k] * g[k], where
# D[j,k] = f_{d-d0}(s_j - tau_k) is pre-computed once.
# ==============================================================================

#' Build the convolution machinery for fast serial interval density evaluation
#'
#' @param a_obs Detection delay Gamma shape
#' @param b_obs Detection delay Gamma rate
#' @param alpha_range Range of alpha values to support
#' @param beta_range Range of beta values to support
#' @return List with D_matrix, tau_grid, s_grid, dtau, and an evaluation function
build_serial_density_engine <- function(a_obs, b_obs, alpha_range, beta_range) {
	sd_d <- sqrt(a_obs) / b_obs

	# Determine grid ranges to cover all (alpha, beta) combinations
	# Max GI mean and SD occur at max alpha, min beta
	max_gi_mean <- max(alpha_range) / min(beta_range)
	max_gi_sd   <- sqrt(max(alpha_range)) / min(beta_range)

	# tau grid: covers Gamma(alpha, beta) for all (alpha, beta)
	# Start at small positive value to avoid Inf from dgamma(0, shape<1, rate)
	tau_grid <- seq(1e-4, max_gi_mean + 8 * max_gi_sd, length.out = 401)
	dtau <- diff(tau_grid)[1]

	# s grid: covers f_S support for all (alpha, beta)
	s_lo <- -8 * sd_d
	s_hi <- max_gi_mean + 8 * (max_gi_sd + sd_d)
	s_grid <- seq(s_lo, s_hi, length.out = 501)

	# Pre-tabulate f_{d_j - d_0} on a fine grid
	ddiff_extent <- max(abs(s_lo), s_hi) + max(tau_grid)
	ddiff_grid <- seq(-ddiff_extent, ddiff_extent, length.out = 2001)
	ddiff_vals <- dgamma_diff_gi(ddiff_grid, a_obs, b_obs, a_obs, b_obs)
	ddiff_interp <- approxfun(ddiff_grid, ddiff_vals, rule = 2, yleft = 0, yright = 0)

	# Pre-compute D matrix: D[j, k] = f_{d-d0}(s_j - tau_k)
	D_matrix <- outer(s_grid, tau_grid, function(s, t) ddiff_interp(s - t))

	# Return engine
	list(
		D_matrix  = D_matrix,
		tau_grid  = tau_grid,
		s_grid    = s_grid,
		dtau      = dtau,
		eval = function(alpha, beta, obs) {
			# Compute serial interval density for given (alpha, beta)
			g_vec <- dgamma(tau_grid, shape = alpha, rate = beta)
			# Guard against Inf/NaN from dgamma (e.g. shape < 1 near tau=0)
			g_vec[!is.finite(g_vec)] <- 0
			f_s_vals <- as.numeric(D_matrix %*% g_vec) * dtau
			f_s <- approxfun(s_grid, pmax(f_s_vals, 0), rule = 2, yleft = 0, yright = 0)

			# Evaluate log-likelihood
			f_vals <- f_s(obs)
			f_vals[f_vals < .Machine$double.xmin] <- .Machine$double.xmin
			sum(log(f_vals))
		}
	)
}

# ==============================================================================
# 3. 2D grid posterior for (alpha, beta)
# ==============================================================================

#' Compute the iid posterior on a 2D grid over (alpha, beta)
compute_posterior_ab <- function(serial_intervals, alpha_grid, beta_grid, engine) {
	grid <- expand.grid(alpha = alpha_grid, beta = beta_grid)
	grid$loglik <- vapply(seq_len(nrow(grid)), function(idx) {
		engine$eval(grid$alpha[idx], grid$beta[idx], serial_intervals)
	}, numeric(1))

	# Normalize to posterior (flat prior, log-sum-exp)
	# Replace NaN with -Inf
	grid$loglik[is.nan(grid$loglik)] <- -Inf
	finite_ll <- grid$loglik[is.finite(grid$loglik)]
	if (length(finite_ll) == 0) {
		grid$posterior <- 0
		return(grid[, c("alpha", "beta", "loglik", "posterior")])
	}
	max_ll <- max(finite_ll)
	grid$posterior <- exp(grid$loglik - max_ll)
	grid$posterior[!is.finite(grid$posterior)] <- 0
	grid$posterior <- grid$posterior / sum(grid$posterior)

	grid[, c("alpha", "beta", "loglik", "posterior")]
}

#' Extract summary statistics from 2D posterior
posterior_summary_ab <- function(post_df) {
	# Marginal for alpha
	alpha_marg <- aggregate(posterior ~ alpha, data = post_df, FUN = sum)
	alpha_marg$posterior <- alpha_marg$posterior / sum(alpha_marg$posterior)
	alpha_cdf <- cumsum(alpha_marg$posterior)
	alpha_mean <- sum(alpha_marg$alpha * alpha_marg$posterior)
	alpha_ci_lo <- alpha_marg$alpha[which.min(abs(alpha_cdf - 0.025))]
	alpha_ci_hi <- alpha_marg$alpha[which.min(abs(alpha_cdf - 0.975))]

	# Marginal for beta
	beta_marg <- aggregate(posterior ~ beta, data = post_df, FUN = sum)
	beta_marg$posterior <- beta_marg$posterior / sum(beta_marg$posterior)
	beta_cdf <- cumsum(beta_marg$posterior)
	beta_mean <- sum(beta_marg$beta * beta_marg$posterior)
	beta_ci_lo <- beta_marg$beta[which.min(abs(beta_cdf - 0.025))]
	beta_ci_hi <- beta_marg$beta[which.min(abs(beta_cdf - 0.975))]

	# Posterior mean of GI mean = alpha/beta
	post_df$gi_mean <- post_df$alpha / post_df$beta
	gi_mean_post <- sum(post_df$gi_mean * post_df$posterior)

	# GI mean CI (conservative: from marginal corners)
	gi_mean_ci_lo <- alpha_ci_lo / beta_ci_hi
	gi_mean_ci_hi <- alpha_ci_hi / beta_ci_lo

	list(
		alpha_mean  = alpha_mean,
		alpha_ci_lo = alpha_ci_lo,
		alpha_ci_hi = alpha_ci_hi,
		beta_mean   = beta_mean,
		beta_ci_lo  = beta_ci_lo,
		beta_ci_hi  = beta_ci_hi,
		gi_mean     = gi_mean_post,
		gi_mean_ci_lo = gi_mean_ci_lo,
		gi_mean_ci_hi = gi_mean_ci_hi
	)
}

# ==============================================================================
# 4. ICC and effective sample size
# ==============================================================================

compute_icc_gi <- function(psi, alpha, beta, a_obs, b_obs) {
	var_gi <- alpha / beta^2
	var_d  <- a_obs / b_obs^2
	((1 - psi) * var_gi + var_d) / (var_gi + 2 * var_d)
}

compute_neff <- function(clusters, icc) {
	cluster_sizes <- vapply(clusters, length, integer(1))
	sum(cluster_sizes / (1 + (cluster_sizes - 1) * icc))
}

# ==============================================================================
# 5. Coverage simulation study
# ==============================================================================

cat("--- Coverage simulation study ---\n")

# Parameters: Omicron
pars <- parslist[[2]]
alpha_true <- pars$alpha
beta_true  <- pars$beta
R0         <- pars$R0
a_obs      <- 4
b_obs      <- 1

cat(sprintf("  Pathogen: %s (alpha=%.3f, beta=%.4f, R0=%d)\n",
    pars$pathogen, alpha_true, beta_true, R0))
cat(sprintf("  Detection delays: a_obs=%d, b_obs=%d\n", a_obs, b_obs))

# Grid for (alpha, beta) — centered on true values with wide margins
alpha_grid <- seq(max(0.5, alpha_true * 0.4), alpha_true * 2.2, length.out = 50)
beta_grid  <- seq(max(0.05, beta_true * 0.4), beta_true * 2.2, length.out = 50)

psi_vals     <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0)
K            <- 100
n_replicates <- 200

cache_file <- file.path("output", "g_identifiability_results.csv")

if (file.exists(cache_file)) {
	cat("  Loading cached results from", cache_file, "\n")
	results_df <- read_csv(cache_file, show_col_types = FALSE)
} else {

	# Build the serial density engine (pre-compute D matrix once)
	cat("  Building serial density engine...\n")
	engine <- build_serial_density_engine(a_obs, b_obs, range(alpha_grid), range(beta_grid))
	cat("  Engine ready.\n")

	results_list <- vector("list", length(psi_vals) * n_replicates)
	result_idx <- 0

	for (psi_true in psi_vals) {
		icc <- compute_icc_gi(psi_true, alpha_true, beta_true, a_obs, b_obs)
		cat(sprintf("  psi = %.1f (ICC = %.3f)\n", psi_true, icc))

		for (rep in seq_len(n_replicates)) {
			set.seed(10000 * which(psi_vals == psi_true) + rep)

			clusters <- simulate_clusters(
				K, psi_true, R0, alpha_true, beta_true, a_obs, b_obs, p_asc = 1
			)

			all_si <- unlist(clusters)
			n_obs <- length(all_si)
			n_clusters <- length(clusters)

			if (n_obs < 5) {
				result_idx <- result_idx + 1
				results_list[[result_idx]] <- data.frame(
					psi_true = psi_true, replicate = rep,
					n_obs = n_obs, n_clusters = n_clusters,
					icc = icc, neff = NA,
					alpha_mean = NA, alpha_ci_lo = NA, alpha_ci_hi = NA,
					alpha_covers = NA,
					beta_mean = NA, beta_ci_lo = NA, beta_ci_hi = NA,
					beta_covers = NA,
					gi_mean = NA, gi_mean_ci_lo = NA, gi_mean_ci_hi = NA,
					gi_mean_covers = NA
				)
				next
			}

			neff <- compute_neff(clusters, icc)

			post <- compute_posterior_ab(all_si, alpha_grid, beta_grid, engine)
			summ <- posterior_summary_ab(post)

			result_idx <- result_idx + 1
			results_list[[result_idx]] <- data.frame(
				psi_true    = psi_true,
				replicate   = rep,
				n_obs       = n_obs,
				n_clusters  = n_clusters,
				icc         = icc,
				neff        = neff,
				alpha_mean  = summ$alpha_mean,
				alpha_ci_lo = summ$alpha_ci_lo,
				alpha_ci_hi = summ$alpha_ci_hi,
				alpha_covers = (alpha_true >= summ$alpha_ci_lo & alpha_true <= summ$alpha_ci_hi),
				beta_mean   = summ$beta_mean,
				beta_ci_lo  = summ$beta_ci_lo,
				beta_ci_hi  = summ$beta_ci_hi,
				beta_covers = (beta_true >= summ$beta_ci_lo & beta_true <= summ$beta_ci_hi),
				gi_mean     = summ$gi_mean,
				gi_mean_ci_lo = summ$gi_mean_ci_lo,
				gi_mean_ci_hi = summ$gi_mean_ci_hi,
				gi_mean_covers = (alpha_true / beta_true >= summ$gi_mean_ci_lo &
				                  alpha_true / beta_true <= summ$gi_mean_ci_hi)
			)

			if (rep %% 50 == 0) {
				cat(sprintf("    psi=%.1f, rep %d/%d done\n", psi_true, rep, n_replicates))
			}
		}
	}

	results_df <- bind_rows(results_list)
	write_csv(results_df, cache_file)
	cat(sprintf("  Saved %d results to %s\n", nrow(results_df), cache_file))
}

# ==============================================================================
# 6. Figures
# ==============================================================================

cat("--- Generating figures ---\n")

# --------------------------------------------------------------------------
# Figure 1: Example 2D posteriors at psi = 0.1 vs psi = 1.0
# --------------------------------------------------------------------------

cat("  Figure 1: Example posteriors\n")

# Build engine if not already available
if (!exists("engine")) {
	engine <- build_serial_density_engine(a_obs, b_obs, range(alpha_grid), range(beta_grid))
}

example_psi <- c(0.1, 1.0)
example_posteriors <- list()

for (psi_ex in example_psi) {
	set.seed(2)
	clusters_ex <- simulate_clusters(
		K, psi_ex, R0, alpha_true, beta_true, a_obs, b_obs, p_asc = 1
	)
	all_si_ex <- unlist(clusters_ex)

	post_ex <- compute_posterior_ab(all_si_ex, alpha_grid, beta_grid, engine)
	post_ex$psi <- psi_ex
	example_posteriors <- c(example_posteriors, list(post_ex))
}

post_ex_df <- bind_rows(example_posteriors)
post_ex_df$psi_label <- factor(
	post_ex_df$psi,
	levels = c(0.1, 1.0),
	labels = c("psi = 0.1 (punctuated)", "psi = 1.0 (smooth, iid)")
)

fig1 <- ggplot(post_ex_df, aes(x = alpha, y = beta)) +
	geom_tile(aes(fill = posterior)) +
	geom_contour(aes(z = posterior), color = "white", alpha = 0.5, bins = 8) +
	geom_point(data = data.frame(alpha = alpha_true, beta = beta_true),
	           aes(x = alpha, y = beta),
	           color = "red", size = 3, shape = 4, stroke = 1.5) +
	facet_wrap(~ psi_label) +
	scale_fill_viridis_c(option = "magma") +
	labs(
		x = expression(alpha ~ "(GI shape)"),
		y = expression(beta ~ "(GI rate)"),
		fill = "Posterior",
		title = sprintf("IID posterior for GI parameters (%s, K = %d clusters)",
		                pars$pathogen, K),
		subtitle = "Red cross = true value. Posteriors look similar despite psi = 0.1 data being less informative"
	) +
	theme_minimal(base_size = 13)
save_fig(fig1, "g_identifiability_posteriors", width = 12, height = 5)

# --------------------------------------------------------------------------
# Figure 2: Coverage vs psi (main result)
# --------------------------------------------------------------------------

cat("  Figure 2: Coverage vs psi\n")
coverage_df <- results_df %>%
	filter(!is.na(alpha_covers)) %>%
	group_by(psi_true) %>%
	summarise(
		alpha_coverage = mean(alpha_covers),
		beta_coverage  = mean(beta_covers),
		gi_mean_coverage = mean(gi_mean_covers, na.rm = TRUE),
		n_reps = n(),
		mean_nobs = mean(n_obs),
		mean_neff = mean(neff, na.rm = TRUE),
		mean_icc  = mean(icc),
		.groups = "drop"
	)

coverage_long <- coverage_df %>%
	pivot_longer(
		cols = c(alpha_coverage, beta_coverage, gi_mean_coverage),
		names_to = "parameter",
		values_to = "coverage"
	) %>%
	mutate(parameter = case_when(
		parameter == "alpha_coverage"    ~ "alpha (GI shape)",
		parameter == "beta_coverage"     ~ "beta (GI rate)",
		parameter == "gi_mean_coverage"  ~ "alpha/beta (GI mean)"
	))

coverage_long <- coverage_long %>%
	mutate(se = sqrt(coverage * (1 - coverage) / n_reps))

fig2 <- ggplot(coverage_long, aes(x = psi_true, y = coverage, color = parameter)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2.5) +
	geom_errorbar(aes(ymin = coverage - 1.96 * se, ymax = pmin(1, coverage + 1.96 * se)),
	              width = 0.03, linewidth = 0.5) +
	geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
	scale_x_continuous(breaks = psi_vals) +
	ylim(0.5, 1.0) +
	labs(
		x = expression(psi),
		y = "Empirical coverage of 95% CI",
		color = "Parameter",
		title = sprintf("IID posterior coverage vs psi (%s, K = %d, %d replicates)",
		                pars$pathogen, K, n_replicates),
		subtitle = "Dashed line = nominal 95%. Coverage drops at low psi due to ignored within-cluster correlation"
	) +
	theme_minimal(base_size = 13) +
	theme(legend.position = "right")
save_fig(fig2, "g_identifiability_coverage", width = 10, height = 6)

# --------------------------------------------------------------------------
# Figure 3: Effective sample size comparison
# --------------------------------------------------------------------------

cat("  Figure 3: Effective sample size\n")
neff_summary <- results_df %>%
	filter(!is.na(neff)) %>%
	group_by(psi_true) %>%
	summarise(
		mean_nobs = mean(n_obs),
		mean_neff = mean(neff),
		neff_ratio = mean(neff) / mean(n_obs),
		alpha_coverage = mean(alpha_covers, na.rm = TRUE),
		.groups = "drop"
	)

# Theoretical ICC curve
psi_fine <- seq(0.05, 1, length.out = 100)
icc_fine <- sapply(psi_fine, compute_icc_gi,
                   alpha = alpha_true, beta = beta_true, a_obs = a_obs, b_obs = b_obs)

# Average cluster size
avg_m <- mean(results_df$n_obs / results_df$n_clusters, na.rm = TRUE)
neff_ratio_theory <- 1 / (1 + (avg_m - 1) * icc_fine)

fig3 <- ggplot() +
	geom_line(data = data.frame(psi = psi_fine, ratio = neff_ratio_theory),
	          aes(x = psi, y = ratio), linewidth = 0.8, color = "steelblue") +
	geom_point(data = neff_summary,
	           aes(x = psi_true, y = neff_ratio),
	           color = "darkred", size = 3) +
	geom_point(data = neff_summary,
	           aes(x = psi_true, y = (alpha_coverage - 0.5) / 0.5),
	           color = "orange", size = 2.5, shape = 17) +
	scale_y_continuous(
		name = expression(N[eff] / N),
		sec.axis = sec_axis(~ . * 0.5 + 0.5, name = "alpha coverage")
	) +
	labs(
		x = expression(psi),
		title = sprintf("Effective sample size and coverage (%s)", pars$pathogen),
		subtitle = "Blue = theoretical N_eff/N, red dots = observed, orange triangles = alpha coverage"
	) +
	theme_minimal(base_size = 13)
save_fig(fig3, "g_identifiability_neff", width = 9, height = 6)

# --------------------------------------------------------------------------
# Print summary table
# --------------------------------------------------------------------------

cat("\n--- Coverage summary ---\n")
cat(sprintf("  %-6s  %-8s  %-8s  %-10s  %-10s  %-10s  %-10s\n",
    "psi", "ICC", "N_eff/N", "α cover", "β cover", "mean cover", "n_obs"))
for (i in seq_len(nrow(coverage_df))) {
	r <- coverage_df[i, ]
	cat(sprintf("  %-6.1f  %-8.3f  %-8.3f  %-10.3f  %-10.3f  %-10.3f  %-10.1f\n",
	    r$psi_true, r$mean_icc,
	    r$mean_neff / r$mean_nobs,
	    r$alpha_coverage, r$beta_coverage, r$gi_mean_coverage,
	    r$mean_nobs))
}

# ==============================================================================
# 7. MI-based effective sample size
# ==============================================================================
#
# Compute the actual mutual information between GI parameters and clustered
# data, using the cluster-aware likelihood (no Gaussian approximation).
# Compare with iid MI to get the MI-based n_eff for alpha and beta separately.
# ==============================================================================

#' Build engine for density of Gamma(shape, rate) + Gamma(a_obs, b_obs)
#' Precomputes D[j,k] = dgamma(x_j - t_k, a_obs, b_obs), then for each
#' (shape, rate), f_+(x) = dt * D %*% dgamma(t, shape, rate).
build_sum_engine <- function(a_obs, b_obs, shape_range, rate_range) {
	sd_d <- sqrt(a_obs) / b_obs
	max_mean <- max(shape_range) / min(rate_range) + a_obs / b_obs
	max_sd <- sqrt(max(shape_range)) / min(rate_range) + sd_d

	t_grid <- seq(1e-4, max(shape_range) / min(rate_range) + 6 * sqrt(max(shape_range)) / min(rate_range),
	              length.out = 301)
	dt <- diff(t_grid)[1]
	x_grid <- seq(0, max_mean + 8 * max_sd, length.out = 401)

	D <- outer(x_grid, t_grid, function(x, t) dgamma(x - t, shape = a_obs, rate = b_obs))

	list(
		t_grid = t_grid, x_grid = x_grid, dt = dt, D = D,
		eval = function(shape, rate) {
			g <- dgamma(t_grid, shape = shape, rate = rate)
			g[!is.finite(g)] <- 0
			f_vals <- as.numeric(D %*% g) * dt
			approxfun(x_grid, pmax(f_vals, 0), rule = 2, yleft = 0, yright = 0)
		}
	)
}

#' Compute cluster-aware posterior on (alpha, beta) grid
#' Uses the decomposition s_j = delta + w_j where delta = l - d0 (shared)
#' and w_j = eps_j + d_j (independent per sibling).
#' Vectorized: evaluates f_plus on all serial intervals × delta grid at once.
compute_cluster_posterior <- function(clusters, alpha_grid, beta_grid,
                                      a_obs, b_obs, psi, sum_eng) {
	grid <- expand.grid(alpha = alpha_grid, beta = beta_grid)

	# Filter to non-empty clusters
	clusters <- clusters[sapply(clusters, length) > 0]
	K <- length(clusters)
	if (K == 0) {
		grid$loglik <- 0
		grid$posterior <- 1 / nrow(grid)
		return(grid[, c("alpha", "beta", "loglik", "posterior")])
	}

	# Pre-compute cluster structure (shared across all grid points)
	all_si <- unlist(clusters)
	cluster_idx <- rep(seq_len(K), sapply(clusters, length))
	N_total <- length(all_si)

	# Delta grid for integrating out shift = l - d0
	d0_mean <- a_obs / b_obs
	d0_sd <- sqrt(a_obs) / b_obs
	if (psi < 1 - 1e-10) {
		max_l <- max(alpha_grid) * (1 - psi) / min(beta_grid)
		max_l_sd <- sqrt(max(alpha_grid) * (1 - psi)) / min(beta_grid)
		delta_hi <- max_l + 6 * max_l_sd
	} else {
		delta_hi <- 1
	}
	delta_lo <- -(d0_mean + 6 * d0_sd)
	n_delta <- 151L
	delta_grid <- seq(delta_lo, delta_hi, length.out = n_delta)
	d_delta <- diff(delta_grid)[1]

	# Pre-compute the big matrix: all_si_minus_delta[i, k] = all_si[i] - delta_grid[k]
	# Dimensions: N_total × n_delta
	si_minus_delta <- outer(all_si, delta_grid, "-")

	# For f_shift: y grid and weights
	y_grid <- seq(1e-4, d0_mean + 8 * d0_sd, length.out = 151)
	dy <- diff(y_grid)[1]
	h_vec <- dgamma(y_grid, shape = a_obs, rate = b_obs)
	delta_plus_y <- outer(delta_grid, y_grid, "+")

	grid$loglik <- vapply(seq_len(nrow(grid)), function(idx) {
		a <- grid$alpha[idx]
		b <- grid$beta[idx]

		# f_plus: density of Gamma(psi*a, b) + Gamma(a_obs, b_obs)
		f_plus <- sum_eng$eval(psi * a, b)

		# f_shift: density of Gamma((1-psi)*a, b) - Gamma(a_obs, b_obs)
		if (psi < 1 - 1e-10) {
			G_shift <- dgamma(delta_plus_y, shape = (1 - psi) * a, rate = b)
			f_shift_vals <- dy * as.numeric(G_shift %*% h_vec)
		} else {
			f_shift_vals <- dgamma(-delta_grid, shape = a_obs, rate = b_obs)
		}
		f_shift_vals[f_shift_vals < .Machine$double.xmin] <- .Machine$double.xmin
		log_f_shift <- log(f_shift_vals)

		# Evaluate f_plus on ALL shifted values at once (vectorized)
		fv_all <- f_plus(as.vector(si_minus_delta))
		fv_all[fv_all < .Machine$double.xmin] <- .Machine$double.xmin
		log_fv <- matrix(log(fv_all), nrow = N_total, ncol = n_delta)

		# Sum log f_plus per cluster using rowsum (fast grouped sum)
		log_f_sums <- rowsum(log_fv, cluster_idx)  # K × n_delta

		# Add log f_shift to each row
		log_integrand <- sweep(log_f_sums, 2, log_f_shift, "+")

		# Log-sum-exp per cluster (vectorized)
		row_max <- apply(log_integrand, 1, function(x) {
			fx <- x[is.finite(x)]
			if (length(fx) == 0) -Inf else max(fx)
		})
		if (any(!is.finite(row_max))) return(-Inf)
		log_Lc <- row_max + log(rowSums(exp(log_integrand - row_max))) + log(d_delta)

		sum(log_Lc)
	}, numeric(1))

	# Normalize
	grid$loglik[is.nan(grid$loglik)] <- -Inf
	finite_ll <- grid$loglik[is.finite(grid$loglik)]
	if (length(finite_ll) == 0) {
		grid$posterior <- 0
	} else {
		max_ll <- max(finite_ll)
		grid$posterior <- exp(grid$loglik - max_ll)
		grid$posterior[!is.finite(grid$posterior)] <- 0
		grid$posterior <- grid$posterior / sum(grid$posterior)
	}
	grid[, c("alpha", "beta", "loglik", "posterior")]
}


#' Compute marginal posterior standard deviations from 2D grid posterior
posterior_sds <- function(post_df) {
	am <- aggregate(posterior ~ alpha, data = post_df, FUN = sum)
	am$posterior <- am$posterior / sum(am$posterior)
	E_a <- sum(am$alpha * am$posterior)
	sd_a <- sqrt(max(sum(am$alpha^2 * am$posterior) - E_a^2, 0))

	bm <- aggregate(posterior ~ beta, data = post_df, FUN = sum)
	bm$posterior <- bm$posterior / sum(bm$posterior)
	E_b <- sum(bm$beta * bm$posterior)
	sd_b <- sqrt(max(sum(bm$beta^2 * bm$posterior) - E_b^2, 0))

	c(sd_alpha = sd_a, sd_beta = sd_b)
}

#' Create a closure for the cluster-aware log-likelihood
#'
#' Precomputes shared quantities (delta grid, si_minus_delta, etc.) once,
#' then returns a function loglik(c(alpha, beta)) for use with optim.
make_cluster_loglik_fn <- function(clusters, a_obs, b_obs, psi, sum_eng,
                                    alpha_range = c(0.5, 5), beta_range = c(0.1, 1)) {
	clusters <- clusters[sapply(clusters, length) > 0]
	K <- length(clusters)
	if (K == 0) return(function(par) 0)

	all_si <- unlist(clusters)
	cluster_idx <- rep(seq_len(K), sapply(clusters, length))
	N_total <- length(all_si)

	d0_mean <- a_obs / b_obs
	d0_sd <- sqrt(a_obs) / b_obs
	if (psi < 1 - 1e-10) {
		max_l <- max(alpha_range) * (1 - psi) / min(beta_range)
		max_l_sd <- sqrt(max(alpha_range) * (1 - psi)) / min(beta_range)
		delta_hi <- max_l + 6 * max_l_sd
	} else {
		delta_hi <- 1
	}
	delta_lo <- -(d0_mean + 6 * d0_sd)
	n_delta <- 151L
	delta_grid <- seq(delta_lo, delta_hi, length.out = n_delta)
	d_delta <- diff(delta_grid)[1]

	si_minus_delta <- outer(all_si, delta_grid, "-")

	y_grid <- seq(1e-4, d0_mean + 8 * d0_sd, length.out = 151)
	dy <- diff(y_grid)[1]
	h_vec <- dgamma(y_grid, shape = a_obs, rate = b_obs)
	delta_plus_y <- outer(delta_grid, y_grid, "+")

	function(par) {
		a <- par[1]; b <- par[2]
		if (a <= 0.01 || b <= 0.01) return(-1e20)

		f_plus <- tryCatch(sum_eng$eval(psi * a, b), error = function(e) NULL)
		if (is.null(f_plus)) return(-1e20)

		if (psi < 1 - 1e-10) {
			shape_shift <- (1 - psi) * a
			if (shape_shift < 0.001) return(-1e20)
			G_shift <- dgamma(delta_plus_y, shape = shape_shift, rate = b)
			f_shift_vals <- dy * as.numeric(G_shift %*% h_vec)
		} else {
			f_shift_vals <- dgamma(-delta_grid, shape = a_obs, rate = b_obs)
		}
		f_shift_vals[f_shift_vals < .Machine$double.xmin] <- .Machine$double.xmin
		log_f_shift <- log(f_shift_vals)

		fv_all <- f_plus(as.vector(si_minus_delta))
		fv_all[fv_all < .Machine$double.xmin] <- .Machine$double.xmin
		log_fv <- matrix(log(fv_all), nrow = N_total, ncol = n_delta)

		log_f_sums <- rowsum(log_fv, cluster_idx)
		log_integrand <- sweep(log_f_sums, 2, log_f_shift, "+")

		row_max <- apply(log_integrand, 1, function(x) {
			fx <- x[is.finite(x)]
			if (length(fx) == 0) -Inf else max(fx)
		})
		if (any(!is.finite(row_max))) return(-1e20)
		log_Lc <- row_max + log(rowSums(exp(log_integrand - row_max))) + log(d_delta)

		sum(log_Lc)
	}
}

#' Compute posterior variance via Laplace approximation (optim + Hessian)
posterior_var_laplace <- function(neg_loglik_fn, start = c(2.39, 0.339),
                                   lower = c(0.5, 0.05), upper = c(6, 1.5)) {
	result <- tryCatch(
		optim(start, neg_loglik_fn, method = "L-BFGS-B",
		      lower = lower, upper = upper, hessian = TRUE,
		      control = list(factr = 1e7)),
		error = function(e) NULL
	)

	if (is.null(result) || result$convergence != 0)
		return(c(var_alpha = NA, var_beta = NA, mode_alpha = NA, mode_beta = NA))

	# Hessian of neg log-lik = observed Fisher information
	# Posterior variance approx inv(Fisher info)
	V <- tryCatch(solve(result$hessian), error = function(e) matrix(NA, 2, 2))

	c(var_alpha = V[1,1], var_beta = V[2,2],
	  mode_alpha = result$par[1], mode_beta = result$par[2])
}

# --- Main: effective sample size via Laplace approximation ---

cat("\n--- MI-based effective sample size ---\n")

mi_cache <- file.path("output", "g_mi_neff_results.csv")

if (file.exists(mi_cache)) {
	cat("  Loading cached MI results\n")
	mi_df <- read_csv(mi_cache, show_col_types = FALSE)
} else {
	psi_mi   <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.0)
	n_mc     <- 30
	K_mi     <- 40

	# Serial density engine for iid log-likelihood
	alpha_range_mi <- c(0.5, 5)
	beta_range_mi  <- c(0.1, 1)
	cat("  Building engines for MI computation...\n")
	serial_eng <- build_serial_density_engine(a_obs, b_obs,
	                                          alpha_range_mi, beta_range_mi)
	cat("  Engines ready.\n")

	mi_results <- list()

	for (psi in psi_mi) {
		cat(sprintf("  psi = %.1f: ", psi))

		# Build sum engine for this psi
		psi_shape_range <- c(psi * min(alpha_range_mi), psi * max(alpha_range_mi))
		sum_eng <- build_sum_engine(a_obs, b_obs,
		                            shape_range = psi_shape_range,
		                            rate_range = beta_range_mi)

		var_cl_a <- var_cl_b <- numeric(n_mc)
		var_iid_a <- var_iid_b <- numeric(n_mc)
		mode_iid_a <- mode_iid_b <- numeric(n_mc)
		mode_cl_a <- mode_cl_b <- numeric(n_mc)
		N_total <- numeric(n_mc)

		for (i in seq_len(n_mc)) {
			set.seed(20000 + which(psi_mi == psi) * 1000 + i)
			cls <- simulate_clusters(K_mi, psi, R0, alpha_true, beta_true,
			                         a_obs, b_obs, p_asc = 1)
			all_si <- unlist(cls)
			N_total[i] <- length(all_si)

			if (N_total[i] < 5) {
				var_cl_a[i] <- var_cl_b[i] <- NA
				var_iid_a[i] <- var_iid_b[i] <- NA
				mode_iid_a[i] <- mode_iid_b[i] <- NA
				mode_cl_a[i] <- mode_cl_b[i] <- NA
				next
			}

			# Cluster-aware: Laplace approximation
			ll_fn <- make_cluster_loglik_fn(cls, a_obs, b_obs, psi, sum_eng,
			                                alpha_range = alpha_range_mi,
			                                beta_range = beta_range_mi)
			neg_ll_cl <- function(par) -ll_fn(par)
			res_cl <- posterior_var_laplace(neg_ll_cl)
			var_cl_a[i]  <- res_cl["var_alpha"]
			var_cl_b[i]  <- res_cl["var_beta"]
			mode_cl_a[i] <- res_cl["mode_alpha"]
			mode_cl_b[i] <- res_cl["mode_beta"]

			# IID: Laplace approximation
			all_si_local <- all_si  # capture for closure
			neg_ll_iid <- function(par) -serial_eng$eval(par[1], par[2], all_si_local)
			res_iid <- posterior_var_laplace(neg_ll_iid)
			var_iid_a[i]  <- res_iid["var_alpha"]
			var_iid_b[i]  <- res_iid["var_beta"]
			mode_iid_a[i] <- res_iid["mode_alpha"]
			mode_iid_b[i] <- res_iid["mode_beta"]

			if (i %% 10 == 0) cat(".")
		}
		cat("\n")

		mi_results <- c(mi_results, list(data.frame(
			type = "clustered", psi = psi,
			n = mean(N_total, na.rm = TRUE),
			# Mean posterior variances (Laplace)
			var_cluster_alpha = mean(var_cl_a, na.rm = TRUE),
			var_cluster_beta  = mean(var_cl_b, na.rm = TRUE),
			var_iid_alpha     = mean(var_iid_a, na.rm = TRUE),
			var_iid_beta      = mean(var_iid_b, na.rm = TRUE),
			# Frequentist variance of the iid posterior mode
			var_freq_alpha = var(mode_iid_a, na.rm = TRUE),
			var_freq_beta  = var(mode_iid_b, na.rm = TRUE)
		)))
	}

	mi_df <- bind_rows(mi_results)
	write_csv(mi_df, mi_cache)
	cat(sprintf("  Saved MI results to %s\n", mi_cache))
}

# --- Compute n_eff from Laplace variances ---

cluster_mi <- mi_df %>% filter(type == "clustered")

# n_eff_cluster / N: how many iid observations match the cluster-aware precision
# n_eff = N * Var_iid / Var_cluster  (if cluster is more precise, n_eff > N)
cluster_mi <- cluster_mi %>%
	mutate(
		neff_alpha = n * var_iid_alpha / var_cluster_alpha,
		neff_beta  = n * var_iid_beta / var_cluster_beta,
		# Calibration: how overconfident is the iid model?
		# overcal = Var_freq / Var_iid  (> 1 means overconfident)
		overcal_alpha = var_freq_alpha / var_iid_alpha,
		overcal_beta  = var_freq_beta / var_iid_beta
	)

# --- Print summary ---
cat("\n--- MI-based n_eff summary ---\n")
cat(sprintf("  %-5s  %5s  %9s  %9s  %8s  %9s  %9s  %9s  %9s\n",
    "psi", "N", "SD_cl(a)", "SD_id(a)", "Kish",
    "neff_a/N", "neff_b/N", "ocal_a", "ocal_b"))
for (i in seq_len(nrow(cluster_mi))) {
	r <- cluster_mi[i, ]
	icc <- compute_icc_gi(r$psi, alpha_true, beta_true, a_obs, b_obs)
	avg_m <- r$n / K_mi
	kish <- 1 / (1 + (avg_m - 1) * icc)
	cat(sprintf("  %-5.1f  %5.0f  %9.4f  %9.4f  %8.3f  %9.3f  %9.3f  %9.2f  %9.2f\n",
	    r$psi, r$n,
	    sqrt(r$var_cluster_alpha), sqrt(r$var_iid_alpha),
	    kish, r$neff_alpha / r$n, r$neff_beta / r$n,
	    r$overcal_alpha, r$overcal_beta))
}

# --- Figure 4: n_eff comparison ---
cat("  Figure 4: MI-based effective sample size\n")

psi_fine <- seq(0.05, 1, length.out = 100)
icc_fine <- sapply(psi_fine, compute_icc_gi,
                   alpha = alpha_true, beta = beta_true, a_obs = a_obs, b_obs = b_obs)
avg_m_all <- mean(cluster_mi$n) / K_mi
kish_curve <- data.frame(psi = psi_fine, ratio = 1 / (1 + (avg_m_all - 1) * icc_fine))

mi_plot_df <- cluster_mi %>%
	mutate(neff_alpha_ratio = neff_alpha / n,
	       neff_beta_ratio  = neff_beta / n) %>%
	select(psi, neff_alpha_ratio, neff_beta_ratio) %>%
	pivot_longer(cols = c(neff_alpha_ratio, neff_beta_ratio),
	             names_to = "parameter", values_to = "neff_ratio") %>%
	mutate(parameter = ifelse(parameter == "neff_alpha_ratio",
	                          "alpha (Laplace)", "beta (Laplace)"))

fig4 <- ggplot() +
	geom_line(data = kish_curve, aes(x = psi, y = ratio, linetype = "Kish (design effect)"),
	          linewidth = 0.8, color = "grey40") +
	geom_line(data = mi_plot_df,
	          aes(x = psi, y = neff_ratio, color = parameter), linewidth = 0.7) +
	geom_point(data = mi_plot_df,
	           aes(x = psi, y = neff_ratio, color = parameter), size = 3) +
	geom_hline(yintercept = 1, linetype = "dotted", color = "grey60") +
	scale_linetype_manual(values = c("Kish (design effect)" = "dashed")) +
	scale_color_manual(values = c("alpha (Laplace)" = "steelblue",
	                              "beta (Laplace)" = "darkred")) +
	labs(
		x = expression(psi),
		y = expression(N[eff] / N),
		color = "Parameter", linetype = "",
		title = sprintf("Cluster-aware effective sample size (%s, K=%d)", pars$pathogen, K_mi),
		subtitle = "Cluster-aware posterior precision vs iid (Laplace approx). Dotted = N_eff = N."
	) +
	theme_minimal(base_size = 13) +
	theme(legend.position = "right")
save_fig(fig4, "g_mi_neff", width = 10, height = 6)

cat("\n=== GI parameter identifiability complete ===\n")
