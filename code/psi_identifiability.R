# ==============================================================================
# Power analysis for excluding high psi
# ==============================================================================
#
# Question: if psi_true is near 0, how many clusters K do we need so that
# P(psi > 0.5 | data) < 0.05?
#
# Key optimization: precompute density interpolators once per (pathogen, obs_delay)
# combo, then reuse across all (psi_true, K, replicate) iterations.
#
# Depends on: parslist (from parameters.R), save_fig, simulate_clusters (from utils.R)
# ==============================================================================

cat("=== Power analysis for psi exclusion ===\n")

# ==============================================================================
# 1. Embedded helper functions
# ==============================================================================

#' Density of X + Y where X ~ Gamma(shape1, rate1), Y ~ Gamma(shape2, rate2)
dgamma_sum_psi <- function(x, shape1, rate1, shape2, rate2) {
	if (abs(rate1 - rate2) < 1e-10 * max(rate1, rate2)) {
		return(dgamma(x, shape = shape1 + shape2, rate = rate1))
	}
	sapply(x, function(xi) {
		if (xi <= 0) return(0)
		integrand <- function(t) {
			dgamma(t, shape1, rate1) * dgamma(xi - t, shape2, rate2)
		}
		tryCatch(
			integrate(integrand, lower = 0, upper = xi,
			          rel.tol = 1e-8, abs.tol = 1e-12)$value,
			error = function(e) 0
		)
	})
}

#' Density of X - Y where X ~ Gamma(shape1, rate1), Y ~ Gamma(shape2, rate2)
dgamma_diff_psi <- function(d, shape1, rate1, shape2, rate2) {
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

#' Pre-tabulate a density on a grid and return an approxfun interpolator
make_density_interp_psi <- function(dfun, grid, ...) {
	vals <- dfun(grid, ...)
	approxfun(grid, vals, rule = 2, yleft = 0, yright = 0)
}

#' Log-likelihood for one cluster via delta quadrature
loglik_cluster_psi <- function(s_vec, f_nu_interp, f_delta_interp, delta_grid) {
	dd <- diff(delta_grid)[1]
	m <- length(s_vec)

	f_delta_vals <- f_delta_interp(delta_grid)

	log_integrand <- rep(0, length(delta_grid))
	for (j in seq_len(m)) {
		f_nu_vals <- f_nu_interp(s_vec[j] - delta_grid)
		f_nu_vals[f_nu_vals < .Machine$double.xmin] <- .Machine$double.xmin
		log_integrand <- log_integrand + log(f_nu_vals)
	}
	log_integrand <- log_integrand + log(pmax(f_delta_vals, .Machine$double.xmin))

	max_li <- max(log_integrand)
	if (is.infinite(max_li) && max_li < 0) return(-Inf)
	log(sum(exp(log_integrand - max_li))) + max_li + log(dd)
}

#' Intraclass correlation coefficient for sibling serial intervals
compute_icc_psi <- function(psi, alpha, beta, a_obs, b_obs) {
	var_gi <- alpha / beta^2
	var_d  <- a_obs / b_obs^2
	((1 - psi) * var_gi + var_d) / (var_gi + 2 * var_d)
}

# ==============================================================================
# 2. Precompute density engines for all psi values
# ==============================================================================

#' Precompute f_nu and f_delta interpolators for every psi in psi_grid
#'
#' This is the key optimization: for a given (alpha, beta, a_obs, b_obs),
#' the density interpolators depend only on psi. We compute them once and
#' reuse across all replicates.
#'
#' @param psi_grid Vector of psi values
#' @param alpha,beta Generation interval Gamma parameters
#' @param a_obs,b_obs Detection delay Gamma parameters
#' @return List with: f_nu_interps (list of approxfuns), f_delta_interps (list),
#'         delta_grid (numeric vector)
precompute_density_engines <- function(psi_grid, alpha, beta, a_obs, b_obs) {
	mean_gi <- alpha / beta
	sd_gi   <- sqrt(alpha) / beta
	mean_d  <- a_obs / b_obs
	sd_d    <- sqrt(a_obs) / b_obs

	# Grids for quadrature (same for all psi)
	delta_lo <- -(mean_d + 5 * sd_d)
	delta_hi <- mean_gi + 5 * sd_gi
	delta_grid <- seq(delta_lo, delta_hi, length.out = 401)

	nu_hi <- mean_gi + mean_d + 5 * (sd_gi + sd_d)
	nu_grid <- seq(1e-6, nu_hi, length.out = 501)

	n_cores <- max(1L, parallel::detectCores(logical = FALSE))

	cat(sprintf("    Precomputing density engines for %d psi values on %d cores...\n",
	    length(psi_grid), n_cores))

	# Parallel computation over psi values
	engines <- parallel::mclapply(psi_grid, function(psi) {
		# psi = 0: eps is point mass at 0, so nu = d_j ~ Gamma(a_obs, b_obs)
		if (psi < 1e-10) {
			f_nu <- approxfun(nu_grid, dgamma(nu_grid, shape = a_obs, rate = b_obs),
			                  rule = 2, yleft = 0, yright = 0)
			f_delta <- make_density_interp_psi(
				dgamma_diff_psi, delta_grid,
				shape1 = alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
		# psi = 1: l is point mass at 0, so delta = -d_0, f_delta(d) = f_{d_0}(-d)
		} else if (psi > 1 - 1e-10) {
			f_nu <- make_density_interp_psi(
				dgamma_sum_psi, nu_grid,
				shape1 = alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
			f_delta <- approxfun(delta_grid, dgamma(-delta_grid, shape = a_obs, rate = b_obs),
			                     rule = 2, yleft = 0, yright = 0)
		} else {
			f_nu <- make_density_interp_psi(
				dgamma_sum_psi, nu_grid,
				shape1 = psi * alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
			f_delta <- make_density_interp_psi(
				dgamma_diff_psi, delta_grid,
				shape1 = (1 - psi) * alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
		}

		list(f_nu = f_nu, f_delta = f_delta)
	}, mc.cores = n_cores)

	cat("    Engines ready.\n")

	list(
		f_nu_interps    = lapply(engines, `[[`, "f_nu"),
		f_delta_interps = lapply(engines, `[[`, "f_delta"),
		delta_grid      = delta_grid
	)
}

# ==============================================================================
# 3. Compute posterior over psi given data
# ==============================================================================

#' Compute posterior weights over psi_grid using precomputed engines
#'
#' @param clusters List of numeric vectors (serial intervals per cluster, m >= 2)
#' @param engines Output of precompute_density_engines()
#' @param psi_grid The psi grid corresponding to engines
#' @return Numeric vector of posterior weights (sums to 1) over psi_grid
compute_psi_posterior <- function(clusters, engines, psi_grid) {
	f_nu_interps    <- engines$f_nu_interps
	f_delta_interps <- engines$f_delta_interps
	delta_grid      <- engines$delta_grid

	# Compute log-likelihood for each psi value using cached interpolators
	logliks <- sapply(seq_along(psi_grid), function(idx) {
		total <- 0
		for (k in seq_along(clusters)) {
			ll <- loglik_cluster_psi(
				clusters[[k]],
				f_nu_interps[[idx]],
				f_delta_interps[[idx]],
				delta_grid
			)
			total <- total + ll
		}
		total
	})

	# Normalize to posterior (flat prior)
	max_ll <- max(logliks)
	log_post <- logliks - max_ll
	post <- exp(log_post)
	post / sum(post)  # discrete probability weights
}

#' Compute P(psi > threshold | data) — wrapper for backward compatibility
compute_tail_prob <- function(clusters, engines, psi_grid, threshold) {
	post <- compute_psi_posterior(clusters, engines, psi_grid)
	sum(post[psi_grid > threshold])
}

#' Compute P(lo <= psi <= hi | data)
compute_interval_prob <- function(clusters, engines, psi_grid, lo, hi) {
	post <- compute_psi_posterior(clusters, engines, psi_grid)
	sum(post[psi_grid >= lo & psi_grid <= hi])
}

# ==============================================================================
# 4. Simulation parameters
# ==============================================================================

# obs_delay_scenarios <- list(
# 	list(label = "fast",     a_obs = 4, b_obs = 2),    # mean 2d, var 1
# 	list(label = "moderate", a_obs = 4, b_obs = 1),    # mean 4d, var 4
# 	list(label = "slow",     a_obs = 4, b_obs = 0.5)   # mean 8d, var 16
# )
obs_delay_scenarios <- list(
	list(label = "fast",     a_obs = 4, b_obs = 4),    # mean 1d, var 1/4
	list(label = "moderate", a_obs = 4, b_obs = 2),    # mean 2d, var 1
	list(label = "slow",     a_obs = 4, b_obs = 1)     # mean 4d, var 4
)

psi_true_vals <- c(0, 0.5, 1)
# K_vals        <- c(25, 50, 100, 200, 500, 1000)
K_vals        <- c(10, 25, 50, 100)
n_replicates  <- 200
psi_grid      <- seq(0, 1, length.out = 101)
# For psi=0: exclude psi > 0.2 (tail test)
# For psi=1: exclude psi < 0.8 (tail test)
# For psi=0.5: concentrate in [0.3, 0.7] (interval test)
# All tests use the same tolerance: psi within 0.2 of true value
psi_thresholds <- list(
	"0"   = list(type = "upper_tail", bound = 0.2),
	"0.5" = list(type = "interval",   lo = 0.3, hi = 0.7),
	"1"   = list(type = "lower_tail", bound = 0.8)
)
tail_cutoff    <- 0.05  # "exclude" if relevant tail/interval prob < this
power_target  <- 0.80
p_asc         <- 1.0

cache_file <- file.path("output", "psi_identifiability_results.csv")

# ==============================================================================
# 5. Simulation study
# ==============================================================================

if (file.exists(cache_file)) {
	cat("  Loading cached results from", cache_file, "\n")
	results_df <- read_csv(cache_file, show_col_types = FALSE)
} else {

	n_cores <- max(1L, parallel::detectCores(logical = FALSE))
	cat(sprintf("  Using %d cores for replicate parallelization\n", n_cores))

	results_list <- list()
	result_idx   <- 0

	for (p in parslist) {
		alpha <- p$alpha
		beta  <- p$beta
		R0    <- p$R0

		for (obs in obs_delay_scenarios) {
			a_obs <- obs$a_obs
			b_obs <- obs$b_obs

			cat(sprintf("\n  %s + %s detection (a=%.0f, b=%.1f):\n",
			    p$pathogen, obs$label, a_obs, b_obs))

			# Precompute density engines ONCE for this (pathogen, obs_delay) combo
			engines <- precompute_density_engines(psi_grid, alpha, beta, a_obs, b_obs)

			for (psi_true in psi_true_vals) {
				for (K in K_vals) {
					test_spec <- psi_thresholds[[as.character(psi_true)]]
					cat(sprintf("    psi_true=%.2f, K=%d: ", psi_true, K))

					# Precompute seed values for reproducibility
					path_idx <- which(sapply(parslist, `[[`, "pathogen") == p$pathogen)
					obs_idx  <- which(sapply(obs_delay_scenarios, `[[`, "label") == obs$label)
					psi_idx  <- which(psi_true_vals == psi_true)
					K_idx    <- which(K_vals == K)
					seed_base <- 100000 * path_idx + 10000 * obs_idx +
					             1000 * psi_idx + 10 * K_idx

					# Parallel replicates
					rep_results <- parallel::mclapply(seq_len(n_replicates), function(rep) {
						set.seed(seed_base + rep)

						clusters <- simulate_clusters(
							K, psi_true, R0, alpha, beta, a_obs, b_obs, p_asc
						)
						n_informative <- length(clusters)

						if (n_informative < 2) {
							return(data.frame(
								pathogen      = p$pathogen,
								R0            = R0,
								obs_delay     = obs$label,
								psi_true      = psi_true,
								K             = K,
								replicate     = rep,
								n_informative = n_informative,
								target_prob   = NA
							))
						}

						# Compute posterior and extract relevant probability
						post <- compute_psi_posterior(clusters, engines, psi_grid)
						target_p <- if (test_spec$type == "interval") {
							sum(post[psi_grid >= test_spec$lo & psi_grid <= test_spec$hi])
						} else if (test_spec$type == "upper_tail") {
							1 - sum(post[psi_grid > test_spec$bound])
						} else {  # lower_tail
							1 - sum(post[psi_grid < test_spec$bound])
						}

						data.frame(
							pathogen      = p$pathogen,
							R0            = R0,
							obs_delay     = obs$label,
							psi_true      = psi_true,
							K             = K,
							replicate     = rep,
							n_informative = n_informative,
							target_prob   = target_p
						)
					}, mc.cores = n_cores)

					rep_df <- bind_rows(rep_results)
					for (i in seq_len(nrow(rep_df))) {
						result_idx <- result_idx + 1
						results_list[[result_idx]] <- rep_df[i, ]
					}

					# Progress: show fraction where target region has > 95% mass
					frac_excluded <- mean(rep_df$target_prob > (1 - tail_cutoff), na.rm = TRUE)
					cat(sprintf("power = %.2f\n", frac_excluded))
				}
			}
		}
	}

	results_df <- bind_rows(results_list)
	write_csv(results_df, cache_file)
	cat(sprintf("\n  Saved %d results to %s\n", nrow(results_df), cache_file))
}

# ==============================================================================
# 6. Compute power summaries
# ==============================================================================

results_df <- results_df %>%
	mutate(obs_delay = factor(obs_delay, levels = c("fast", "moderate", "slow")))

# target_prob = P(target region | data)
# For all psi_true values: "excluded" (i.e. identified) when target_prob > 0.95
results_df <- results_df %>%
	mutate(excluded = target_prob > (1 - tail_cutoff))

power_df <- results_df %>%
	filter(!is.na(excluded)) %>%
	group_by(pathogen, R0, obs_delay, psi_true, K) %>%
	summarise(
		power       = mean(excluded),
		n_reps      = n(),
		power_se    = sqrt(power * (1 - power) / n_reps),
		mean_n_inf  = mean(n_informative),
		.groups     = "drop"
	)

pathogen_order <- sapply(parslist, function(p) p$pathogen)
power_df$pathogen <- factor(power_df$pathogen, levels = pathogen_order)

# ==============================================================================
# 7. Figures
# ==============================================================================

cat("--- Generating figures ---\n")

# --------------------------------------------------------------------------
# Figure 1: Power curves — P(exclude psi > 0.5) vs K
# --------------------------------------------------------------------------

cat("  Figure 1: Power curves\n")

obs_labels <- c(fast = "Fast detection (mean 1d)",
                moderate = "Moderate detection (mean 2d)",
                slow = "Slow detection (mean 4d)")

psi_true_labels <- c("0"   = "psi_true = 0 (exclude psi > 0.2)",
                     "0.5" = "psi_true = 0.5 (psi in [0.3, 0.7])",
                     "1"   = "psi_true = 1 (exclude psi < 0.8)")

fig1 <- ggplot(power_df,
               aes(x = K, y = power, color = pathogen)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	geom_errorbar(aes(ymin = pmax(0, power - 1.96 * power_se),
	                  ymax = pmin(1, power + 1.96 * power_se)),
	              width = 0.05, linewidth = 0.4) +
	geom_hline(yintercept = power_target, linetype = "dashed", color = "grey50") +
	scale_x_log10(breaks = K_vals, labels = K_vals) +
	facet_grid(psi_true ~ obs_delay,
	           labeller = labeller(obs_delay = obs_labels,
	                               psi_true = psi_true_labels)) +
	ylim(0, 1) +
	labs(
		x = "Number of index cases (K)",
		y = "Power to identify psi",
		color = "Pathogen",
		title = expression(paste("Power to identify psi (", p[asc], " = 1)"))
	) +
	theme_minimal(base_size = 13) +
	theme(legend.position = "bottom",
	      legend.box = "horizontal")
save_fig(fig1, "psi_identifiability_power", width = 14, height = 10)

# --------------------------------------------------------------------------
# Figure 2: Required K for 80% power (interpolated)
# --------------------------------------------------------------------------

cat("  Figure 2: Required K heatmap\n")

# Interpolate K at which power crosses 80%
required_K <- power_df %>%
	group_by(pathogen, obs_delay, psi_true) %>%
	summarise(
		K_80 = {
			if (max(power) < power_target) {
				Inf  # never reached
			} else if (min(power) >= power_target) {
				min(K)  # already there at smallest K
			} else {
				# Linear interpolation on log(K) scale
				log_K <- log10(K)
				idx_cross <- which(diff(sign(power - power_target)) != 0)[1]
				if (is.na(idx_cross)) {
					Inf
				} else {
					# Interpolate between idx_cross and idx_cross+1
					x1 <- log_K[idx_cross]; x2 <- log_K[idx_cross + 1]
					y1 <- power[idx_cross]; y2 <- power[idx_cross + 1]
					x_interp <- x1 + (power_target - y1) / (y2 - y1) * (x2 - x1)
					round(10^x_interp)
				}
			}
		},
		.groups = "drop"
	)

# Display label
required_K <- required_K %>%
	mutate(K_label = ifelse(is.infinite(K_80), ">1000", as.character(K_80)))

fig2 <- ggplot(required_K %>% filter(is.finite(K_80)),
               aes(x = obs_delay, y = pathogen, fill = log10(K_80))) +
	geom_tile(color = "white", linewidth = 1) +
	geom_text(aes(label = K_label), size = 4.5) +
	geom_text(data = required_K %>% filter(is.infinite(K_80)),
	          aes(x = obs_delay, y = pathogen, label = K_label),
	          size = 4.5, color = "grey50") +
	geom_tile(data = required_K %>% filter(is.infinite(K_80)),
	          aes(x = obs_delay, y = pathogen),
	          fill = "grey90", color = "white", linewidth = 1) +
	geom_text(data = required_K %>% filter(is.infinite(K_80)),
	          aes(x = obs_delay, y = pathogen, label = K_label),
	          size = 4.5, color = "grey50") +
	facet_wrap(~ psi_true, nrow = 1,
	           labeller = labeller(psi_true = psi_true_labels)) +
	scale_fill_viridis_c(option = "plasma", direction = -1,
	                     name = expression(log[10](K))) +
	labs(
		x = "Detection speed",
		y = "Pathogen",
		title = "Required K for 80% power to identify psi"
	) +
	theme_minimal(base_size = 13) +
	theme(legend.position = "right",
	      panel.grid = element_blank())
save_fig(fig2, "psi_identifiability_required_K", width = 12, height = 4.5)

# Print table
cat("\n--- Required K for 80% power ---\n")
cat(sprintf("  %-12s  %-10s  %-10s  K_80\n", "pathogen", "obs_delay", "psi_true"))
for (i in seq_len(nrow(required_K))) {
	r <- required_K[i, ]
	cat(sprintf("  %-12s  %-10s  %-10s  %s\n",
	    r$pathogen, r$obs_delay, r$psi_true, r$K_label))
}

# --------------------------------------------------------------------------
# Figure 3: ICC diagnostic
# --------------------------------------------------------------------------

cat("\n  Figure 3: ICC diagnostic\n")

icc_data <- expand.grid(
	psi = seq(0, 1, length.out = 101),
	pathogen = sapply(parslist, `[[`, "pathogen"),
	obs_delay = sapply(obs_delay_scenarios, `[[`, "label"),
	stringsAsFactors = FALSE
)

icc_data$icc <- NA
for (i in seq_len(nrow(icc_data))) {
	p_idx <- which(sapply(parslist, `[[`, "pathogen") == icc_data$pathogen[i])
	o_idx <- which(sapply(obs_delay_scenarios, `[[`, "label") == icc_data$obs_delay[i])
	pp <- parslist[[p_idx]]
	oo <- obs_delay_scenarios[[o_idx]]
	icc_data$icc[i] <- compute_icc_psi(icc_data$psi[i], pp$alpha, pp$beta, oo$a_obs, oo$b_obs)
}

icc_data$obs_delay <- factor(icc_data$obs_delay, levels = c("fast", "moderate", "slow"))
icc_data$pathogen  <- factor(icc_data$pathogen, levels = pathogen_order)

# Compute ICC range (psi=0 vs psi=1) for annotation
icc_range <- icc_data %>%
	filter(psi %in% c(0, 1)) %>%
	group_by(pathogen, obs_delay) %>%
	summarise(
		icc_at_0 = icc[psi == 0],
		icc_at_1 = icc[psi == 1],
		delta_icc = icc_at_0 - icc_at_1,
		.groups = "drop"
	)

fig3 <- ggplot(icc_data, aes(x = psi, y = icc, color = obs_delay)) +
	geom_line(linewidth = 0.8) +
	facet_wrap(~ pathogen, nrow = 1) +
	scale_color_manual(
		values = c(fast = "#2ca02c", moderate = "#1f77b4", slow = "#d62728"),
		labels = c(fast = "Fast (2d)", moderate = "Moderate (4d)", slow = "Slow (8d)")
	) +
	labs(
		x = expression(psi),
		y = "Intraclass correlation (ICC)",
		color = "Detection delay",
		title = expression(paste("ICC vs ", psi, ": detection noise compresses the signal range")),
		subtitle = "Larger ICC range (from psi=0 to psi=1) means psi is easier to identify"
	) +
	theme_minimal(base_size = 13) +
	theme(legend.position = "bottom")
save_fig(fig3, "psi_identifiability_icc", width = 14, height = 5)

# Print ICC range summary
cat("\n--- ICC signal range ---\n")
cat(sprintf("  %-12s  %-10s  ICC(0)  ICC(1)  Delta\n", "pathogen", "obs_delay"))
for (i in seq_len(nrow(icc_range))) {
	r <- icc_range[i, ]
	cat(sprintf("  %-12s  %-10s  %.3f   %.3f   %.3f\n",
	    r$pathogen, r$obs_delay, r$icc_at_0, r$icc_at_1, r$delta_icc))
}

cat("\n=== Power analysis for psi exclusion complete ===\n")
