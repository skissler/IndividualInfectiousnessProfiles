# ============================================================
# Monte Carlo MI comparison: clustered vs iid serial intervals
# ============================================================
#
# For Omicron parameters, compare the mutual information between
# the GI shape parameter (alpha) and observed serial intervals
# under two designs:
#   1. n iid serial intervals (independent transmission pairs)
#   2. K=50 clusters of siblings (with psi controlling correlation)
#
# MI(alpha; X) = E[log p(X|alpha) - log p(X)]
# estimated purely by Monte Carlo sampling.
#
# Beta is fixed at its true value to keep this a 1D problem.
# ============================================================

library(tidyverse)
source("code/utils.R")
source("code/parameters.R")

# --- Parameters ---
pars <- parslist[[2]]  # Omicron
alpha_true <- pars$alpha
beta_true  <- pars$beta
R0 <- pars$R0
a_obs <- 4; b_obs <- 1
K <- 50

# Prior on alpha (beta fixed)
alpha_lo <- 0.5
alpha_hi <- 5.0

# MC tuning
n_outer <- 100    # joint (alpha, X) draws for MI
n_inner <- 1000   # prior draws for marginal likelihood
R_mc    <- 50     # MC draws for cluster likelihood integral

cat(sprintf("Omicron: alpha=%.3f, beta=%.4f, R0=%d\n", alpha_true, beta_true, R0))
cat(sprintf("Prior: alpha ~ U(%.1f, %.1f), beta fixed\n", alpha_lo, alpha_hi))
cat(sprintf("MC: n_outer=%d, n_inner=%d, R_mc=%d\n\n", n_outer, n_inner, R_mc))

# ============================================================
# 1. Density helpers (one-time precomputation)
# ============================================================

# Density of X-Y where X~Gamma(s1,r1), Y~Gamma(s2,r2)
dgamma_diff_local <- function(d, s1, r1, s2, r2) {
	sapply(d, function(di) {
		f <- function(y) dgamma(di + y, s1, r1) * dgamma(y, s2, r2)
		tryCatch(integrate(f, max(0, -di), Inf, rel.tol = 1e-8, abs.tol = 1e-12)$value,
		         error = function(e) 0)
	})
}

# Pre-build f_S: density of s = Gamma(alpha, beta) + Gamma(a,b) - Gamma(a,b)
# Returns a lookup function: f_S(si_vector, alpha) -> density values
build_f_S_table <- function(a_obs, b_obs, alpha_grid, beta) {
	cat("  Building f_S table...")
	sd_d <- sqrt(a_obs) / b_obs
	max_mean <- max(alpha_grid) / beta
	max_sd <- sqrt(max(alpha_grid)) / beta

	tau_grid <- seq(1e-4, max_mean + 8 * max_sd, length.out = 401)
	dtau <- diff(tau_grid)[1]
	s_lo <- -8 * sd_d
	s_hi <- max_mean + 8 * (max_sd + sd_d)
	s_grid <- seq(s_lo, s_hi, length.out = 501)

	ddiff_ext <- max(abs(s_lo), s_hi) + max(tau_grid)
	ddiff_grid <- seq(-ddiff_ext, ddiff_ext, length.out = 2001)
	ddiff_vals <- dgamma_diff_local(ddiff_grid, a_obs, b_obs, a_obs, b_obs)
	ddiff_fn <- approxfun(ddiff_grid, ddiff_vals, rule = 2, yleft = 0, yright = 0)

	D <- outer(s_grid, tau_grid, function(s, t) ddiff_fn(s - t))

	fns <- vector("list", length(alpha_grid))
	for (i in seq_along(alpha_grid)) {
		g <- dgamma(tau_grid, shape = alpha_grid[i], rate = beta)
		g[!is.finite(g)] <- 0
		fv <- pmax(as.numeric(D %*% g) * dtau, 0)
		fns[[i]] <- approxfun(s_grid, fv, rule = 2, yleft = 0, yright = 0)
	}
	cat(" done.\n")

	function(si, alpha) {
		k <- which.min(abs(alpha_grid - alpha))
		out <- fns[[k]](si)
		out[out < .Machine$double.xmin] <- .Machine$double.xmin
		out
	}
}

# Pre-build f_+: density of Gamma(psi*alpha, beta) + Gamma(a_obs, b_obs)
# Returns a lookup function: get_f_plus(alpha) -> approxfun
build_f_plus_table <- function(a_obs, b_obs, psi, alpha_grid, beta) {
	cat(sprintf("  Building f_+ table (psi=%.1f)...", psi))
	sd_d <- sqrt(a_obs) / b_obs
	shapes <- pmax(psi * alpha_grid, 0.01)
	max_mean <- max(shapes) / beta + a_obs / b_obs
	max_sd <- sqrt(max(shapes)) / beta + sd_d

	t_grid <- seq(1e-4, max(shapes) / beta + 6 * sqrt(max(shapes)) / beta,
	              length.out = 301)
	dt <- diff(t_grid)[1]
	x_grid <- seq(0, max_mean + 8 * max_sd, length.out = 401)

	D <- outer(x_grid, t_grid, function(x, t) dgamma(x - t, shape = a_obs, rate = b_obs))

	fns <- vector("list", length(alpha_grid))
	for (i in seq_along(alpha_grid)) {
		g <- dgamma(t_grid, shape = shapes[i], rate = beta)
		g[!is.finite(g)] <- 0
		fv <- pmax(as.numeric(D %*% g) * dt, 0)
		fns[[i]] <- approxfun(x_grid, fv, rule = 2, yleft = 0, yright = 0)
	}
	cat(" done.\n")

	function(alpha) {
		k <- which.min(abs(alpha_grid - alpha))
		fns[[k]]
	}
}

# ============================================================
# 2. Log-likelihood functions
# ============================================================

loglik_iid <- function(alpha, si, f_S) {
	sum(log(f_S(si, alpha)))
}

loglik_cluster <- function(alpha, clusters, psi, get_f_plus) {
	f_plus <- get_f_plus(alpha)
	shape_l <- max((1 - psi) * alpha, 0.001)

	ll <- 0
	for (k in seq_along(clusters)) {
		si <- clusters[[k]]
		m <- length(si)
		if (m == 0) next

		l  <- rgamma(R_mc, shape_l, beta_true)
		d0 <- rgamma(R_mc, a_obs, b_obs)
		shift <- -l + d0

		shifted <- outer(si, shift, "+")  # m x R_mc
		fv <- f_plus(as.vector(shifted))
		fv[fv < .Machine$double.xmin] <- .Machine$double.xmin
		log_f <- matrix(log(fv), m, R_mc)
		log_prod <- colSums(log_f)

		mx <- max(log_prod)
		ll <- ll + mx + log(mean(exp(log_prod - mx)))
	}
	ll
}

# ============================================================
# 3. MI estimation
# ============================================================

compute_mi <- function(simulate_fn, loglik_fn, label = "") {
	contribs <- numeric(n_outer)

	for (i in seq_len(n_outer)) {
		alpha_i <- runif(1, alpha_lo, alpha_hi)
		data_i  <- simulate_fn(alpha_i)

		ll_true <- loglik_fn(alpha_i, data_i)

		alpha_j <- runif(n_inner, alpha_lo, alpha_hi)
		ll_j <- sapply(alpha_j, function(a) loglik_fn(a, data_i))

		ok <- is.finite(ll_j)
		if (sum(ok) < 10) { contribs[i] <- NA; next }
		max_ll <- max(ll_j[ok])
		log_marg <- max_ll + log(mean(exp(ll_j[ok] - max_ll)))

		contribs[i] <- ll_true - log_marg
		if (i %% 20 == 0) cat(".")
	}
	cat("\n")

	list(mi = mean(contribs, na.rm = TRUE),
	     se = sd(contribs, na.rm = TRUE) / sqrt(sum(!is.na(contribs))),
	     raw = contribs)
}

# ============================================================
# 4. Main comparison
# ============================================================

alpha_table <- seq(alpha_lo, alpha_hi, length.out = 200)
f_S <- build_f_S_table(a_obs, b_obs, alpha_table, beta_true)

# --- IID MI at various n (only need to compute once) ---
n_vals <- c(20, 50, 100, 200, 300)
mi_iid <- list()

for (n in n_vals) {
	cat(sprintf("IID n=%d: ", n))
	sim_fn <- function(alpha) {
		tau <- rgamma(n, alpha, beta_true)
		d1  <- rgamma(n, a_obs, b_obs)
		d0  <- rgamma(n, a_obs, b_obs)
		tau + d1 - d0
	}
	ll_fn <- function(alpha, data) loglik_iid(alpha, data, f_S)
	mi_iid[[as.character(n)]] <- compute_mi(sim_fn, ll_fn)
}

# --- Clustered MI at various psi ---
psi_vals <- c(0.3, 0.5, 0.7, 1.0)
mi_cluster <- list()

for (psi in psi_vals) {
	cat(sprintf("\nClustered psi=%.1f (K=%d): ", psi, K))

	get_f_plus <- build_f_plus_table(a_obs, b_obs, psi, alpha_table, beta_true)

	sim_fn <- function(alpha) {
		simulate_clusters(K, psi, R0, alpha, beta_true, a_obs, b_obs, p_asc = 1)
	}
	ll_fn <- function(alpha, data) loglik_cluster(alpha, data, psi, get_f_plus)
	mi_cluster[[as.character(psi)]] <- compute_mi(sim_fn, ll_fn)
}

# ============================================================
# 5. Results
# ============================================================

cat("\n========================================\n")
cat("MI comparison: clustered vs iid\n")
cat("========================================\n\n")

cat("IID serial intervals:\n")
for (n in n_vals) {
	r <- mi_iid[[as.character(n)]]
	cat(sprintf("  n = %3d:  MI = %.3f +/- %.3f nats\n", n, r$mi, r$se))
}

cat(sprintf("\nClustered (K=%d index cases, expected n ~ %d):\n", K, K * R0))
for (psi in psi_vals) {
	r <- mi_cluster[[as.character(psi)]]
	cat(sprintf("  psi = %.1f: MI = %.3f +/- %.3f nats\n", psi, r$mi, r$se))
}

# Interpolate n_eff
cat("\nEffective sample sizes:\n")
mi_iid_vals <- sapply(mi_iid, function(r) r$mi)
mi_iid_ns   <- n_vals

for (psi in psi_vals) {
	mi_cl <- mi_cluster[[as.character(psi)]]$mi
	if (mi_cl >= min(mi_iid_vals) && mi_cl <= max(mi_iid_vals)) {
		neff_fn <- approxfun(mi_iid_vals, mi_iid_ns, rule = 2)
		neff <- neff_fn(mi_cl)
		cat(sprintf("  psi = %.1f: n_eff = %.0f  (out of ~%d total)\n",
		    psi, neff, K * R0))
	} else if (mi_cl > max(mi_iid_vals)) {
		cat(sprintf("  psi = %.1f: n_eff > %d  (cluster MI exceeds iid range)\n",
		    psi, max(n_vals)))
	} else {
		cat(sprintf("  psi = %.1f: n_eff < %d\n", psi, min(n_vals)))
	}
}

cat("\nDone.\n")
