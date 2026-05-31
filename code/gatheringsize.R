# ==============================================================================
# Gathering size restrictions: three-scenario comparison
# ==============================================================================
#
# Compares establishment probability under:
#   1. "Restricted" — truncated contact process (true R0_new, true k_new)
#   2. "R0-only"   — original contact process, R0 scaled down to R0_new
#   3. "Poisson"   — no contact process, Poisson(R0_new) offspring
#
# Expected ordering: P_est_R0only < P_est_restricted < P_est_Poisson,
# with gaps larger for small psi (punctuated profiles).
# ==============================================================================

library(tidyverse)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

set.seed(42)

# --- Parameters ---
k_c_gs     <- 1                          # moderate contact heterogeneity
lambda_gs  <- 1                          # Poisson switching rate
c_max_vals <- c(1.5, 2, 3, 5, Inf)      # gathering size caps
psi_vals   <- c(0, 0.25, 0.5, 0.75, 1)
c_max_main <- 2                          # focal restriction for main figure

# ==============================================================================
# Start pathogen loop
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
Tgen     <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

cat(sprintf("\n===== Gathering size: %s (T=%.2f, alpha=%.2f, R0=%.1f) =====\n",
            pathogen, Tgen, alpha, R0))

# --- Summary table: analytical mu_T and R0_new for each c_max ---
cat("\n  Summary table: mu_T, R0_new by c_max\n")
cat(sprintf("  %-8s  %-8s  %-8s\n", "c_max", "mu_T", "R0_new"))

for (c_max_val in c_max_vals) {
	mu_T   <- mu_truncated_gamma(k_c_gs, c_max_val)
	R0_new <- R0 * mu_T
	cat(sprintf("  %-8s  %-8.4f  %-8.3f\n",
	    ifelse(is.infinite(c_max_val), "Inf", sprintf("%.1f", c_max_val)),
	    mu_T, R0_new))
}

# ==============================================================================
# Three-scenario epidemic simulations
# ==============================================================================

total_sims <- length(c_max_vals) * length(psi_vals) * 3 * nsim_small
epi_results <- vector("list", total_sims)
epi_idx     <- 0L

for (c_max_val in c_max_vals) {

	mu_T   <- mu_truncated_gamma(k_c_gs, c_max_val)
	R0_new <- R0 * mu_T

	cat(sprintf("\n  c_max = %s (mu_T = %.4f, R0_new = %.3f)\n",
	    ifelse(is.infinite(c_max_val), "Inf", sprintf("%.1f", c_max_val)),
	    mu_T, R0_new))

	for (psi_val in psi_vals) {

		# --- Scenario 1: Restricted (truncated contact process) ---
		gfun_restricted <- gen_inf_attempts_gammapoisson_contacts(
			Tgen, R0, alpha, psi_val, k_c_gs, lambda_gs, c_max = c_max_val)

		# --- Scenario 2: R0-only (original contacts, scaled R0) ---
		gfun_r0only <- gen_inf_attempts_gammapoisson_contacts(
			Tgen, R0_new, alpha, psi_val, k_c_gs, lambda_gs)

		# --- Scenario 3: Poisson (no contact process) ---
		gfun_poisson <- gen_inf_attempts_gamma(Tgen, R0_new, alpha, psi_val)

		scenarios <- list(
			list(name = "Restricted", gfun = gfun_restricted),
			list(name = "R0-only",    gfun = gfun_r0only),
			list(name = "Poisson",    gfun = gfun_poisson)
		)

		for (sc in scenarios) {
			for (sim in seq_len(nsim_small)) {
				tinf <- sim_stochastic_fast(n = popsize,
				                            gen_inf_attempts = sc$gfun,
				                            maxinf = establishment_threshold)
				n_infected <- sum(is.finite(tinf))

				epi_idx <- epi_idx + 1L
				epi_results[[epi_idx]] <- tibble(
					c_max      = c_max_val,
					psi        = psi_val,
					scenario   = sc$name,
					sim        = sim,
					n_infected = n_infected,
					established = as.integer(n_infected >= establishment_threshold))
			}
		}

		cat(sprintf("    psi = %.2f done\n", psi_val))
	}
}

epi_df <- bind_rows(epi_results)

# ==============================================================================
# Compute P(establishment) summary
# ==============================================================================

pest_summary <- epi_df %>%
	group_by(c_max, psi, scenario) %>%
	summarise(
		p_est  = mean(established),
		n_sims = n(),
		se     = sqrt(p_est * (1 - p_est) / n_sims),
		.groups = "drop"
	) %>%
	mutate(scenario = factor(scenario, levels = c("Poisson", "Restricted", "R0-only")))

cat(sprintf("\n  %s: P(establishment) summary\n", pathogen))
pest_summary %>%
	select(c_max, psi, scenario, p_est, se) %>%
	print(n = Inf)

# ==============================================================================
# Verification: simulated mean offspring vs analytical R0_new
# ==============================================================================

cat(sprintf("\n  %s: R0 cross-check (simulated mean offspring vs analytical R0_new)\n", pathogen))
for (c_max_val in c_max_vals) {
	mu_T   <- mu_truncated_gamma(k_c_gs, c_max_val)
	R0_new <- R0 * mu_T

	gfun_check <- gen_inf_attempts_gammapoisson_contacts(
		Tgen, R0, alpha, 1, k_c_gs, lambda_gs, c_max = c_max_val)
	n_offspring <- sapply(rep(0, 5000), function(t) length(gfun_check(t)))
	R0_sim <- mean(n_offspring)

	cat(sprintf("    c_max = %-5s: R0_new(theory) = %.3f, R0_new(sim) = %.3f\n",
	    ifelse(is.infinite(c_max_val), "Inf", sprintf("%.1f", c_max_val)),
	    R0_new, R0_sim))
}

# ==============================================================================
# Main figure: three scenarios at fixed c_max, across psi
# ==============================================================================

pest_main <- pest_summary %>%
	filter(c_max == c_max_main)

mu_T_main   <- mu_truncated_gamma(k_c_gs, c_max_main)
R0_new_main <- R0 * mu_T_main

fig_gathering_pest <- ggplot(pest_main,
		aes(x = factor(psi), y = p_est, fill = scenario)) +
	geom_col(position = position_dodge(width = 0.7), width = 0.6) +
	geom_errorbar(aes(ymin = pmax(p_est - 1.96 * se, 0),
	                  ymax = pmin(p_est + 1.96 * se, 1)),
	              position = position_dodge(width = 0.7), width = 0.2) +
	scale_fill_manual(
		values = c("Poisson" = "#2196F3", "Restricted" = "#FF9800", "R0-only" = "#4CAF50"),
		name = "Scenario") +
	theme_classic() +
	labs(x = expression(psi),
	     y = expression(P(establishment)),
	     title = sprintf("%s: gathering cap c_max = %g (R0: %.1f -> %.2f)",
	                     pathogen, c_max_main, R0, R0_new_main))

save_fig(fig_gathering_pest, sprintf("fig_gathering_pest_%s", pathogen), width = 8, height = 5)
cat(sprintf("  Saved fig_gathering_pest_%s\n", pathogen))

# ==============================================================================
# Secondary figure: P(establishment) vs c_max for each psi (restricted only)
# ==============================================================================

pest_cmax <- pest_summary %>%
	filter(scenario == "Restricted", is.finite(c_max))

pest_baseline <- pest_summary %>%
	filter(scenario == "Restricted", is.infinite(c_max)) %>%
	select(psi, p_est_baseline = p_est)

fig_gathering_cmax <- ggplot(pest_cmax,
		aes(x = c_max, y = p_est, color = factor(psi), group = factor(psi))) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	geom_errorbar(aes(ymin = pmax(p_est - 1.96 * se, 0),
	                  ymax = pmin(p_est + 1.96 * se, 1)),
	              width = 0.1) +
	geom_hline(data = pest_baseline,
	           aes(yintercept = p_est_baseline, color = factor(psi)),
	           linetype = "dashed", alpha = 0.5) +
	scale_color_viridis_d(option = "plasma", name = expression(psi)) +
	theme_classic() +
	labs(x = expression(c[max]),
	     y = expression(P(establishment)),
	     title = sprintf("%s: restricted model by gathering cap", pathogen))

save_fig(fig_gathering_cmax, sprintf("fig_gathering_cmax_%s", pathogen), width = 8, height = 5)
cat(sprintf("  Saved fig_gathering_cmax_%s\n", pathogen))

} # end pathogen loop