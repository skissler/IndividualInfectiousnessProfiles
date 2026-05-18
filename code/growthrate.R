library(tidyverse)
source("code/utils.R")
source("code/parameters.R")

# ==============================================================================
# Infinite-population growth rate analysis
#
# Runs infinite-population stochastic simulations (no susceptible depletion)
# to estimate exponential growth rates under different punctuation (psi) values,
# and compares empirical growth rates to the Lotka-Euler theoretical prediction.
# ==============================================================================

# ==============================================================================
# 1. Global parameters
# ==============================================================================

psivals <- c(0, 0.5, 1)

# Growth rate estimation parameters
min_growth_threshold <- 100
growth_window_days   <- 7

# Simulation parameters
nsim_infpop          <- 5000
max_plot_sims_infpop <- 250

# ==============================================================================
# 2. Loop over pathogens
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
Tgen     <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

cat(sprintf("\n===== %s: growth rate analysis =====\n", pathogen))

# Theoretical growth rate: solve R0 * (beta/(beta+r))^alpha = 1  (Lotka-Euler)
r_malthusian <- uniroot(
	function(r) R0 * (beta / (beta + r))^alpha - 1,
	interval = c(0, 10))$root

cat(sprintf("  %s: theoretical growth rate = %.4f\n", pathogen, r_malthusian))

# Compute max_cases from growth window: enough cases to cover establishment
# (min_growth_threshold) plus growth_window_days of exponential growth,
# with a 3x safety margin for stochastic variation.
max_cases <- ceiling(min_growth_threshold * exp(r_malthusian * (growth_window_days + 1)) * 3)

# ==============================================================================
# 2a. Infinite-population stochastic simulations (cached)
# ==============================================================================

cache_infpop <- load_cache_infpop(pathogen, nsim_infpop, max_cases, psivals)

if (!is.null(cache_infpop)) {
	infpop_summary_df <- cache_infpop$summary
	infpop_plot_df    <- cache_infpop$plot
} else {
	summary_list_ip <- vector("list", nsim_infpop * length(psivals))
	plot_list_ip    <- vector("list", max_plot_sims_infpop * length(psivals))
	summindex_ip <- 0L
	plotindex_ip <- 0L

	for (sim in 1:nsim_infpop) {
	    for (psi in psivals) {
	        infection_times <- sim_infinite_pop(
	            max_cases = max_cases,
	            gen_inf_attempts = gen_inf_attempts_gamma(Tgen, R0, alpha, psi))

	        n_infected <- length(infection_times)

	        # Compute growth rate
	        growthrate <- compute_growth_rate(infection_times, min_growth_threshold, growth_window_days)

	        # Store summary row
	        summindex_ip <- summindex_ip + 1L
	        summary_list_ip[[summindex_ip]] <- tibble(
	            sim = sim, psi = psi,
	            n_infected = n_infected,
	            growthrate = growthrate
	        )

	        # Keep full trajectory for plotting subset
	        if (sim <= max_plot_sims_infpop && n_infected > 0) {
	            plotindex_ip <- plotindex_ip + 1L
	            plot_list_ip[[plotindex_ip]] <- tibble(
	                tinf = infection_times,
	                cuminf = seq_along(infection_times),
	                sim = sim,
	                psi = psi
	            )
	        }
	    }
	    if (sim %% 100 == 0) cat(sprintf("  %s [infpop]: sim %d/%d\n", pathogen, sim, nsim_infpop))
	}

	infpop_summary_df <- bind_rows(summary_list_ip[1:summindex_ip])
	infpop_plot_df    <- bind_rows(plot_list_ip[seq_len(plotindex_ip)])

	# Save cache
	write_csv(infpop_summary_df,
	          cache_path_infpop_summary(pathogen, nsim_infpop, max_cases))
	write_csv(infpop_plot_df,
	          cache_path_infpop_plot(pathogen, nsim_infpop, max_cases))
	cat(sprintf("  %s: infpop simulations saved\n", pathogen))

	infpop_summary_df <- infpop_summary_df %>% mutate(psi = factor(psi))
	infpop_plot_df    <- infpop_plot_df %>% mutate(psi = factor(psi))
}

# ==============================================================================
# 2b. Growth rate analysis and figures
# ==============================================================================

infpop_growthrate_df <- infpop_summary_df %>%
	filter(!is.na(growthrate))

infpop_growthrate_table <- infpop_growthrate_df %>%
	group_by(psi) %>%
	summarise(mean = mean(growthrate), sd = sd(growthrate), .groups = "drop")

cat(sprintf("  %s [infpop]: theoretical growth rate = %.4f\n", pathogen, r_malthusian))
print(infpop_growthrate_table)

# Histogram of growth rates by psi
fig_growthrate_infpop_hists <- ggplot(infpop_growthrate_df, aes(x = growthrate)) +
	geom_histogram(aes(y = after_stat(density)), bins = 40,
	               fill = "white", col = "darkgrey") +
	geom_density(adjust = 2) +
	geom_vline(xintercept = r_malthusian, col = "blue", lty = "dashed", linewidth = 0.8) +
	# geom_vline(data = infpop_growthrate_table, aes(xintercept = mean),
	           # col = "red", linewidth = 0.8) +
	theme_classic() +
	facet_wrap(~psi, nrow = 1) +
	labs(x = "Empirical growth rate (1/day)", y = "Density",
	     title = paste0(pathogen, " (infinite pop)"))

save_fig(fig_growthrate_infpop_hists, paste0("fig_growthrate_infpop_hists_", pathogen))

# Overlaid histogram of growth rates by psi
fig_growthrate_infpop_hists_overlay <- ggplot(infpop_growthrate_df,
                                             aes(x = growthrate, fill = psi, col = psi)) +
	geom_histogram(aes(y = after_stat(density)), bins = 40,
	               alpha = 0.3, position = "identity") +
	geom_density(adjust = 2, linewidth = 0.8, fill = NA) +
	geom_vline(xintercept = r_malthusian, col = "black", lty = "dashed", linewidth = 0.8) +
	theme_classic() +
	labs(x = "Empirical growth rate (1/day)", y = "Density",
	     fill = expression(psi), col = expression(psi),
	     title = paste0(pathogen, " (infinite pop)"))

save_fig(fig_growthrate_infpop_hists_overlay, paste0("fig_growthrate_infpop_hists_overlay_", pathogen))

# Growth rate lines: daily incidence on log scale in the growth window
# Find the first full calendar day after reaching min_growth_threshold
start_days <- infpop_plot_df %>%
	filter(cuminf == min_growth_threshold) %>%
	transmute(sim, psi, start_day = floor(tinf) + 1L)

# Count ALL infections per day within the fixed-length window
daily_counts <- infpop_plot_df %>%
	mutate(day = floor(tinf)) %>%
	inner_join(start_days, by = c("sim", "psi")) %>%
	filter(day >= start_day, day < start_day + growth_window_days) %>%
	group_by(sim, psi, day) %>%
	summarise(count = n(), .groups = "drop")

# Build complete grid and fill zero-count days
infpop_growth_incidence <- start_days %>%
	group_by(sim, psi) %>%
	reframe(day = seq(start_day, start_day + growth_window_days - 1)) %>%
	left_join(daily_counts, by = c("sim", "psi", "day")) %>%
	replace_na(list(count = 0)) %>%
	group_by(sim, psi) %>%
	mutate(day0 = day - min(day)) %>%
	ungroup()

# Fit Poisson GLM per psi to get empirical intercept for reference line
refline_df <- infpop_growth_incidence %>%
	group_by(psi) %>%
	summarise(
		intercept = coef(glm(count ~ day0, family = poisson))[1] / log(10),
		slope = r_malthusian / log(10),
		.groups = "drop")

fig_growthrate_infpop_lines <- infpop_growth_incidence %>%
	filter(count > 0) %>%
	ggplot(aes(x = day0, y = count, group = factor(sim))) +
		geom_line(alpha = 0.1, linewidth = 0.3, col = "grey") +
		geom_point(alpha = 0.2, size = 0.3, col = "grey") +
		geom_abline(data = refline_df, aes(intercept = intercept, slope = slope),
		            col = "blue", linewidth = 0.8, lty = "dashed", inherit.aes = FALSE) +
		scale_y_log10() +
		theme_classic() +
		facet_wrap(~psi, nrow = 1) +
		labs(x = sprintf("Days since case %d", min_growth_threshold),
		     y = "Daily incidence",
		     title = paste0(pathogen, " (infinite pop)"))

save_fig(fig_growthrate_infpop_lines, paste0("fig_growthrate_infpop_lines_", pathogen))

cat(sprintf("  %s: figures saved.\n", pathogen))

} # end pathogen loop
