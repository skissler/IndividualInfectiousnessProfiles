library(tidyverse)
library(odin)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# ==============================================================================
# Stochastic epidemic simulations: uncontrolled dynamics
#
# Corresponds to writeup sections:
#   - "The impact of individual infectiousness profile on uncontrolled
#      epidemic dynamics"
#
# Runs nsim stochastic epidemics for spike (psi=0), smooth (psi=1), and
# intermediate (psi=0.5) individual infectiousness profiles, and:
#   - Plots epidemic trajectories (daily incidence and cumulative incidence)
#   - Calculates extinction probabilities and final size distributions
# ==============================================================================

# ==============================================================================
# 1. Set parameters
# ==============================================================================

psivals <- c(0, 0.5, 1)

# ==============================================================================
# 2. Loop over pathogens
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
Tgen     <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

cat(sprintf("\n===== %s: T=%.2f, alpha=%.2f, R0=%.1f =====\n",
            pathogen, Tgen, alpha, R0))

# ==============================================================================
# 2a. Mean-field solution
# ==============================================================================

# Run the renewal equation model
ren_out <- renewal_epidemic(R0, alpha, Tgen, popsize)

# Aggregate renewal equation output to daily new-infection counts
ren_daily <- ren_out %>%
	mutate(day=floor(t)) %>%
	group_by(day) %>%
	summarise(cuminf=max(cuminf)) %>%
	arrange(day) %>%
	mutate(newinf=cuminf-lag(cuminf)) %>%
	mutate(newinf=case_when(is.na(newinf)~cuminf, TRUE~newinf))

# ==============================================================================
# 2b. Stochastic simulations (cached, incremental processing)
# ==============================================================================

cache <- load_cache(pathogen, nsim, popsize, psivals)

if (!is.null(cache)) {
	sim_summary_df <- cache$summary
	plot_df        <- cache$plot
} else {
	# Run simulations with incremental processing:
	#   - Compute per-sim summary immediately, discard raw infection times
	#   - Keep full trajectories only for sim <= max_plot_sims (for plotting)
	summary_list <- vector("list", nsim * length(psivals))
	plot_list    <- vector("list", max_plot_sims * length(psivals))
	summindex <- 0L
	plotindex <- 0L

	for(sim in 1:nsim){
	    for(psi in psivals){
	        tinf <- sim_stochastic_fast(n=popsize,
	                                    gen_inf_attempts=gen_inf_attempts_gamma(Tgen, R0, alpha, psi))

	        # Extract sorted finite infection times
	        infection_times <- sort(tinf[is.finite(tinf)])
	        final_size <- length(infection_times)
	        established <- as.integer(final_size >= establishment_threshold)

	        # Compute time to establishment threshold
	        establishment_time <- if(final_size >= establishment_threshold) infection_times[establishment_threshold] else NA_real_

	        # Store summary row
	        summindex <- summindex + 1L
	        summary_list[[summindex]] <- tibble(
	            sim = sim,
	            psi = psi,
	            established = established,
	            final_size = final_size,
	            establishment_time = establishment_time
	        )

	        # Keep full trajectories for a subset of established epidemics for plotting 
	        if(sim <= max_plot_sims && established == 1){
	            plotindex <- plotindex + 1L
	            plot_list[[plotindex]] <- tibble(
	                sim = sim,
	                psi = psi,
	                tinf = infection_times,
	                cuminf = seq_along(infection_times),
	            )
	        }
	    }
	    if(sim %% 100 == 0) cat(sprintf("  %s: sim %d/%d\n", pathogen, sim, nsim))
	}

	sim_summary_df <- bind_rows(summary_list[1:summindex])
	plot_df        <- bind_rows(plot_list[seq_len(plotindex)])

	# Save to cache
	write_csv(sim_summary_df, cache_path_summary(pathogen, nsim, popsize))
	write_csv(plot_df, cache_path_plot(pathogen, nsim, popsize))
	cat(sprintf("  %s: simulations saved\n", pathogen))

	sim_summary_df <- sim_summary_df %>% mutate(psi = factor(psi))
	plot_df        <- plot_df %>% mutate(psi = factor(psi))
}

# ==============================================================================
# 2c. Aggregate plot subset to daily resolution (for trajectory plots)
# ==============================================================================

lastday <- ceiling(max(plot_df$tinf))

dayjoin <- expand_grid(
	psi=unique(plot_df$psi),
	sim=unique(plot_df$sim),
	day=0:lastday)

dailyinf_df <- plot_df %>%
	mutate(day = floor(tinf)) %>%
	group_by(psi, sim, day) %>%
	summarise(ninf = n(), .groups = "drop") %>%
	right_join(dayjoin, by = c("psi", "sim", "day")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(psi, sim) %>%
	arrange(day, .by_group = TRUE) %>%
	mutate(cuminf = cumsum(ninf))

# ==============================================================================
# 2d. Key metrics (from summary)
# ==============================================================================

pest_table <- sim_summary_df %>%
	group_by(psi) %>%
	summarise(pestablished=mean(established), .groups="drop")

fs_table <- sim_summary_df %>%
	filter(established==1) %>%
	group_by(psi) %>%
	summarise(fs_mean=mean(final_size), fs_sd=sd(final_size), .groups="drop")

est_time_table <- sim_summary_df %>%
	filter(established==1, !is.na(establishment_time)) %>%
	group_by(psi) %>%
	summarise(
		et_mean = mean(establishment_time),
		et_sd   = sd(establishment_time),
		et_q05  = quantile(establishment_time, 0.05),
		et_q95  = quantile(establishment_time, 0.95),
		.groups = "drop")

cat(sprintf("  %s: P(established) by psi:\n", pathogen))
print(pest_table)
cat(sprintf("  %s: Final size summary:\n", pathogen))
print(fs_table)
cat(sprintf("  %s: Time to %d cases (established epidemics):\n", pathogen, establishment_threshold))
print(est_time_table)

# ==============================================================================
# 2e. Figures — epidemic trajectories (plot subset only)
# ==============================================================================

# Cumulative stochastic epidemic curves (grey) with ODE overlay (black)
fig_cuminf_overlay <- plot_df %>%
	ggplot(aes(x=tinf, y=cuminf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_out, t<=lastday),
			aes(x=t, y=cuminf*popsize, group=1),
			alpha=0.8, linewidth=1, col="black") +
		theme_classic() +
		labs(x="Time (days)", y="Cumulative infections", title = pathogen) +
		facet_wrap(~psi, nrow=1)

save_fig(fig_cuminf_overlay, paste0("fig_cuminf_overlay_", pathogen))

# Cumulative curves with time-to-threshold milestone annotations
milestone_df <- plot_df %>%
	filter(cuminf >= establishment_threshold) %>%
	group_by(sim, psi) %>%
	slice(1) %>%
	group_by(psi) %>%
	summarise(
		mean_t = mean(tinf),
		q05_t  = quantile(tinf, 0.05),
		q95_t  = quantile(tinf, 0.95),
		.groups  = "drop")

# cat(sprintf("  %s: Time to reach 5% of the population infected:\n", pathogen))
# print(milestone_df, n=Inf)

fig_cuminf_milestones <- plot_df %>%
	ggplot(aes(x=tinf, y=cuminf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_out, t<=lastday),
			aes(x=t, y=cuminf*popsize, group=1),
			alpha=0.8, linewidth=1, col="black") +
		geom_segment(data=milestone_df,
			aes(x=q05_t, xend=q05_t,
			    y=establishment_threshold*0.85, yend=establishment_threshold*1.15, group=1),
			col="red", linewidth=0.6) +
		geom_segment(data=milestone_df,
			aes(x=q95_t, xend=q95_t,
			    y=establishment_threshold*0.85, yend=establishment_threshold*1.15, group=1),
			col="red", linewidth=0.6) +
		geom_point(data=milestone_df,
			aes(x=mean_t, y=establishment_threshold, group=1),
			col="red", size=2.5) +
		theme_classic() +
		labs(x="Time (days)", y="Cumulative infections", title = pathogen) +
		facet_wrap(~psi, nrow=1)

save_fig(fig_cuminf_milestones, paste0("fig_cuminf_milestones_", pathogen))

# Daily stochastic incidence curves (grey) with ODE overlay (black)
fig_inf_overlay <- dailyinf_df %>%
	ggplot(aes(x=day, y=ninf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_daily, day<=lastday),
			aes(x=day, y=newinf*popsize, group=1),
			alpha=0.8, linewidth=1, col="black") +
		theme_classic() +
		labs(x="Time (days)", y="Daily new infections", title = pathogen) +
		facet_wrap(~psi, nrow=1)

save_fig(fig_inf_overlay, paste0("fig_inf_overlay_", pathogen))

cat(sprintf("  %s: figures saved.\n", pathogen))

} # end pathogen loop
