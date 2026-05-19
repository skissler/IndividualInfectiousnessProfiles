# ==============================================================================
# Empirical estimation of psi from OutbreakTrees transmission tree data
# ==============================================================================
#
# Uses the cluster-based likelihood from psi_identifiability.R to estimate psi
# for real pathogens using published transmission trees with symptom onset data.
#
# Data source: OutbreakTrees database (Taube et al., PLOS Biology 2022)
#   https://outbreaktrees.ecology.uga.edu
#   https://github.com/DrakeLab/taube-transmission-trees
#
# Depends on: save_fig (from utils.R)
# ==============================================================================

cat("=== Empirical psi estimation from OutbreakTrees ===\n")

library(igraph)

# ==============================================================================
# 1. Incubation period parameters (Gamma shape & rate) from published estimates
# ==============================================================================
#
# We parameterize incubation periods as Gamma(a_obs, b_obs) with
# mean = a_obs/b_obs, var = a_obs/b_obs^2.
#
# Given published (mean, sd), we moment-match:
#   a_obs = mean^2 / var,  b_obs = mean / var
#
# Sources:
#   COVID-19: Li et al. (2020) NEJM, mean 5.2d, sd 2.8d (early Wuhan)
#   MERS: Assiri et al. (2013) Lancet ID, mean 5.2d, sd 2.5d
#   Measles: Lessler et al. (2009) Lancet ID, mean 11.5d, sd 2.2d
#   Ebola: WHO Ebola Response Team (2014) NEJM, mean 11.4d, sd 6.5d
#   Pneumonic plague: Gani & Leach (2004) EID, mean 3.5d, sd 1.5d
#   Norovirus: Lee et al. (2013) BMC ID, mean 1.2d, sd 0.5d
#   Nipah: Nikolay et al. (2019) NEJM, mean 9.0d, sd 4.0d
#   Smallpox: Nishiura & Eichner (2007) Int J Hyg, mean 12.5d, sd 2.2d
#   Hepatitis A: CDC, mean 28d, sd 7d (range 15-50)
#   H1N1 influenza: Lessler et al. (2009), mean 1.4d, sd 0.5d
#   H1N1/H3N2 influenza: same as H1N1

incubation_params <- list(
	"COVID-19"         = list(mean = 5.2,  sd = 2.8),
	"MERS"             = list(mean = 5.2,  sd = 2.5),
	"Measles"          = list(mean = 11.5, sd = 2.2),
	"Ebola"            = list(mean = 11.4, sd = 6.5),
	"Pneumonic plague" = list(mean = 3.5,  sd = 1.5),
	"Norovirus"        = list(mean = 1.2,  sd = 0.5),
	"Nipah virus"      = list(mean = 9.0,  sd = 4.0),
	"Smallpox"         = list(mean = 12.5, sd = 2.2),
	"Hepatitis A"      = list(mean = 28.0, sd = 7.0),
	"H1N1"             = list(mean = 1.4,  sd = 0.5),
	"H1N1, H3N2"       = list(mean = 1.4,  sd = 0.5)
)

# Moment-match to Gamma(a_obs, b_obs)
for (d in names(incubation_params)) {
	ip <- incubation_params[[d]]
	incubation_params[[d]]$a_obs <- ip$mean^2 / ip$sd^2
	incubation_params[[d]]$b_obs <- ip$mean / ip$sd^2
}

# ==============================================================================
# 2. Extract serial interval clusters from OutbreakTrees
# ==============================================================================

cat("--- Loading and parsing OutbreakTrees data ---\n")

dat <- readRDS("data/outbreaktrees_data.RDS")
has_onset <- grepl("symptom_onset", dat$Attributes)
onset_dat <- dat[has_onset, ]

cat(sprintf("  %d trees with symptom onset data (out of %d total)\n",
    sum(has_onset), nrow(dat)))

#' Parse onset values to numeric days (relative to earliest onset in tree)
parse_onset <- function(onset_raw) {
	n <- length(onset_raw)
	result <- rep(NA_real_, n)

	# Skip non-informative values
	skip <- is.na(onset_raw) |
		grepl("asymptomatic|unclear|before|after|largely|diagnosed",
		      onset_raw, ignore.case = TRUE)
	working <- onset_raw
	working[skip] <- NA

	# Try 1: already numeric (relative days)
	nums <- suppressWarnings(as.numeric(working))
	if (sum(!is.na(nums)) > 0) {
		result[!is.na(nums)] <- nums[!is.na(nums)]
		return(result)
	}

	# Try 2: "Day X" format
	day_pattern <- "^Day ([0-9]+)$"
	day_matches <- grepl(day_pattern, working)
	if (any(day_matches, na.rm = TRUE)) {
		for (i in which(!is.na(working) & day_matches)) {
			result[i] <- as.numeric(sub(day_pattern, "\\1", working[i]))
		}
		if (sum(!is.na(result)) > 0) return(result)
	}

	# Try 3: M/D date format (handle Dec->Jan year boundary)
	dates <- as.Date(working, format = "%m/%d")
	if (sum(!is.na(dates)) > sum(!is.na(working)) * 0.3) {
		month_nums <- as.numeric(format(dates, "%m"))
		has_late <- any(month_nums >= 10, na.rm = TRUE)
		has_early <- any(month_nums <= 3, na.rm = TRUE)
		if (has_late && has_early) {
			early_idx <- which(!is.na(dates) & month_nums <= 3)
			dates[early_idx] <- dates[early_idx] + 365
		}
		ref <- min(dates, na.rm = TRUE)
		result[!is.na(dates)] <- as.numeric(dates[!is.na(dates)] - ref)
		return(result)
	}

	result
}

#' Extract clusters of serial intervals from a transmission tree
#' Returns list of numeric vectors (serial intervals grouped by infector)
extract_clusters <- function(tree) {
	tree <- upgrade_graph(tree)
	onset_raw <- vertex_attr(tree, "onset")
	onset_num <- parse_onset(onset_raw)
	names(onset_num) <- V(tree)$name

	el <- as_edgelist(tree)
	infectors <- unique(el[, 1])

	clusters <- list()
	for (inf in infectors) {
		infectees <- el[el[, 1] == inf, 2]
		onset_inf <- onset_num[inf]
		onset_infectees <- onset_num[infectees]
		if (is.na(onset_inf)) next
		valid <- !is.na(onset_infectees)
		if (sum(valid) < 1) next
		si <- onset_infectees[valid] - onset_inf
		clusters[[length(clusters) + 1]] <- si
	}
	clusters
}

# Process all trees with onset data
all_disease_clusters <- list()
for (i in seq_len(nrow(onset_dat))) {
	disease <- as.character(onset_dat$Disease[i])
	tree_id <- as.character(onset_dat$id[i])

	clusters <- tryCatch(
		extract_clusters(onset_dat$tree[[i]]),
		error = function(e) list()
	)
	if (length(clusters) == 0) next

	# Sanity check: skip trees with clearly broken SI
	si_all <- unlist(clusters)
	if (abs(median(si_all)) > 100) next

	if (!(disease %in% names(all_disease_clusters))) {
		all_disease_clusters[[disease]] <- list(clusters = list(), tree_ids = character(0))
	}
	all_disease_clusters[[disease]]$clusters <- c(
		all_disease_clusters[[disease]]$clusters, clusters
	)
	all_disease_clusters[[disease]]$tree_ids <- c(
		all_disease_clusters[[disease]]$tree_ids, tree_id
	)
}

# ----------------------------------------------------------------------
# 2b. Supplement with additional published datasets
# ----------------------------------------------------------------------

# --- Tianjin COVID-19 clusters (Ganyani et al. 2020 Eurosurveillance) ---
cat("  Supplementing with Tianjin COVID-19 data (Ganyani et al. 2020)...\n")

tj <- tryCatch(read.csv("data/tianjin_covid_data.csv", stringsAsFactors = FALSE),
               error = function(e) NULL)
if (!is.null(tj)) {
	tj_pairs <- list()
	for (i in seq_len(nrow(tj))) {
		if (is.na(tj$Infection_source[i])) next
		infector_ids <- regmatches(tj$Infection_source[i],
		                           gregexpr("TJ[0-9]+", tj$Infection_source[i]))[[1]]
		for (inf_id in infector_ids) {
			tj_pairs[[length(tj_pairs) + 1]] <- data.frame(
				infector = inf_id, infectee = tj$case_id[i],
				stringsAsFactors = FALSE
			)
		}
	}
	if (length(tj_pairs) > 0) {
		tj_pair_df <- do.call(rbind, tj_pairs)
		onset_map <- setNames(tj$symptom_onset, tj$case_id)
		tj_pair_df$onset_inf <- as.Date(onset_map[tj_pair_df$infector], format = "%d/%m/%Y")
		tj_pair_df$onset_ife <- as.Date(onset_map[tj_pair_df$infectee], format = "%d/%m/%Y")
		tj_pair_df$si <- as.numeric(tj_pair_df$onset_ife - tj_pair_df$onset_inf)
		tj_valid <- tj_pair_df[!is.na(tj_pair_df$si), ]

		tj_clusters <- split(tj_valid$si, tj_valid$infector)
		tj_clusters <- unname(tj_clusters)

		if (!("COVID-19" %in% names(all_disease_clusters))) {
			all_disease_clusters[["COVID-19"]] <- list(clusters = list(), tree_ids = character(0))
		}
		all_disease_clusters[["COVID-19"]]$clusters <- c(
			all_disease_clusters[["COVID-19"]]$clusters, tj_clusters
		)
		all_disease_clusters[["COVID-19"]]$tree_ids <- c(
			all_disease_clusters[["COVID-19"]]$tree_ids, "tianjin_ganyani_2020"
		)
		cat(sprintf("    Added %d Tianjin clusters (%d with >=2 offspring), %d SI\n",
		    length(tj_clusters), sum(sapply(tj_clusters, length) >= 2),
		    sum(sapply(tj_clusters, length))))
	}
}

# --- Hart et al. (2022) eLife: transmission pairs ---
cat("  Supplementing with Hart et al. (2022) transmission pairs...\n")

hart <- tryCatch(readxl::read_xlsx("data/hart_transmission_pairs.xlsx"),
                 error = function(e) NULL)
if (!is.null(hart)) {
	# Zhang dataset: 35 pairs from one superspreader
	zhang <- hart[hart$Dataset == "Zhang", ]
	zhang_si <- zhang$t_s2 - zhang$t_s1
	zhang_cluster <- list(zhang_si)

	# Other datasets: group by onset day as putative infector
	other_clusters <- list()
	for (ds in setdiff(unique(hart$Dataset), "Zhang")) {
		sub <- hart[hart$Dataset == ds, ]
		sub$si <- sub$t_s2 - sub$t_s1
		by_infector <- split(sub$si, sub$t_s1)
		other_clusters <- c(other_clusters, unname(by_infector))
	}

	hart_clusters <- c(zhang_cluster, other_clusters)

	if (!("COVID-19" %in% names(all_disease_clusters))) {
		all_disease_clusters[["COVID-19"]] <- list(clusters = list(), tree_ids = character(0))
	}
	all_disease_clusters[["COVID-19"]]$clusters <- c(
		all_disease_clusters[["COVID-19"]]$clusters, hart_clusters
	)
	all_disease_clusters[["COVID-19"]]$tree_ids <- c(
		all_disease_clusters[["COVID-19"]]$tree_ids, "hart_2022_elife"
	)
	cat(sprintf("    Added %d Hart clusters (%d with >=2 offspring), %d SI\n",
	    length(hart_clusters), sum(sapply(hart_clusters, length) >= 2),
	    sum(sapply(hart_clusters, length))))
}

# Summary
cat("\n--- Extracted clusters by disease (with supplements) ---\n")
cat(sprintf("  %-20s  %5s  %5s  %5s  %8s\n",
    "Disease", "Trees", "Clust", "Multi", "TotalSI"))
for (d in names(all_disease_clusters)) {
	dc <- all_disease_clusters[[d]]
	n_trees <- length(unique(dc$tree_ids))
	n_clust <- length(dc$clusters)
	n_multi <- sum(sapply(dc$clusters, length) >= 2)
	n_si <- sum(sapply(dc$clusters, length))
	cat(sprintf("  %-20s  %5d  %5d  %5d  %8d\n",
	    d, n_trees, n_clust, n_multi, n_si))
}

# ==============================================================================
# 3. Density functions (same as psi_identifiability.R)
# ==============================================================================

#' Density of X + Y where X ~ Gamma(shape1, rate1), Y ~ Gamma(shape2, rate2)
dgamma_sum_emp <- function(x, shape1, rate1, shape2, rate2) {
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
dgamma_diff_emp <- function(d, shape1, rate1, shape2, rate2) {
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

make_density_interp_emp <- function(dfun, grid, ...) {
	vals <- dfun(grid, ...)
	approxfun(grid, vals, rule = 2, yleft = 0, yright = 0)
}

#' Log-likelihood for one cluster via delta quadrature
loglik_cluster_emp <- function(s_vec, f_nu_interp, f_delta_interp, delta_grid) {
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

# ==============================================================================
# 4. Precompute density engines and estimate psi for each disease
# ==============================================================================

#' Precompute f_nu and f_delta interpolators for a psi grid
precompute_engines_emp <- function(psi_grid, alpha, beta, a_obs, b_obs) {
	mean_gi <- alpha / beta
	sd_gi   <- sqrt(alpha) / beta
	mean_d  <- a_obs / b_obs
	sd_d    <- sqrt(a_obs) / b_obs

	delta_lo <- -(mean_d + 5 * sd_d)
	delta_hi <- mean_gi + 5 * sd_gi
	delta_grid <- seq(delta_lo, delta_hi, length.out = 401)

	nu_hi <- mean_gi + mean_d + 5 * (sd_gi + sd_d)
	nu_grid <- seq(1e-6, nu_hi, length.out = 501)

	n_cores <- max(1L, parallel::detectCores(logical = FALSE))

	engines <- parallel::mclapply(psi_grid, function(psi) {
		if (psi < 1e-10) {
			f_nu <- approxfun(nu_grid, dgamma(nu_grid, shape = a_obs, rate = b_obs),
			                  rule = 2, yleft = 0, yright = 0)
			f_delta <- make_density_interp_emp(
				dgamma_diff_emp, delta_grid,
				shape1 = alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
		} else if (psi > 1 - 1e-10) {
			f_nu <- make_density_interp_emp(
				dgamma_sum_emp, nu_grid,
				shape1 = alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
			f_delta <- approxfun(delta_grid, dgamma(-delta_grid, shape = a_obs, rate = b_obs),
			                     rule = 2, yleft = 0, yright = 0)
		} else {
			f_nu <- make_density_interp_emp(
				dgamma_sum_emp, nu_grid,
				shape1 = psi * alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
			f_delta <- make_density_interp_emp(
				dgamma_diff_emp, delta_grid,
				shape1 = (1 - psi) * alpha, rate1 = beta,
				shape2 = a_obs, rate2 = b_obs
			)
		}
		list(f_nu = f_nu, f_delta = f_delta)
	}, mc.cores = n_cores)

	list(
		f_nu_interps    = lapply(engines, `[[`, "f_nu"),
		f_delta_interps = lapply(engines, `[[`, "f_delta"),
		delta_grid      = delta_grid
	)
}

#' Compute posterior over psi given observed clusters
compute_psi_posterior_emp <- function(clusters, engines, psi_grid) {
	f_nu_interps    <- engines$f_nu_interps
	f_delta_interps <- engines$f_delta_interps
	delta_grid      <- engines$delta_grid

	logliks <- sapply(seq_along(psi_grid), function(idx) {
		total <- 0
		for (k in seq_along(clusters)) {
			ll <- loglik_cluster_emp(
				clusters[[k]],
				f_nu_interps[[idx]],
				f_delta_interps[[idx]],
				delta_grid
			)
			total <- total + ll
		}
		total
	})

	max_ll <- max(logliks)
	log_post <- logliks - max_ll
	post <- exp(log_post)
	post / sum(post)
}

# ==============================================================================
# 5. Main estimation loop
# ==============================================================================

cat("\n--- Estimating psi for each disease ---\n")

# We need GI parameters (alpha, beta) for each disease.
# Strategy: moment-match from published serial interval mean & SD,
# adjusting for incubation period contribution.
# SI mean ~ GI mean (since E[d_j - d_0] = 0 for symmetric incubation)
# SI var  = GI var + 2 * incubation var
# => GI var = SI var - 2 * incubation var
# => alpha = GI_mean^2 / GI_var,  beta = GI_mean / GI_var

psi_grid <- seq(0, 1, length.out = 101)

# Minimum multi-offspring clusters required
min_multi_clusters <- 5

# Collect results
posterior_results <- list()

for (disease in names(all_disease_clusters)) {
	if (!(disease %in% names(incubation_params))) {
		cat(sprintf("  Skipping %s (no incubation parameters)\n", disease))
		next
	}

	dc <- all_disease_clusters[[disease]]
	clusters <- dc$clusters
	multi_clusters <- clusters[sapply(clusters, length) >= 2]

	if (length(multi_clusters) < min_multi_clusters) {
		cat(sprintf("  Skipping %s (only %d multi-offspring clusters)\n",
		    disease, length(multi_clusters)))
		next
	}

	ip <- incubation_params[[disease]]
	a_obs <- ip$a_obs
	b_obs <- ip$b_obs

	# Estimate GI parameters from observed serial intervals
	si_all <- unlist(clusters)
	si_mean <- mean(si_all)
	si_var  <- var(si_all)
	incub_var <- ip$sd^2

	gi_mean <- max(si_mean, 0.5)  # E[SI] = E[GI] when incubation is symmetric
	gi_var_raw <- si_var - 2 * incub_var
	# Floor: if incubation variance explains most of SI variance,
	# use a conservative minimum (at least 10% of SI var, or 1 day^2)
	gi_var <- max(gi_var_raw, si_var * 0.1, 1.0)

	alpha_gi <- gi_mean^2 / gi_var
	beta_gi  <- gi_mean / gi_var

	if (gi_var_raw < 0) {
		cat(sprintf("  Note: %s GI var floored (SI var=%.1f < 2*incub var=%.1f)\n",
		    disease, si_var, 2 * incub_var))
	}

	cat(sprintf("\n  %s:\n", disease))
	cat(sprintf("    Data: %d clusters (%d with >=2 offspring), %d serial intervals\n",
	    length(clusters), length(multi_clusters), length(si_all)))
	cat(sprintf("    SI mean=%.1f, sd=%.1f | Incubation mean=%.1f, sd=%.1f\n",
	    si_mean, sd(si_all), ip$mean, ip$sd))
	cat(sprintf("    Estimated GI: mean=%.1f, var=%.1f => alpha=%.2f, beta=%.3f\n",
	    gi_mean, gi_var, alpha_gi, beta_gi))

	# Precompute density engines
	cat("    Building density engines...\n")
	engines <- precompute_engines_emp(psi_grid, alpha_gi, beta_gi, a_obs, b_obs)
	cat("    Computing posterior...\n")

	# Use ALL clusters (including singletons — they inform GI but not psi correlation)
	# The likelihood handles singletons correctly (no within-cluster correlation term)
	post <- compute_psi_posterior_emp(clusters, engines, psi_grid)

	# Summary statistics
	post_mean <- sum(psi_grid * post)
	post_cdf  <- cumsum(post)
	ci_lo <- psi_grid[which.min(abs(post_cdf - 0.025))]
	ci_hi <- psi_grid[which.min(abs(post_cdf - 0.975))]
	post_mode <- psi_grid[which.max(post)]

	cat(sprintf("    Posterior: mode=%.2f, mean=%.2f, 95%% CI=[%.2f, %.2f]\n",
	    post_mode, post_mean, ci_lo, ci_hi))

	posterior_results[[disease]] <- data.frame(
		disease    = disease,
		psi_grid   = psi_grid,
		posterior  = post,
		n_clusters = length(clusters),
		n_multi    = length(multi_clusters),
		n_si       = length(si_all),
		alpha_gi   = alpha_gi,
		beta_gi    = beta_gi,
		a_obs      = a_obs,
		b_obs      = b_obs,
		post_mode  = post_mode,
		post_mean  = post_mean,
		ci_lo      = ci_lo,
		ci_hi      = ci_hi
	)
}

# ==============================================================================
# 6. Figures
# ==============================================================================

cat("\n--- Generating figures ---\n")

posterior_df <- bind_rows(posterior_results)

# Order diseases by posterior mode of psi
disease_order <- posterior_df %>%
	group_by(disease) %>%
	summarise(mode = psi_grid[which.max(posterior)], .groups = "drop") %>%
	arrange(mode) %>%
	pull(disease)
posterior_df$disease <- factor(posterior_df$disease, levels = disease_order)

# Annotation data for summary stats
annot_df <- posterior_df %>%
	group_by(disease) %>%
	summarise(
		post_mode = psi_grid[which.max(posterior)],
		ci_lo     = first(ci_lo),
		ci_hi     = first(ci_hi),
		n_multi   = first(n_multi),
		n_si      = first(n_si),
		.groups   = "drop"
	)
annot_df$disease <- factor(annot_df$disease, levels = disease_order)
annot_df$label <- sprintf("n=%d SI, %d clusters", annot_df$n_si, annot_df$n_multi)

# --------------------------------------------------------------------------
# Figure 1: Posterior distributions for psi, one panel per disease
# --------------------------------------------------------------------------

cat("  Figure 1: Posterior distributions\n")

fig1 <- ggplot(posterior_df, aes(x = psi_grid, y = posterior)) +
	geom_line(linewidth = 0.9, color = "steelblue") +
	geom_ribbon(aes(ymin = 0, ymax = posterior), fill = "steelblue", alpha = 0.2) +
	geom_vline(data = annot_df, aes(xintercept = ci_lo),
	           linetype = "dashed", color = "grey50", linewidth = 0.4) +
	geom_vline(data = annot_df, aes(xintercept = ci_hi),
	           linetype = "dashed", color = "grey50", linewidth = 0.4) +
	geom_text(data = annot_df,
	          aes(x = 0.95, y = Inf, label = label),
	          hjust = 1, vjust = 1.5, size = 3, color = "grey40") +
	facet_wrap(~ disease, scales = "free_y", ncol = 3) +
	xlim(0, 1) +
	labs(
		x = expression(psi),
		y = "Posterior density",
		title = expression("Posterior distribution of" ~ psi ~ "from transmission tree data"),
		subtitle = "Dashed lines = 95% credible interval. Data from OutbreakTrees (Taube et al. 2022)"
	) +
	theme_minimal(base_size = 12) +
	theme(strip.text = element_text(face = "bold"))

n_diseases <- length(disease_order)
fig_height <- ceiling(n_diseases / 3) * 3.5
save_fig(fig1, "psi_empirical_posteriors", width = 14, height = fig_height)

# --------------------------------------------------------------------------
# Figure 2: Summary dot plot with credible intervals
# --------------------------------------------------------------------------

cat("  Figure 2: Summary dot plot\n")

annot_df$disease <- factor(annot_df$disease, levels = rev(disease_order))

fig2 <- ggplot(annot_df, aes(y = disease, x = post_mode)) +
	geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
	               height = 0.3, linewidth = 0.6, color = "steelblue") +
	geom_point(size = 3, color = "steelblue") +
	geom_text(aes(x = ci_hi + 0.03, label = label),
	          hjust = 0, size = 3, color = "grey40") +
	xlim(0, 1.15) +
	labs(
		x = expression(psi ~ "(posterior mode with 95% CI)"),
		y = NULL,
		title = expression("Estimated" ~ psi ~ "across diseases"),
		subtitle = "Data from OutbreakTrees (Taube et al. 2022). Bars = 95% credible intervals."
	) +
	theme_minimal(base_size = 13) +
	theme(panel.grid.major.y = element_blank())
save_fig(fig2, "psi_empirical_summary", width = 10, height = max(4, n_diseases * 0.6))

# --------------------------------------------------------------------------
# Print summary table
# --------------------------------------------------------------------------

cat("\n--- Psi estimation summary ---\n")
cat(sprintf("  %-20s  %5s  %5s  %8s  %8s  %15s  %10s\n",
    "Disease", "Multi", "TotSI", "GI_alpha", "GI_beta", "psi 95% CI", "psi mode"))
cat(paste(rep("-", 85), collapse = ""), "\n")

for (d in disease_order) {
	r <- annot_df[annot_df$disease == d, ]
	p <- posterior_results[[d]]
	cat(sprintf("  %-20s  %5d  %5d  %8.2f  %8.3f  [%4.2f, %4.2f]      %5.2f\n",
	    d, r$n_multi, r$n_si,
	    p$alpha_gi[1], p$beta_gi[1],
	    r$ci_lo, r$ci_hi, r$post_mode))
}

# Save results
results_out <- posterior_df %>%
	select(disease, psi_grid, posterior, n_clusters, n_multi, n_si,
	       alpha_gi, beta_gi, a_obs, b_obs, post_mode, post_mean, ci_lo, ci_hi)
write_csv(results_out, file.path("output", "psi_empirical_results.csv"))
cat(sprintf("\n  Saved results to output/psi_empirical_results.csv\n"))

cat("\n=== Empirical psi estimation complete ===\n")
