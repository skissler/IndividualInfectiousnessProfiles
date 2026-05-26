library(tidyverse)
library(odin)
library(Rcpp)

# Compile Rcpp simulation engine (priority-queue based, ~50-100x faster)
sourceCpp("code/src/sim_stochastic_rcpp.cpp")

# ==============================================================================
# Deterministic models 
# ==============================================================================

#' SEIR ODE model
#'
#' Generates a standard SEIR epidemic trajectory, initialized with a single
#' person in the infectious (I) state. 
#'
#' @param R0    Basic reproduction number (default 2)
#' @param nu    E->I rate, 1/mean latent period (default 1/2)
#' @param gamma I->R rate, 1/mean infectious period (default 1/3)
#' @param N     Population size, used to set initial conditions (default 1000)
seir <- odin({
  ## Derivatives
  deriv(S) <- -beta*S*I
  deriv(E) <- beta*S*I - nu*E
  deriv(I) <- nu*E - gamma*I
  deriv(R) <- gamma*I
  deriv(cuminf) <- beta*S*I

  ## Initial conditions
  initial(S) <- 1-1/N
  initial(E) <- 0
  initial(I) <- 1/N
  initial(R) <- 0
  initial(cuminf) <- 1/N

  ## Parameters with default values 
  R0 <- user(2)
  nu <- user(1/2)
  gamma <- user(1/3)
  beta <- R0*gamma
  N <- user(1000)
})

#' Renewal equation epidemic model
#'
#' Generates a deterministic epidemic trajectory using a discrete approximation
#' to the renewal equation framework, with a Gamma-distributed infectiousness 
#' profile and initalized with a single infectious person. 
#'
#' @param R0     Basic reproduction number
#' @param alpha  Shape parameter of Gamma generation interval distribution
#' @param Tgen  Mean generation interval (days)
#' @param N      Population size, used to set initial conditions
#' @param dt     Time step for discretization (default 0.01 days)
#' @param tmax   Maximum simulation length (default 100 days)

renewal_epidemic <- function(R0, alpha, Tgen, N, dt = 0.01, tmax = 100) {

	# Set rate parameter of the generation interval distribution
  beta <- alpha / Tgen

  # Initalize time vectors
  times <- seq(0, tmax, by = dt)
  nt <- length(times)

  # Initialize incidence and cumulative infections
  inc <- numeric(nt)
  cuminf <- numeric(nt)
  cuminf[1] <- 1 / N

  # Generation interval weights (with dt baked in for the convolution)
  g <- dgamma(times, shape = alpha, rate = beta) * dt
  # Index case kernel (no dt — it's a single person, not a rate)
  g_index <- dgamma(times, shape = alpha, rate = beta) / N

  for (i in 2:nt) {
    S <- 1 - cuminf[i-1]
    # Index case + renewal convolution
    inc[i] <- S * R0 * (g_index[i] + sum(g[2:i] * inc[(i-1):1]))
    cuminf[i] <- cuminf[i-1] + inc[i] * dt
  }

  tibble(t = times, inc = inc, cuminf = cuminf, S = 1 - cuminf)
}

# ==============================================================================
# Stochastic simulation model
# ==============================================================================

#' Infection attempt generator: stepwise SEIR profile
#'
#' Returns a closure that, given an infection time, draws a latent period
#' E_j ~ Exp(1/e_dur) and infectious period I_j ~ Exp(1/i_dur), then generates
#' Poisson(nu_j) attempts uniformly distributed over the infectious window 
#' [E_j, E_j+I_j], where nu_j = beta*I_j. Corresponds to the classical SEIR
#' compartmental model.
#'
#' @param e_dur Mean latent (E) period (days).
#' @param i_dur Mean infectious (I) period (days).
#' @param R0 Basic reproduction number (expected number of infection attempts).
#' @return Function(tinf) returning a numeric vector of absolute attempt times.
gen_inf_attempts_stepwise <- function(e_dur, i_dur, R0){
	beta <- R0/i_dur
	function(tinf){
		t1 <- rexp(1, 1/e_dur)
		t2 <- t1 + rexp(1, 1/i_dur)
		nattempts <- rpois(1, (t2 - t1)*beta) 
		if(nattempts == 0L) return(numeric(0))
		attempt_times <- tinf + sort(runif(nattempts, min=t1, max=t2))
		return(attempt_times)
	}
}

#' Infection attempt generator: smooth SEIR profile
#'
#' Returns a closure that draws Poisson(R0) attempts, each at an independent
#' time E_j + I_j where E_j ~ Exp(1/e_dur) and I_j ~ Exp(1/i_dur). All
#' offspring times are iid draws from the generation interval distribution
#' (convolution of two exponentials). 
#'
#' @param e_dur Mean latent (E) period (days).
#' @param i_dur Mean infectious (I) period (days).
#' @param R0 Basic reproduction number (expected number of infection attempts).
#' @return Function(tinf) returning a numeric vector of absolute attempt times.
gen_inf_attempts_smooth <- function(e_dur, i_dur, R0){
	function(tinf){
		nattempts <- rpois(1,R0)
		if(nattempts == 0L) return(numeric(0))
		attempt_times <- tinf + sort(rexp(nattempts, 1/e_dur) + rexp(nattempts, 1/i_dur))
		return(attempt_times)
	}
}

#' Infection attempt generator: spike SEIR profile
#'
#' Returns a closure that draws Poisson(R0) attempts, all occurring at the
#' same instant E_j + I_j, where E_j ~ Exp(1/e_dur) and I_j ~ Exp(1/i_dur) 
#' are drawn once and shared across all offspring. 
#'
#' @param e_dur Mean latent (E) period (days).
#' @param i_dur Mean infectious (I) period (days).
#' @param R0 Basic reproduction number (expected number of infection attempts).
#' @return Function(tinf) returning a numeric vector of absolute attempt times.
gen_inf_attempts_spike <- function(e_dur, i_dur, R0){
	function(tinf){
		nattempts <- rpois(1,R0)
		if(nattempts == 0L) return(numeric(0))
		attempt_times <- tinf + rep(rexp(1, 1/e_dur) + rexp(1, 1/i_dur), nattempts)
		return(attempt_times)
	}
}

#' Raw infection attempt generator (internal helper)
#'
#' Creates a closure that generates raw (unsorted) infection attempt times and
#' the mode time for a single person. The three-branch logic (spike / smooth /
#' general psi) lives here once; the public gen_inf_attempts_gamma* functions
#' wrap this with their own post-processing.
#'
#' @param n_rate Poisson rate for the number of proposals (R0 or z_max).
#' @param alpha Shape parameter for the population generation interval.
#' @param beta Rate parameter for the population generation interval.
#' @param psi Punctuation parameter (0 = spike, 1 = smooth).
#' @return A function(tinf) returning list(attempts, mode_time).
.make_raw_generator <- function(n_rate, alpha, beta, psi) {
	if (psi < 1e-6) { # spike: all proposals at the same instant
		function(tinf) {
			n <- rpois(1L, n_rate)
			if (n == 0L) return(list(attempts = numeric(0), mode_time = NA_real_))
			shift <- rgamma(1, shape = alpha, rate = beta)
			list(attempts = rep(tinf + shift, n), mode_time = tinf + shift)
		}
	} else if (psi > 1 - 1e-6) { # smooth: independent generation intervals
		mode_psi <- max(0, (alpha - 1) / beta)
		function(tinf) {
			n <- rpois(1L, n_rate)
			if (n == 0L) return(list(attempts = numeric(0), mode_time = NA_real_))
			list(attempts = tinf + rgamma(n, shape = alpha, rate = beta),
			     mode_time = tinf + mode_psi)
		}
	} else { # general: shared shift + independent jitter
		shiftshape <- (1 - psi) * alpha
		mode_psi <- max(0, (psi * alpha - 1) / beta)
		function(tinf) {
			n <- rpois(1L, n_rate)
			if (n == 0L) return(list(attempts = numeric(0), mode_time = NA_real_))
			shift <- rgamma(1, shape = shiftshape, rate = beta)
			list(attempts = tinf + shift + rgamma(n, shape = psi * alpha, rate = beta),
			     mode_time = tinf + shift + mode_psi)
		}
	}
}

#' Infection attempt generator for gamma profiles
#'
#' Generates attempted infection times for a single person, themselves infected
#' at time tinf, using the Gamma-distributed individual infectiousness profile
#' with punctuation parameter psi.
#'
#' @param Tgen Mean generation interval
#' @param R0 Basic reproduction number
#' @param alpha Shape parameter for the Gamma-distributed population
#'   infectiousness profile (A(tau) ~ Gamma(alpha, beta) where
#'   beta = alpha/Tgen)
#' @param psi Parameter governing punctuation of the individual
#'   infectiousness profile (psi \in [0, 1], with psi -> 0 giving sharp
#'   infectiousness profiles and psi -> 1 giving smooth ones equivalent to
#'   the population generation interval distribution)
#' @return A function(t_inf) that returns sorted infection attempt times
gen_inf_attempts_gamma <- function(Tgen, R0, alpha, psi) {
	stopifnot(psi >= 0, psi <= 1)
	raw <- .make_raw_generator(R0, alpha, alpha / Tgen, psi)
	function(tinf) sort(raw(tinf)$attempts)
}

#' Infection attempt generator: gamma profile with time-varying deterministic #' contact rate
#'
#' Like \code{gen_inf_attempts_gamma}, but modulates infection attempts by a
#' time-varying contact rate z(t). Uses thinning: proposes Poisson(z_max) 
#' attempts from the gamma profile, then independently accepts each with 
#' probability z(t)/z_max, where t is the absolute time of the proposal.
#'
#' @param Tgen Mean generation interval (days).
#' @param z Function of absolute time returning the contact rate (must satisfy
#'   0 <= z(t) <= z_max for all t).
#' @param z_max Upper bound on z(t), used as the proposal rate.
#' @param alpha Shape parameter for the population generation interval
#'   (Gamma(alpha, alpha/Tgen)).
#' @param psi Punctuation parameter (psi in [0, 1]).
#' @return Function(tinf) returning a numeric vector of absolute attempt times.
gen_inf_attempts_gamma_contacts <- function(Tgen, z, z_max, alpha, psi) {
	stopifnot(psi >= 0, psi <= 1, z_max > 0)
	raw <- .make_raw_generator(z_max, alpha, alpha / Tgen, psi)
	function(tinf) {
		res <- raw(tinf)
		if (length(res$attempts) == 0L) return(numeric(0))
		keep <- runif(length(res$attempts)) < z(res$attempts) / z_max
		if (!any(keep)) return(numeric(0))
		sort(res$attempts[keep])
	}
}

#' Truncated Gamma sampler (rejection method)
#'
#' Draws n samples from Gamma(shape, rate) truncated above at upper.
#' Uses simple rejection sampling. If upper = Inf, falls back to rgamma().
#'
#' @param n Number of draws.
#' @param shape Shape parameter of the Gamma distribution.
#' @param rate Rate parameter of the Gamma distribution.
#' @param upper Upper truncation point (default Inf = no truncation).
#' @return Numeric vector of length n.
rtgamma <- function(n, shape, rate, upper = Inf) {
	stopifnot(upper>=0)
	if (is.infinite(upper)) return(rgamma(n, shape, rate))
	x <- numeric(n)
	for (i in seq_len(n)) {
		repeat {
			draw <- rgamma(1, shape, rate)
			if (draw <= upper) { x[i] <- draw; break }
		}
	}
	x
}

#' Infection attempt generator: gamma profile with Gamma/Poisson stochastic
#' contact rate
#'
#' Like \code{gen_inf_attempts_gamma_contacts}, but generates a fresh random
#' piecewise-constant contact trajectory per person. Contact levels are iid
#' Gamma(k_c, k_c) and switch at Poisson(lambda) arrivals.
#'
#' @param Tgen Mean generation interval (days).
#' @param R0 Basic reproduction number.
#' @param alpha Shape parameter for the population generation interval.
#' @param psi Punctuation parameter (psi in [0, 1]).
#' @param k_c Shape (and rate) of Gamma distribution for contact levels.
#' @param lambda Rate of the Poisson switching process (switches per day).
#' @param c_max Upper truncation point for contact levels (gathering size cap).
#'   Default Inf (no truncation). When finite, contact levels are drawn from
#'   Gamma(k_c, k_c) truncated at c_max, reducing R0 by a factor mu_T.
#' @return Function(tinf) returning a numeric vector of absolute attempt times.
gen_inf_attempts_gammapoisson_contacts <- function(Tgen, R0, alpha, psi, k_c, lambda, c_max = Inf) {
	stopifnot(psi >= 0, psi <= 1, k_c > 0, lambda > 0)
	beta <- alpha / Tgen
	traj_duration <- 5 * Tgen  # covers >99.99% of Gamma(alpha, beta) mass

	function(tinf) {
		# Generate per-person piecewise-constant contact trajectory
		n_switches <- rpois(1, lambda * traj_duration)
		switch_offsets <- if (n_switches > 0) sort(runif(n_switches, 0, traj_duration)) else numeric(0)
		levels <- rtgamma(n_switches + 1, k_c, k_c, upper = c_max)
		breaks <- tinf + c(0, switch_offsets)
		z_max_person <- R0 * max(levels)

		# Generate raw proposals at this person's z_max
		raw <- .make_raw_generator(z_max_person, alpha, beta, psi)
		res <- raw(tinf)
		if (length(res$attempts) == 0L) return(numeric(0))

		# Thin using per-person contact trajectory
		seg_idx <- findInterval(res$attempts, breaks)
		z_vals <- R0 * levels[seg_idx]
		keep <- runif(length(res$attempts)) < z_vals / z_max_person
		if (!any(keep)) return(numeric(0))
		sort(res$attempts[keep])
	}
}


#' Mean-reduction factor for truncated Gamma contact levels
#'
#' Computes mu_T = E[X | X <= c_max] / E[X] for X ~ Gamma(k_c, k_c),
#' where E[X] = 1. Uses the identity E[X | X <= c] = (k_c/k_c) *
#' pgamma(c, k_c+1, k_c) / pgamma(c, k_c, k_c).
#'
#' @param k_c Shape (and rate) of the Gamma contact-level distribution.
#' @param c_max Upper truncation point (gathering size cap).
#' @return Scalar mu_T in (0, 1].
mu_truncated_gamma <- function(k_c, c_max) {
	if (is.infinite(c_max)) return(1)
	pgamma(c_max, k_c + 1, k_c) / pgamma(c_max, k_c, k_c)
}

# --- Helpers for detection-and-isolation protocols ---

#' Determine isolation time under regular screening
#'
#' Draws a random test phase uniformly in [0, Delta), then tests at regular
#' intervals within the detectability window [mode_time - d_pre, mode_time -
#' d_pre + w]. Each test detects with probability p_sens. If detected, isolation
#' occurs after an Exp(lambda_act) delay. Returns Inf if undetected.
#'
#' @param mode_time Absolute time of the infectiousness profile mode.
#' @param d_pre Days before mode that the detectability window begins.
#' @param w Total width of the detectability window (d_pre + d_post).
#' @param Delta Gap between tests (days).
#' @param p_sens Test sensitivity (probability of detection per test).
#' @param lambda_act Rate of exponential delay from detection to isolation.
#' @return Scalar: absolute isolation time, or Inf if undetected.
.detect_screening <- function(mode_time, d_pre, w, Delta, p_sens, lambda_act){
	window_start <- mode_time - d_pre
	window_end <- window_start + w
	U <- runif(1, 0, Delta)
	first_test <- window_start + U
	if(first_test > window_end) return(Inf)
	test_times <- seq(first_test, window_end, by=Delta)
	if(p_sens >= 1){
		tau_det <- test_times[1]
	} else {
		detected <- runif(length(test_times)) < p_sens
		if(!any(detected)) return(Inf)
		tau_det <- test_times[which(detected)[1]]
	}
	tau_det + rexp(1, lambda_act)
}

#' Determine isolation time under symptom-triggered isolation
#'
#' With probability p_sym the individual is symptomatic; symptom onset occurs
#' at mode_time + N(mu_sym, sigma_sym). Isolation follows after an
#' Exp(lambda_act) delay. Returns Inf if asymptomatic.
#'
#' @param mode_time Absolute time of the infectiousness profile mode.
#' @param mu_sym Mean symptom onset time relative to mode (days).
#' @param sigma_sym SD of symptom onset time (0 = deterministic).
#' @param lambda_act Rate of exponential delay from symptom onset to isolation.
#' @param p_sym Probability of being symptomatic.
#' @return Scalar: absolute isolation time, or Inf if asymptomatic.
.detect_symptoms <- function(mode_time, mu_sym, sigma_sym, lambda_act, p_sym){
	if(runif(1) > p_sym) return(Inf)
	delta_sym <- if(sigma_sym > 0) rnorm(1, mu_sym, sigma_sym) else mu_sym
	mode_time + delta_sym + rexp(1, lambda_act)
}

#' Remove or thin infection attempts after isolation
#'
#' Keeps all attempts at or before tau_iso. If eta = 1, discards all later
#' attempts (perfect isolation). If 0 < eta < 1, each post-isolation attempt
#' is independently removed with probability eta (imperfect isolation).
#'
#' @param attempt_times Numeric vector of absolute infection attempt times.
#' @param tau_iso Absolute isolation time (Inf if not isolated).
#' @param eta Isolation effectiveness: probability of blocking each
#'   post-isolation attempt (1 = perfect, 0 = no effect).
#' @return Sorted numeric vector of surviving attempt times.
.thin_attempts <- function(attempt_times, tau_iso, eta){
	if(!is.finite(tau_iso) || length(attempt_times) == 0L) return(attempt_times)
	pre <- attempt_times <= tau_iso
	if(eta >= 1){
		attempt_times <- attempt_times[pre]
	} else {
		npost <- sum(!pre)
		if(npost > 0L){
			post_keep <- runif(npost) >= eta
			attempt_times <- c(attempt_times[pre], attempt_times[!pre][post_keep])
		}
	}
	if(length(attempt_times) == 0L) return(numeric(0))
	sort(attempt_times)
}

#' Infection attempt generator with regular screening and isolation
#'
#' Like gen_inf_attempts_gamma, but after generating attempts, determines
#' whether the individual is detected via regular screening (tests every
#' Delta days with random phase, detection within a window around peak
#' infectiousness) and removes post-isolation attempts.
#'
#' @param Tgen Mean generation interval
#' @param R0 Basic reproduction number
#' @param alpha Shape parameter for population infectiousness profile
#' @param psi Punctuation parameter (0 = spike, 1 = smooth)
#' @param d_pre Days before peak that detectability window begins
#' @param d_post Days after peak that detectability window ends
#' @param Delta Gap between tests (days)
#' @param lambda_act Rate of exponential delay from detection to isolation
#'   (Inf = immediate isolation, the default)
#' @param p_sens Test sensitivity (probability of detection per test)
#' @param eta Isolation effectiveness (1 = perfect, 0 = none)
#' @return A function(t_inf) that returns sorted infection attempt times
gen_inf_attempts_gamma_screening <- function(Tgen, R0, alpha, psi,
                                             d_pre, d_post, Delta,
                                             lambda_act=Inf, p_sens=1, eta=1) {
	stopifnot(psi >= 0, psi <= 1)
	raw <- .make_raw_generator(R0, alpha, alpha / Tgen, psi)
	w <- d_pre + d_post
	function(tinf) {
		res <- raw(tinf)
		if (length(res$attempts) == 0L) return(numeric(0))
		tau_iso <- .detect_screening(res$mode_time, d_pre, w, Delta, p_sens, lambda_act)
		.thin_attempts(res$attempts, tau_iso, eta)
	}
}

#' Infection attempt generator with symptom-triggered isolation
#'
#' Like gen_inf_attempts_gamma, but after generating attempts, determines
#' whether the individual develops symptoms (probability p_sym) at a time
#' anchored to peak infectiousness, and removes post-isolation attempts.
#' Subsumes the fixed-isolation case: use sigma_sym=0, lambda_act=Inf,
#' mu_sym=tau_offset, p_sym=p_adhere for deterministic isolation at a fixed
#' offset relative to peak.
#'
#' @param Tgen Mean generation interval
#' @param R0 Basic reproduction number
#' @param alpha Shape parameter for population infectiousness profile
#' @param psi Punctuation parameter (0 = spike, 1 = smooth)
#' @param mu_sym Mean symptom onset time relative to peak infectiousness
#'   (positive = symptoms lag peak; negative = symptoms precede peak)
#' @param sigma_sym SD of symptom onset time (0 = deterministic, the default)
#' @param lambda_act Rate of exponential delay from symptom onset to isolation
#'   (Inf = immediate isolation, the default)
#' @param p_sym Probability of being symptomatic
#' @param eta Isolation effectiveness (1 = perfect, 0 = none)
#' @return A function(t_inf) that returns sorted infection attempt times
gen_inf_attempts_gamma_symptoms <- function(Tgen, R0, alpha, psi,
                                            mu_sym, sigma_sym=0,
                                            lambda_act=Inf, p_sym=1, eta=1) {
	stopifnot(psi >= 0, psi <= 1)
	raw <- .make_raw_generator(R0, alpha, alpha / Tgen, psi)
	function(tinf) {
		res <- raw(tinf)
		if (length(res$attempts) == 0L) return(numeric(0))
		tau_iso <- .detect_symptoms(res$mode_time, mu_sym, sigma_sym, lambda_act, p_sym)
		.thin_attempts(res$attempts, tau_iso, eta)
	}
}

#' Stochastic epidemic simulation (R fallback)
#'
#' Pure-R implementation using sorted vector queue. Kept as fallback;
#' the default sim_stochastic_fast() uses the Rcpp priority-queue version.
#'
#' @param n Population size (default 1000)
#' @param gen_inf_attempts Infection attempt time generator function,
#'   capturing the individual infectiousness profile, with signature
#'   function(t_inf) -> numeric vector.
#' @param maxinf Max number of infections to simulate (for early stopping)
#' @return Numeric vector of length n containing infection times (Inf if the
#'   person remains uninfected at the end of the simulation)
sim_stochastic_fast_r <- function(n=1000, gen_inf_attempts, maxinf=Inf){

	# Initialize infection times
	tinf_vec <- rep(Inf, n)

	# Seed infection
	indexcase <- sample.int(n,1)
	tinf_vec[indexcase] <- 0
	n_infected <- 1L
	queue <- gen_inf_attempts(0)
	qi <- 1L
	qn <- length(queue)

	# Process events in chronological order
	while(qi <= qn){
		t_attempt <- queue[qi]
		qi <- qi + 1L

		target <- sample.int(n,1)
		if(tinf_vec[target]==Inf){
			tinf_vec[target] <- t_attempt
			n_infected <- n_infected + 1L
			if(n_infected >= maxinf) return(tinf_vec)
			new_attempts <- gen_inf_attempts(t_attempt)
			if(length(new_attempts) > 0L){
				# Merge new events with unprocessed remainder of queue
				remaining <- if(qi <= qn) queue[qi:qn] else numeric(0)
				queue <- sort.int(c(remaining, new_attempts))
				qi <- 1L
				qn <- length(queue)
			}
		}
	}
	return(tinf_vec)
}

#' Stochastic epidemic simulation
#'
#' Generates the infection times for a population (size n) given individual
#' infectiousness profiles specified by gen_inf_attempts. Uses the Rcpp
#' priority-queue implementation for speed.
#'
#' @param n Population size (default 1000)
#' @param gen_inf_attempts Infection attempt time generator function,
#'   capturing the individual infectiousness profile, with signature
#'   function(t_inf) -> numeric vector.
#' @param maxinf Max number of infections to simulate (for early stopping)
#' @return Numeric vector of length n containing infection times (Inf if the
#'   person remains uninfected at the end of the simulation)
sim_stochastic_fast <- function(n=1000, gen_inf_attempts, maxinf=Inf){
	sim_stochastic_rcpp(n, gen_inf_attempts, as.integer(if(is.infinite(maxinf)) -1L else maxinf))
}

#' Stochastic epidemic simulation with infinite population (R fallback)
#'
#' Like sim_stochastic_fast_r, but assumes an infinite susceptible population
#' so every infection attempt succeeds. Useful for estimating exponential
#' growth rates without susceptible depletion effects.
#'
#' @param max_cases Maximum number of infections before stopping (default 1000)
#' @param tmax Maximum simulation time; stops processing attempts beyond this
#'   (default Inf)
#' @param gen_inf_attempts Infection attempt time generator function,
#'   capturing the individual infectiousness profile, with signature
#'   function(t_inf) -> numeric vector.
#' @return Sorted numeric vector of infection times (length <= max_cases),
#'   starting with 0 for the index case
sim_infinite_pop_r <- function(max_cases=1000, tmax=Inf, gen_inf_attempts){

	tinf_vec <- numeric(max_cases)
	tinf_vec[1] <- 0  # index case
	n_infected <- 1L

	queue <- gen_inf_attempts(0)
	qi <- 1L
	qn <- length(queue)

	while(qi <= qn && n_infected < max_cases){
		t_attempt <- queue[qi]
		qi <- qi + 1L

		if(t_attempt > tmax) break

		# Every attempt succeeds (infinite susceptible population)
		n_infected <- n_infected + 1L
		tinf_vec[n_infected] <- t_attempt

		# Only generate new attempts if we still need more cases
		if(n_infected < max_cases){
			new_attempts <- gen_inf_attempts(t_attempt)
			if(length(new_attempts) > 0L){
				remaining <- if(qi <= qn) queue[qi:qn] else numeric(0)
				queue <- sort.int(c(remaining, new_attempts))
				qi <- 1L
				qn <- length(queue)
			}
		}
	}

	return(tinf_vec[1:n_infected])
}

#' Stochastic epidemic simulation with infinite population
#'
#' Like sim_stochastic_fast, but assumes an infinite susceptible population
#' so every infection attempt succeeds. Uses the Rcpp priority-queue
#' implementation for speed.
#'
#' @param max_cases Maximum number of infections before stopping (default 1000)
#' @param tmax Maximum simulation time (default Inf)
#' @param gen_inf_attempts Infection attempt time generator function,
#'   capturing the individual infectiousness profile, with signature
#'   function(t_inf) -> numeric vector.
#' @return Sorted numeric vector of infection times (length <= max_cases),
#'   starting with 0 for the index case
sim_infinite_pop <- function(max_cases=1000, tmax=Inf, gen_inf_attempts){
	sim_infinite_pop_rcpp(max_cases, tmax, gen_inf_attempts)
}

#' Overlap integral I* for the Gamma/Poisson contact model
#'
#' Computes I*(psi*alpha, lambda/beta) via numerical integration of the
#' spectral formula. Used to obtain the theoretical dispersion parameter
#' k = k_c / I*.
#'
#' @param psi_alpha Product of punctuation parameter and generation-interval shape.
#' @param r Ratio lambda/beta (switching rate / generation-interval rate).
#' @return Scalar value of I* in (0, 1].
.I_star <- function(psi_alpha, r) {
	if (psi_alpha < 1e-10) return(1)  # spike limit: I* = 1
	integrand <- function(omega) {
		(1 + omega^2)^(-psi_alpha) * (r / (r^2 + omega^2)) / pi
	}
	2 * integrate(integrand, 0, Inf, rel.tol = 1e-10)$value
}

#' Theoretical dispersion parameter for Gamma/Poisson contact model
#'
#' Computes k = k_c / I*(psi*alpha, lambda/beta) using the spectral overlap
#' integral.
#'
#' @param k_c Shape (and rate) of Gamma contact-level distribution.
#' @param psi Punctuation parameter.
#' @param alpha Shape parameter for the population generation interval.
#' @param lambda Poisson switching rate.
#' @param beta Rate parameter for the population generation interval (alpha/Tgen).
#' @return Scalar theoretical dispersion parameter k.
k_theory_gammapoisson <- function(k_c, psi, alpha, lambda, beta) {
	k_c / .I_star(psi * alpha, lambda / beta)
}

#' Extinction probability for a Poisson(R0) branching process
#'
#' Computes the smallest root of q = exp(R0*(q-1)) by fixed-point iteration.
#' Returns 1 if R0 <= 1 (certain extinction).
#'
#' @param R0 Basic reproduction number.
#' @return Scalar extinction probability in [0, 1].
extinction_prob <- function(R0) {
	if (R0 <= 1) return(1)
	q <- 0.5
	for (i in 1:1000) {
		q_new <- exp(R0 * (q - 1))
		if (abs(q_new - q) < 1e-12) break
		q <- q_new
	}
	q
}

# Helper: save a ggplot as both .png and .pdf in figures/
save_fig <- function(plot, name, width = 10, height = 5) {
	ggsave(file.path("figures", paste0(name, ".pdf")), plot, width = width, height = height)
	ggsave(file.path("figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 300)
}

# ==============================================================================
# Cache helpers: summary + plot trajectories (two small files)
# ==============================================================================

# Helper: construct cache file paths
cache_path_summary <- function(pathogen, nsim, popsize) {
	file.path("output", sprintf("sim_summary_%s_n%d_s%d.csv", pathogen, popsize, nsim))
}

cache_path_plot <- function(pathogen, nsim, popsize) {
	file.path("output", sprintf("plot_trajectories_%s_n%d_s%d.csv", pathogen, popsize, nsim))
}

# Helper: load cached simulations, or return NULL if missing/stale
#
# Returns a list with $summary (one row per sim*psi) and $plot (individual
# infection times for first max_plot_sims simulations). Returns NULL if
# cache is missing or doesn't contain all required psi values.
load_cache <- function(pathogen, nsim, popsize, psivals) {
	summary_file <- cache_path_summary(pathogen, nsim, popsize)
	plot_file    <- cache_path_plot(pathogen, nsim, popsize)

	if (!file.exists(summary_file) || !file.exists(plot_file)) return(NULL)

	summary_df <- read_csv(summary_file, show_col_types = FALSE)
	cached_psi <- sort(unique(summary_df$psi))
	if (!all(psivals %in% cached_psi)) {
		cat(sprintf("  %s: cache missing psi values %s, re-running\n",
		    pathogen, paste(setdiff(psivals, cached_psi), collapse = ", ")))
		return(NULL)
	}

	plot_df <- read_csv(plot_file, show_col_types = FALSE)

	cat(sprintf("  %s: loading cached simulations from %s\n", pathogen, summary_file))
	list(
		summary = summary_df %>% mutate(psi = factor(psi)),
		plot    = plot_df    %>% mutate(psi = factor(psi))
	)
}

# ==============================================================================
# Infinite-population cache helpers
# ==============================================================================

cache_path_infpop_summary <- function(pathogen, nsim, max_cases) {
	file.path("output", sprintf("infpop_summary_%s_m%d_s%d.csv", pathogen, max_cases, nsim))
}

cache_path_infpop_plot <- function(pathogen, nsim, max_cases) {
	file.path("output", sprintf("infpop_plot_%s_m%d_s%d.csv", pathogen, max_cases, nsim))
}

load_cache_infpop <- function(pathogen, nsim, max_cases, psivals) {
	summary_file <- cache_path_infpop_summary(pathogen, nsim, max_cases)
	plot_file    <- cache_path_infpop_plot(pathogen, nsim, max_cases)

	if (!file.exists(summary_file) || !file.exists(plot_file)) return(NULL)

	summary_df <- read_csv(summary_file, show_col_types = FALSE)
	cached_psi <- sort(unique(summary_df$psi))
	if (!all(psivals %in% cached_psi)) {
		cat(sprintf("  %s: infpop cache missing psi values %s, re-running\n",
		    pathogen, paste(setdiff(psivals, cached_psi), collapse = ", ")))
		return(NULL)
	}

	plot_df <- read_csv(plot_file, show_col_types = FALSE)

	cat(sprintf("  %s: loading cached infpop simulations from %s\n", pathogen, summary_file))
	list(
		summary = summary_df %>% mutate(psi = factor(psi)),
		plot    = plot_df    %>% mutate(psi = factor(psi))
	)
}

# ==============================================================================
# Growth rate helper
# ==============================================================================

#' Compute epidemic growth rate from infection times via Poisson GLM
#'
#' Given a sorted vector of infection times, identifies the first full
#' calendar day after the cumulative case count reaches min_threshold,
#' then counts ALL infections on each of the next window_days calendar
#' days and fits a log-linear Poisson GLM to estimate the exponential
#' growth rate.
#'
#' Because the window is defined by complete calendar days (not by
#' cumulative case thresholds), there is no boundary truncation and no
#' need to drop the first or last day.
#'
#' @param infection_times Sorted vector of infection times (finite values only)
#' @param min_threshold  Cumulative case count that marks establishment (default 100)
#' @param window_days    Number of full calendar days in the estimation window (default 7)
#' @return Scalar growth rate (coefficient on day), or NA if too few data points
compute_growth_rate <- function(infection_times, min_threshold = 100, window_days = 7) {
	stopifnot(all(is.finite(infection_times)), !is.unsorted(infection_times))
	n <- length(infection_times)
	if (n < min_threshold) return(NA_real_)

	# First full calendar day after reaching min_threshold cases
	start_day <- floor(infection_times[min_threshold]) + 1L
	end_day   <- start_day + window_days - 1L

	# Check if simulation data covers the full window
	if (floor(infection_times[n]) < end_day) {
		warning(sprintf(
			"compute_growth_rate: simulation ends on day %d but window extends to day %d (max_cases may be too low)",
			floor(infection_times[n]), end_day), call. = FALSE)
		return(NA_real_)
	}

	# Count ALL infections on each day of the window
	days <- floor(infection_times)
	day_seq <- start_day:end_day
	day_counts <- table(factor(days[days >= start_day & days <= end_day],
	                           levels = day_seq))
	counts <- as.integer(day_counts)
	day0 <- seq(0L, window_days - 1L)

	if (sum(counts > 0) < 3) return(NA_real_)

	tryCatch(
		coef(glm(counts ~ day0, family = poisson))[2],
		error = function(e) NA_real_
	)
}

#' Simulate transmission clusters
#'
#' @param K Number of index cases (clusters)
#' @param psi_true True punctuation parameter
#' @param R0 Basic reproduction number
#' @param alpha,beta Generation interval Gamma parameters
#' @param a_obs,b_obs Detection delay Gamma parameters
#' @param p_asc Ascertainment probability
#' @return List of numeric vectors, one per informative cluster (m >= 2),
#'   containing observed serial intervals
simulate_clusters <- function(K, psi_true, R0, alpha, beta, a_obs, b_obs, p_asc=1) {

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