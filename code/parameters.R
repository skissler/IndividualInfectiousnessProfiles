parslist <- list(
	list(pathogen="influenza",
	     Tgen=3.2,
	     Tvar=2.1^2,
	     R0=2,
	     doi="10.1016/j.epidem.2025.100815"),
	list(pathogen="omicron",
	     Tgen=2.39*2.95,
	     Tvar=2.39*2.95^2,
	     R0=6,
	     doi="10.1016/j.lanepe.2022.100446"),
	list(pathogen="measles",
	     Tgen=12.2,
	     Tvar=3.62^2,
	     R0=12,
	     doi="10.1016/j.sciencedirect.com/S0022519311003146")
)

# Compute moment-matched Gamma shape and rate for each pathogen
# Gamma(alpha, beta): mean = alpha/beta, var = alpha/beta^2
#   => alpha = mean^2 / var, beta = mean / var
#
# Malthusian growth rate from the Euler-Lotka equation
#   1 = R0 * integral_0^infty exp(-r*tau) * g(tau) dtau
# with g ~ Gamma(alpha, beta), the MGF gives
#   1 = R0 * (beta / (beta + r))^alpha
#   => r = beta * (R0^(1/alpha) - 1)
#

for (i in seq_along(parslist)) {
	parslist[[i]]$alpha  <- parslist[[i]]$Tgen^2 / parslist[[i]]$Tvar
	parslist[[i]]$beta   <- parslist[[i]]$Tgen   / parslist[[i]]$Tvar
	parslist[[i]]$r      <- parslist[[i]]$beta *
	                        (parslist[[i]]$R0^(1 / parslist[[i]]$alpha) - 1)
	parslist[[i]]$p_est  <- 1 - extinction_prob(parslist[[i]]$R0)
}

# Print summary
for (p in parslist) {
	cat(sprintf("%-12s  Tgen=%.2f  Tvar=%.2f  alpha=%.2f  beta=%.3f  r=%.3f  p_est=%.3f\n",
	    p$pathogen, p$Tgen, p$Tvar, p$alpha, p$beta, p$r, p$p_est))
}

