
# PERIODIC CONTACTS 
# - Define the contact function (e.g. `make_contact_fn_periodic()`)
#    * Returns R_0 * (1 - c_amp * cos(2*pi*t/c_per)), where R_0, c_amp, and 
#      c_per are input parameters 
# - Create a grid of psi/c_amp pairs. Start with psi in [0,1] and A in [0,1]. 
#    * Start with a 20x20 grid; we can refine later. 
# - For each psi/c_amp pair, calclulate k using the theoretical result
#    * Specifically, I think this should be: 
#      k = (2 / c_amp^2) * (1 + omega^2 / beta^2)^(psi * alpha)
#      where omega = 2*pi/c_per; where alpha and beta are the parameters of the
#      generation interval distribution (Gamma(alpha, beta)); and psi is the 
#      punctuation parameter. 
# - For each psi/c_amp pair, simulate secondary infections for 10000 index
#   cases 
#    * I think this can be done by simply running 
#      gen_inf_attempts_gamma_contacts() with the appropriate inputs. 
# - Using the simulated secondary infections, calculate k empirically for each 
#   psi/c_amp pair 
#    * I think this calculation is k = m_s^2 / (v_s - m_s), where 
#      m_s is the mean number of offspring and v_s is the variance in the number
#      of offspring. This equation only holds when v_s > m_s; otherwise, it 
#      should return Inf. 
# - Generate a heat map with the empirical k values for each psi/c_amp pair as 
#   the patches. Overlay contours reflecting the theoretical k values across 
#   psi/c_amp. 
# - Simulate 5,000 epidemics in 10,000 people for psi=0, psi=0.5, and psi=1, 
#   holding contact amplitude fixed at 0.5. Calculate the number of epidemics 
#   that don't establish, i.e., that fail to reach 500 infections. Report these 
#   failure proportions in a table. 
#    * I think this can be done using sim_stochastic_fast() with
#      gen_inf_attempts_gamma_contacts() and the appropriate inputs. 

# GAMMA CONTACTS 
# - Define the contact function (e.g. `make_contact_fn_gamma()`)
#    * This is genuinely new. See lines 87-94 of writeup/model.tex for a 
#      description of this contact model. Basically, for each person, contact 
#      levels change at arrivals of a Poisson process with rate lambda, so that
#      the amount of time spent at each level is exponentially distributed with
#      mean 1/lambda. The contact levels are drawn independently from a 
#      Gamma(k_c, k_c) distribution. As before, if a level is 
#      l ~ Gamma(k_c, k_c), then the function should output R_0*l for all the 
#      levels l. 
# - Create a grid of psi/c_amp pairs. Start with psi in [0,1] and k_c in 
#   [0.1,1000], using a log scale for k_c. 
#    * Start with a 20x20 grid; we can refine later. 
# - For each psi/c_amp pair, simulate secondary infections for 10000 index 
#   cases 
# - Using the simulated secondary infections, calculate k empirically for each 
#   psi/c_amp pair 
# - Generate a heat map with the empirical k values for each psi/c_amp pair as 
#   the patches. 
# - Simulate 5,000 epidemics in 10,000 people for psi=0, psi=0.5, and psi=1, 
#   holding contact amplitude fixed at 0.5. Calculate the number of epidemics 
#   that don't establish, i.e., that fail to reach 500 infections. Report these 
#   failure proportions in a table. 


