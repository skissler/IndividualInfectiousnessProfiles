// psiinference.stan
// Full latent variable model for inferring the punctuation parameter psi
// from transmission cluster serial interval data.
//
// Generative model for cluster k with m_k siblings:
//   l_k        ~ Gamma((1-psi)*alpha, beta)    shared latent period
//   d0_k       ~ Gamma(a_obs, b_obs)            index case detection delay
//   eps_{k,j}  ~ Gamma(psi*alpha, beta)          individual jitter
//   d_{k,j}    ~ Gamma(a_obs, b_obs)             sibling detection delay
//   s_{k,j}    = l_k + eps_{k,j} + d_{k,j} - d0_k   observed serial interval
//
// We observe s and infer psi. The free latent variables are l, d0, and eps.
// The sibling detection delay d is determined by the constraint:
//   d_{k,j} = s_{k,j} - l_k - eps_{k,j} + d0_k
// and must be positive for the model to be valid.

data {
  int<lower=1> K;                  // number of clusters
  array[K] int<lower=2> m;        // number of siblings per cluster
  int<lower=1> N;                  // total number of siblings (sum of m)
  vector[N] s;                     // observed serial intervals
  array[N] int<lower=1,upper=K> cluster;  // cluster membership index

  // Known parameters
  real<lower=0> alpha;             // generation interval Gamma shape
  real<lower=0> beta_gi;           // generation interval Gamma rate
  real<lower=0> a_obs;             // detection delay Gamma shape
  real<lower=0> b_obs;             // detection delay Gamma rate
}

parameters {
  real<lower=0, upper=1> psi;      // punctuation parameter
  vector<lower=0>[K] l;            // shared latent period per cluster
  vector<lower=0>[K] d0;           // index case detection delay per cluster
  vector<lower=0>[N] eps;          // individual jitter per sibling
}

transformed parameters {
  // Sibling detection delay, derived from the observation constraint
  vector[N] d;
  for (n in 1:N)
    d[n] = s[n] - l[cluster[n]] - eps[n] + d0[cluster[n]];
}

model {
  // Shape parameters, clamped to avoid degenerate Gamma(0, .)
  real shape_l   = fmax((1 - psi) * alpha, 1e-4);
  real shape_eps = fmax(psi * alpha, 1e-4);

  // Flat prior on psi (implicit uniform on (0,1))

  // Priors on latent variables
  l   ~ gamma(shape_l, beta_gi);
  d0  ~ gamma(a_obs, b_obs);
  eps ~ gamma(shape_eps, beta_gi);

  // Likelihood: d is the derived sibling detection delay
  for (n in 1:N) {
    if (d[n] <= 0)
      reject("derived d is non-positive: d = ", d[n]);
    target += gamma_lpdf(d[n] | a_obs, b_obs);
  }
}
