#include <Rcpp.h>
#include <queue>
#include <vector>
#include <functional>

// Use R's sample.int for RNG compatibility with the pure-R implementation
static Rcpp::Function r_sample_int("sample.int");

// [[Rcpp::export]]
Rcpp::NumericVector sim_stochastic_rcpp(int n, Rcpp::Function gen_inf_attempts,
                                        int maxinf = -1) {

  // If maxinf not specified (or negative), default to n (full population)
  if (maxinf < 0) maxinf = n;

  // Initialize infection times to Inf
  Rcpp::NumericVector tinf_vec(n, R_PosInf);

  // Min-heap priority queue for event times
  std::priority_queue<double, std::vector<double>, std::greater<double>> pq;

  // Seed infection: pick random index case using R's sample.int (1-based)
  int indexcase = Rcpp::as<int>(r_sample_int(n, 1)) - 1; // 0-based
  tinf_vec[indexcase] = 0.0;
  int n_infected = 1;

  // Generate initial infection attempts from index case
  Rcpp::NumericVector init_attempts = gen_inf_attempts(0.0);
  for (int i = 0; i < init_attempts.size(); i++) {
    pq.push(init_attempts[i]);
  }

  // Process events in chronological order (min-heap gives smallest first)
  while (!pq.empty()) {
    double t_attempt = pq.top();
    pq.pop();

    // Pick a random target using R's sample.int (1-based -> 0-based)
    int target = Rcpp::as<int>(r_sample_int(n, 1)) - 1;

    if (tinf_vec[target] == R_PosInf) {
      // Successful infection
      tinf_vec[target] = t_attempt;
      n_infected++;

      if (n_infected >= maxinf) break;

      // Generate new infection attempts from this newly infected person
      Rcpp::NumericVector new_attempts = gen_inf_attempts(t_attempt);
      for (int i = 0; i < new_attempts.size(); i++) {
        pq.push(new_attempts[i]);
      }
    }
  }

  return tinf_vec;
}

// [[Rcpp::export]]
Rcpp::NumericVector sim_infinite_pop_rcpp(int max_cases, double tmax,
                                          Rcpp::Function gen_inf_attempts) {

  std::vector<double> tinf_vec;
  tinf_vec.reserve(max_cases);

  // Min-heap priority queue for event times
  std::priority_queue<double, std::vector<double>, std::greater<double>> pq;

  // Seed: index case at t=0
  tinf_vec.push_back(0.0);

  // Generate initial infection attempts from index case
  Rcpp::NumericVector init_attempts = gen_inf_attempts(0.0);
  for (int i = 0; i < init_attempts.size(); i++) {
    pq.push(init_attempts[i]);
  }

  // Process events in chronological order
  while (!pq.empty() && (int)tinf_vec.size() < max_cases) {
    double t_attempt = pq.top();
    pq.pop();

    if (t_attempt > tmax) break;

    // Every attempt succeeds (infinite susceptible population)
    tinf_vec.push_back(t_attempt);

    // Only generate new attempts if we still need more cases
    if ((int)tinf_vec.size() < max_cases) {
      Rcpp::NumericVector new_attempts = gen_inf_attempts(t_attempt);
      for (int i = 0; i < new_attempts.size(); i++) {
        pq.push(new_attempts[i]);
      }
    }
  }

  return Rcpp::wrap(tinf_vec);
}
