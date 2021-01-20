// all functions from EpiNow2 (epiforecasts.io/EpiNow2)
// discretised truncated lognormal pmf
vector discretised_lognormal_pmf(int[] y, real mu, real sigma, int max_val) {
  int n = num_elements(y);
  vector[n] pmf;
  real small = 1e-5;
  real c_sigma = sigma < small ? small : sigma;
  real c_mu = mu < small ? small : mu;
  vector[n] adj_y = to_vector(y) + small;
  vector[n] upper_y = (log(adj_y + 1) - c_mu) / c_sigma;
  vector[n] lower_y = (log(adj_y) - c_mu) / c_sigma;
  real max_cdf = normal_cdf((log(max_val + small) - c_mu) / c_sigma, 0.0, 1.0);
  real min_cdf = normal_cdf((log(small) - c_mu) / c_sigma, 0.0, 1.0);
  real trunc_cdf = max_cdf - min_cdf;
  for (i in 1:n) {
    pmf[i] = (normal_cdf(upper_y[i], 0.0, 1.0) - normal_cdf(lower_y[i], 0.0, 1.0)) /
    trunc_cdf;
  }
  return(pmf);
}

// convolve a pdf and case vector
vector convolve(vector cases, vector rev_pmf) {
    int t = num_elements(cases);
    int max_pmf = num_elements(rev_pmf);
    vector[t] conv_cases = rep_vector(1e-5, t);
    for (s in 1:t) {
        conv_cases[s] += dot_product(cases[max(1, (s - max_pmf + 1)):s],
                                     tail(rev_pmf, min(max_pmf, s)));
    }
   return(conv_cases);
  }


// convolve latent infections to reported (but still unobserved) cases
vector convolve_to_report(vector infections,
                          real[] delay_mean,
                          real[] delay_sd,
                          int[] max_delay,
                          int seeding_time) {
  int t = num_elements(infections);
  vector[t - seeding_time] reports;
  vector[t] unobs_reports = infections;
  int delays = num_elements(delay_mean);
  if (delays) {
    for (s in 1:delays) {
      vector[max_delay[s]] pmf = rep_vector(1e-5, max_delay[s]);
      int delay_indexes[max_delay[s]];
      for (i in 1:max_delay[s]) {
        delay_indexes[i] = max_delay[s] - i;
      }
      pmf = pmf + discretised_lognormal_pmf(delay_indexes, delay_mean[s],
                                            delay_sd[s], max_delay[s]);
      unobs_reports = convolve(unobs_reports, pmf);
    }
    reports = unobs_reports[(seeding_time + 1):t];
  }else{
    reports = infections[(seeding_time + 1):t];
  }
  return(reports);
}
