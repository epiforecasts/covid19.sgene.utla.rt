library(data.table)
library(brms)

convolution_model <- function(formula, data, max_conv = 30, dry_run = FALSE, ...) {
  
  # order data
  data <- as.data.table(data)
  data <- setorder(data, loc, date)
    
  # define custom stan code
  make_convolution_stan <- function(data, max_conv) {
    
    # custom family
    conv_nb <- function() {
      custom_family(
        "conv_nb", dpars = c("mu", "phi"),
        links = c("logit", "log"),
        lb = c(0, 0),
        type = "int",
        vars = c("conv_primary[n]")
      )
    }
    
    # get time per location and indexing
    locs_t <- copy(data)[, .(.N), by = loc]$N
    locs <- length(unique(data$loc))
    lt <- cumsum(locs_t)
    li <- 1
    if (locs > 1) {
      li <- c(li, 1 + lt[1:(locs - 1)])
      
    }
    
    conv_nb_lik <- "
real conv_nb_lpmf(int y, real mu, real phi, real conv_primary) {
    real scaled_primary = mu * conv_primary;
    return  neg_binomial_2_lpmf(y | scaled_primary, phi);
                            }
real conv_nb_rng(int y, real mu, real phi, real conv_primary) {
    real scaled_primary = mu * conv_primary;
    return  neg_binomial_2_rng(scaled_primary, phi);
                            }
"
    
    epinow2_funcs <- "
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
}"
  
  stan_functions <- c(
    stanvar(block = "functions", scode = conv_nb_lik),
    stanvar(block = "functions", scode = epinow2_funcs)
  )      
  
  stan_data <- c(
    stanvar(block = "data",
            scode = "  int primary[N];",
            x = data$primary,
            name = "primary"),
    stanvar(block = "data",
            scode = "  int locs;",
            x = locs,
            name = "locs"),
    stanvar(block = "data",
            scode = "  int li[locs];",
            x = as.array(li),
            name = "li"),
    stanvar(block = "data",
            scode = "  int lt[locs];",
            x =  as.array(lt),
            name = "lt"),
    stanvar(block = "data",
            scode = "  int conv_max[1];",
            x =  as.array(max_conv),
            name = "conv_max")
  )
  
  stan_parameters <- c(
    stanvar(block = "parameters",
            scode = "  
            real conv_mean[1];
            real conv_sd[1];") )
  
  
  stan_cmodel <- c(
    stanvar(block = "model",
            scode = "  
  vector[N] conv_primary;
  for (s in 1:locs) {
    conv_primary[li[s]:lt[s]] = 
          convolve_to_report(to_vector(primary[li[s]:lt[s]]), conv_mean, conv_sd, conv_max, 0);
    }
      ")
  )

  stanvars <- c(
    stan_functions,
    stan_data, 
    stan_parameters,
    stan_cmodel
  )
  
  if (length(stanvars) == 0) {
    stop("Custom stan code incorrectly defined. This is likely an issue with the input data or parameters")
  }
  
  return(list(family = conv_nb, other = stanvars))
}
  
conv_stan <- make_convolution_stan(data, max_conv = max_conv)

if (dry_run) {
  brm_fn <- make_stancode
}else{
  brm_fn <- brm
}
# fit model
fit <- brm_fn(formula = formula,
              family = conv_stan$family(),
              data = data,
              stanvars = conv_stan$other,
              ...)
return(fit)
}