#' Convolution brm Wrappers
#'
#' @param formula A model formula
#' @param data A data.frame that must contain the date, location (as loc), primary
#' (the data that the outcome is a convolution of) and at least the outcome as 
#' specifed in `formula`.
#' @param conv_mean Numeric vector defining the prior for the log mean of the convolution 
#' (mean and standard deviation).
#' @param conv_sd Numeric vector defining the prior for the standard deviation of the
#' convolution (mean and standard deviation).
#' @param conv_max Integer defining the maximum index to use for the convolution.
#' @param conv_varying A character string defining the variability of the convolution 
#' distribution. By default this is fixed by other supported options include allowing 
#' the convolution distribution to vary across location (using a random effect).
#' @param hold_out_time Integer, number of days to hold out from the likelihood. This 
#' is useful as initially the outcome will depend on primary data outside of the range of
#' the training set. 
#' @param dry_run Logical, defaults to TRUE. For testing purposes should just the `stan`
#' code be output with not fitting done.
#' @param ... Additional parameters passed to `brms::brm`.
#'
#' @return A ""brmsfit" object.
#' @export
#' @importFrom data.table as.data.table setorder
#' @importFrom brms brm make_stancode custom_family stanvar
brm_convolution <- function(formula, data, conv_mean = c(2.5, 1),
                            conv_sd = c(1, 0.5), conv_max = 30, 
                            conv_varying = "fixed", hold_out_time = 28, 
                            dry_run = FALSE, ...) {
  
  conv_varying <- match.arg(conv_varying, choices = c("fixed", "loc"))
  
  # order data
  data <- as.data.table(data)
  data <- setorder(data, loc, date)
  data <- data[, index := 1:.N, by = loc]
               
  # get primary cases
  primary <- data$primary
               
  # filter out held out time
  data <- data[index > hold_out_time]
  
  # define custom stan code
  make_convolution_stan <- function(data, primary, max_conv, conv_varying, ut) {
    
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
    ult <- cumsum(locs_t  + ut)
    lt <- cumsum(locs_t)
    uli <- 1
    li <- 1
    if (locs > 1) {
      uli <- c(uli, 1 + ult[1:(locs - 1)])
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

// convolve a pdf and case vector, return only observed data
vector convolve(vector cases, vector rev_pmf, int ut) {
    int t = num_elements(cases);
    vector[t - ut] obs_cases;
    int max_pmf = num_elements(rev_pmf);
    vector[t] conv_cases = rep_vector(1e-5, t);
    for (s in 1:t) {
        conv_cases[s] += dot_product(cases[max(1, (s - max_pmf + 1)):s],
                                     tail(rev_pmf, min(max_pmf, s)));
    }
    obs_cases = conv_cases[(ut + 1):t];
    return(obs_cases);
  }

vector calc_pmf(real conv_mean, real conv_sd, int conv_max) {
  vector[conv_max] pmf = rep_vector(1e-5, conv_max);
      int conv_indexes[conv_max];
      for (i in 1:conv_max) {
        conv_indexes[i] = conv_max - i;
      }
      pmf = pmf + discretised_lognormal_pmf(conv_indexes, conv_mean, conv_sd, conv_max);
      return(pmf);
}"
  
  stan_functions <- c(
    stanvar(block = "functions", scode = conv_nb_lik),
    stanvar(block = "functions", scode = epinow2_funcs)
  )      
  
  stan_data <- c(
    stanvar(block = "data",
            scode = "  int locs;",
            x = locs,
            name = "locs"),
    stanvar(block = "data",
            scode = "  int ut;",
            x = ut,
            name = "ut"),
    stanvar(block = "data",
            scode = "  int primary[N + ut * locs];",
            x = primary,
            name = "primary"),
    stanvar(block = "data",
            scode = "  int uli[locs];",
            x = as.array(uli),
            name = "uli"),
    stanvar(block = "data",
            scode = "  int ult[locs];",
            x =  as.array(ult),
            name = "ult"),
    stanvar(block = "data",
            scode = "  int li[locs];",
            x = as.array(li),
            name = "li"),
    stanvar(block = "data",
            scode = "  int lt[locs];",
            x =  as.array(lt),
            name = "lt"),
    stanvar(block = "data",
            scode = "  int conv_max;",
            x =  conv_max,
            name = "conv_max")
  )
  
  if (conv_varying %in% "loc") {
    stan_parameters <- c(
      stanvar(block = "parameters",
              scode = "  
  real conv_mean;
  real<lower=0> conv_sd;
  real<lower=0> conv_mean_loc_sd;
  real<lower=0> conv_sd_loc_sd;
  real conv_mean_loc[locs];
  real<lower=0> conv_sd_loc[locs];"))
    
    stan_cmodel <- c(
      stanvar(block = "model",
              scode = paste0("  
  vector[N] conv_primary;
  
  conv_mean ~ normal(", conv_mean[1], ",", conv_mean[2], ");
  conv_mean_loc_sd ~ normal(0, 0.1) T[0,];
  conv_mean_loc ~ normal(conv_mean, conv_mean_loc_sd);
  conv_sd ~ normal(", conv_sd[1], ",", conv_sd[2], ") T[0,];
  conv_sd_loc_sd ~ normal(0, 0.1) T[0,];
  for (s in 1:locs) {
    conv_sd_loc[s] ~ normal(conv_sd, conv_mean_loc_sd) T[0,];
  }
  
  for (s in 1:locs) {
    vector[conv_max] pmf = calc_pmf(conv_mean_loc[s], conv_sd_loc[s], conv_max);
    conv_primary[li[s]:lt[s]] = convolve(to_vector(primary[uli[s]:ult[s]]), pmf, ut);
    }
      "))
    )
  }else{
    stan_parameters <- c(
      stanvar(block = "parameters",
              scode = "  
  real conv_mean;
  real<lower=0> conv_sd;"))
    
    stan_cmodel <- c(
      stanvar(block = "model",
              scode = paste0("  
  vector[N] conv_primary;
  vector[conv_max] pmf;
   
  conv_mean ~ normal(", conv_mean[1], ",", conv_mean[2], ");
  conv_sd ~ normal(", conv_sd[1], ",", conv_sd[2], ") T[0,];
  
  pmf = calc_pmf(conv_mean, conv_sd, conv_max);
  for (s in 1:locs) {
    conv_primary[li[s]:lt[s]] = convolve(to_vector(primary[uli[s]:ult[s]]), pmf, ut);
  }"))
    )
  }


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
  
conv_stan <- make_convolution_stan(data, primary,
                                   max_conv = max_conv, 
                                   conv_varying = conv_varying,
                                   ut = hold_out_time)

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