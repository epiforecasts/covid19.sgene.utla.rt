## Define model ------------------------------------------------------------
## D = c^{+}(1-f)C + c^{-}fC
## where c is the CFR, f is the fraction of cases that are sgtf, C is the number
## of cases by date of infection, and D is the number of deaths by date of
## infection. Effect of the variant on the CFR is assumed to be:
## c^{-} = \alpha c^{+} or c^{-} = \alpha + c^{+}
## leading to:
## D = (1 - (\alpha - 1)f)c^{+}C or D = (c^{+} + \alpha f)C

## define custom negative binomial family including variant factor and cases
variant_nb <- function(additive = FALSE) {
  custom_family(
    "variant_nb", dpars = c("mu", "phi", "alpha", "epsilon"),
    links = c("logit", "log", "identity", "log"),
    lb = c(0, 0, ifelse(!additive, 0, NA), 0),
    type = "int",
    vars = c("f[n]", "cases[n]", "effect[1]")
  )
}

## define stan code to scale cfr by cases and variant fraction
make_stanvars <- function(data, additive = FALSE) {
  stan_funs <- "
real variant_nb_lpmf(int y, real mu, real phi, real alpha, real epsilon,
                     real f, int cases, int effect) {
    real scaled_cases;
    if (effect) {
      scaled_cases = (mu + alpha * f) * cases + epsilon;
    }else {
      scaled_cases = (1 + (alpha - 1) * f) * mu * cases + epsilon;
    }
    return  neg_binomial_2_lpmf(y | scaled_cases, phi);
                            }
real variant_nb_rng(int y, real mu, real phi, real alpha, real epsilon,
                    real f, int cases, int effect) {
    real scaled_cases;
    if (effect) {
      scaled_cases = (mu + alpha * f) * cases + epsilon;
    }else {
      scaled_cases = (1 + (alpha - 1) * f) * mu * cases + epsilon;
    }
    return  neg_binomial_2_rng(scaled_cases, phi);
                            }
"
  stanvars <- c(stanvar(block = "functions", scode = stan_funs),
                stanvar(block = "data",
                        scode = "  real f[N];",
                        x = data$prop_sgtf,
                        name = "f"),
                stanvar(block = "data",
                        scode = "  int cases[N];",
                        x = data$cases,
                        name = "cases"),
                stanvar(block = "data",
                        scode = "  int effect[1];",
                        x = array(as.numeric(additive)),
                        name = "effect")
  )
  return(stanvars)
}

## define model function
variant_model <- function(form, iter = 2000, data = deaths_with_cov,
                     additive = FALSE, ...) {
  # define priors
  if (additive) {
    priors <- c(prior(normal(0, 0.01), class = alpha))
  }else{
    priors <- c(prior(lognormal(0, 0.5), class = alpha))
  }
  message("Fitting: ", as.character(form))
  brm(formula = form,
      family = variant_nb(additive = additive),
      prior = priors,
      data,
      stanvars = make_stanvars(data, additive = additive),
      control = list(adapt_delta = 0.99),
      warmup = 1000, iter = iter, ...)
}


log_lik_variant_nb <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  alpha <- prep$dpars$alpha
  epsilon <- prep$dpars$epsilon
  f <- prep$data$f[i]
  y <- prep$data$Y[i]
  cases <- prep$data$cases[i]
  effect <- prep$data$effect[1]
  variant_nb_lpmf(y, mu, phi, alpha, epsilon, f, cases, effect)
}

posterior_predict_variant_nb <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  alpha <- prep$dpars$alpha
  f <- prep$data$f[i]
  y <- prep$data$Y[i]
  cases <- prep$data$cases[i]
  effect <- prep$data$effect[1]
  variant_nb_rng(mu, phi, alpha, epsilon, f, cases, effect)
}