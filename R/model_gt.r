# Packages ----------------------------------------------------------------
library(here)
library(dplyr)
library(brms)
library(parallel)

# Options -----------------------------------------------------------------
options(mc.cores = detectCores())

# Get data ----------------------------------------------------------------
utla_rt_with_covariates <-
  readRDS(here("data", "utla_rt_with_covariates.rds")) %>%
  filter(week_infection > "2020-10-01") %>%
  mutate(tier = if_else(tier == "none", "_none", tier))



# Add custom functions ----------------------------------------------------
form <- bf(rt_mean ~ (1 + ((1 + afk) * kp) * ((1 + afr) * rp) * ((1 + afG) * Gp))^(1/((1 + afk) * kp)),
           rp ~ 1,
           kp ~ 1,
           Gp ~ 1,
           afk ~ 1,
           afG ~ 1,
           afr ~ 1, 
           nl = TRUE)

fit <- make_stancode(formula = form, family = gaussian(), 
           data = utla_rt_with_covariates)

# Add custom family -------------------------------------------------------
stan_funs <- "
real var_gt_lpdf(real y, real mu, real sigma, real nu, real alpha,
                          real f) {
                          
    real R = (1 - f) * pow((1 + k * r * G), 1 / k);
    R = R + (1 - f) * pow((1 + k * r * G), 1 / k)
    real combined_mu = (1 + (alpha - 1) * f) * mu;
    return student_t_lpdf(y | nu, combined_mu, sigma); 
                            }
real var_gt_rng(real mu, real sigma, real nu, real alpha, real f) {
    real combined_mu = (1 + (alpha - 1) * f) * mu;
    return student_t_rng(nu, combined_mu, sigma);
  }
"
stanvars <- stanvar(block = "functions", scode = stan_funs)

# Set up shared priors ----------------------------------------------------
priors <- c(prior(gamma(2, 0.1), class = nu),
            prior(lognormal(0, 1), class = alpha),
            prior(student_t(3, 0, 0.5), class = sigma))

# Set up model ------------------------------------------------------------
base_model <- function(form, iter = 2000, ...) {
  brm(formula = form,
      family = add_var_student,
      stanvars = stanvars,
      warmup = 500, iter = iter, ...)
}

# Fit models --------------------------------------------------------------
# filter for target
dynamic_data <- data %>%
  rename_with(~ sub(paste0("_", gt, "_gt"), "", .x)) %>%
  filter(!is.na(rt_mean))

# Dynamic model -----------------------------------------------------------
# set model settings and priors
dynamic_model <- function(form, iter = 2000, ...) {
  base_model(form = form, data = dynamic_data,
             prior = c(priors,
                       prior(student_t(3, 0, 0.5), class = "b")),
             control = list(adapt_delta = 0.95, max_treedepth = 12),
             iter = iter, ...)
}
parks + transit_stations + workplaces + residential)
}

dynamic_models[["interventions_random_region"]] <-
  as.formula(rt_mean | vreal(prop_sgtf) ~ tier +  (1 | utla_name) +
               nhser_name + retail_and_recreation + grocery_and_pharmacy +
               parks + transit_stations + workplaces + residential)




