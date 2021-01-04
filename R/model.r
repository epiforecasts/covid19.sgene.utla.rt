# Packages ----------------------------------------------------------------
library(here)
library(dplyr)
library(brms)

# Options -----------------------------------------------------------------
options(mc.cores = 4)

# Get data ----------------------------------------------------------------
ltla_rt_with_covariates <- readRDS(here("data", "ltla_rt_with_covariates.rds"))

# filter for type
ltla_rt_with_covariates <- ltla_rt_with_covariates %>% 
  filter(generation_time %in% "short")

# add small amount of noise to 0 measurement error
ltla_rt_with_covariates <- ltla_rt_with_covariates %>%
  mutate(prop_variant_sd = ifelse(prop_variant_sd == 0, 1e-5, prop_variant_sd))

static_rt <- ltla_rt_with_covariates %>%
  filter(week_infection == max(week_infection))

# Add custom family -------------------------------------------------------
add_var_student <- custom_family(
  "add_var_student", dpars = c("mu", "sigma", "nu", "alpha"),
  links = c("log", "identity", "identity", "identity"),
  lb = c(NA, 0, 1, NA),
  type = "real",
  vars = "vreal1[n]"
)

stan_funs <- "
real add_var_student_lpdf(real y, real mu, real sigma, real nu, real alpha,
                          real f) {
    real combined_mu = (1 + alpha * f) * mu;
    return student_t_lpdf(y | nu, combined_mu, sigma);
                            }
real add_var_student_rng(real mu, real sigma, real nu, real alpha, real f) {
    real combined_mu = (1 + alpha * f) * mu;
    return student_t_rng(nu, combined_mu, sigma);
  }
"
stanvars <- stanvar(block = "functions", scode = stan_funs)

# Set up shared priors ----------------------------------------------------
priors <- c(prior(gamma(2, .1), class = nu),
            prior(student_t(3, 0, 0.5), class = alpha),
            prior(student_t(3, 0, 0.5), class = sigma))

# Set up model ------------------------------------------------------------
base_model <- function(form, iter = 2000, ...) {
  brm(formula = form,
      family = add_var_student,
      stanvars = stanvars, 
      warmup = 500, iter = iter, ...)
}
## Static model ------------------------------------------------------------
# set model settings and priors
static_model <- function(form, ...) {
   base_model(form = form, data = static_rt, 
              control = list(adapt_delta = 0.95), ...)
}

static <- list()
static[["intercept"]] <- static_model(rt_mean | vreal(prop_variant) ~ 1,
                                      prior = priors)
static[["region"]] <- static_model(rt_mean | vreal(prop_variant) ~ nhser_name,
                                   prior = c(priors,
                                             prior(student_t(3, 0, 0.5), class = "b")))

# Dynamic model -----------------------------------------------------------
# set model settings and priors
dynamic_model <- function(form, iter = 2000, ...) {
  base_model(form = form, data = ltla_rt_interventions,
             prior = c(priors,
                       prior(student_t(3, 0, 0.5), class = "b")),
             control = list(adapt_delta = 0.95, max_treedepth = 12),
             iter = iter, ...)
}

# fit models
dynamic <- list()
dynamic[["interventions_only"]] <- 
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier)

dynamic[["interventions"]] <- 
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_random"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + (1 | ltla_name) + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + nhser_name + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_random_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier +  (1 | ltla_name) + nhser_name + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_time_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + s(time, k = 10) + nhser_name + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_time_by_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + s(time, k = 10, by = nhser_name) + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_time_by_random_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + s(time, k = 10, by = nhser_name) + 
                 (1 | ltla_name) + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_independent_time_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + factor(time):nhser_name + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_independent_time_random_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + factor(time):nhser_name +
                  (1 | ltla_name) + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

# Save fits ---------------------------------------------------------------
output_path <- here("output")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

saveRDS(list(models = list(static = static, dynamic = dynamic),
             data = list(static = static_rt, dynamic = ltla_rt_interventions)),
        file.path(output_path, "sgene_fits.rds"))
