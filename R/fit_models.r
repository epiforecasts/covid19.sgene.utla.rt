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
  mutate(tier = if_else(tier == "none", "_none", tier)) %>% 
  mutate(sgtf_samples = as.integer(prop_sgtf * samples),
         id = 1:n(), sgtf = NA_real_)

# model with uncertainty ------------------------------------------------

source(here("R", "variant_rt.r"))

# no uncertainty on sgtf
fit_cert <- variant_rt(
  log_rt = ~ 1,
  data = utla_rt_with_covariates %>%
    rename(rt_mean = rt_mean_long_gt, rt_sd = rt_sd_long_gt)%>%
    mutate(sgtf = prop_sgtf) %>%
    mutate(
      sgtf = ifelse(sgtf <= 0, 1e-5, ifelse(prop_sgtf >= 1, 0.9999, prop_sgtf))
    ),
  brm_fn = brm
)

summary(fit_cert)

# uncertainty on sgtf
fit_ <- variant_rt(
  log_rt = ~ 1,
  data = utla_rt_with_covariates %>%
    rename(rt_mean = rt_mean_long_gt, rt_sd = rt_sd_long_gt),
  brm_fn = brm
)

summary(fit)


# Define custom family ------------------------------------------------
add_var_student <- custom_family(
  "add_var_student", dpars = c("mu", "sigma", "nu", "alpha"),
  links = c("log", "identity", "identity", "identity"),
  lb = c(NA, 0, 1, 0),
  type = "real",
  vars = "vreal1[n]"
)

stan_funs <- "
real add_var_student_lpdf(real y, real mu, real sigma, real nu, real alpha,
                          real f) {
    real combined_mu = (1 + (alpha - 1) * f) * mu;
    return student_t_lpdf(y | nu, combined_mu, sigma);
                            }
real add_var_student_rng(real mu, real sigma, real nu, real alpha, real f) {
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
fit_models <- function(gt, data, main_only = TRUE, parallel = TRUE) {
  # filter for target
  dynamic_data <- data %>%
    rename_with(~ sub(paste0("_", gt, "_gt"), "", .x)) %>%
    filter(!is.na(rt_mean))
  static_data <- dynamic_data %>%
    filter(week_infection == max(week_infection))

  ##Static model ------------------------------------------------------------
  # set model settings and priors
  static_model <- function(form, ...) {
    base_model(form = form, data = static_data,
               control = list(adapt_delta = 0.95),
               ...)
  }
  # fit models
  static <- list()
  if (!main_only) {
    static[["intercept"]] <- static_model(rt_mean | vreal(prop_sgtf) ~ 1,
                                          prior = priors)
    static[["region"]] <-
      static_model(rt_mean | vreal(prop_sgtf) ~ nhser_name,
                   prior = c(priors,
                             prior(student_t(3, 0, 0.5), class = "b")))
  }

  # Dynamic model -----------------------------------------------------------
  # set model settings and priors
  dynamic_model <- function(form, iter = 2000, ...) {
    base_model(form = form, data = dynamic_data,
               prior = c(priors,
                         prior(student_t(3, 0, 0.5), class = "b")),
               control = list(adapt_delta = 0.95, max_treedepth = 12),
               iter = iter, ...)
  }
  # fit models

  dynamic_models <- list()
  if (!main_only) {
    dynamic_models[["interventions_only"]] <-
      as.formula(rt_mean | vreal(prop_sgtf) ~ tier)

    dynamic_models[["interventions"]] <-
      as.formula(rt_mean | vreal(prop_sgtf) ~ tier +
                      retail_and_recreation + grocery_and_pharmacy +
                      parks + transit_stations + workplaces + residential)

    dynamic_models[["interventions_random"]] <-
      as.formula(rt_mean | vreal(prop_sgtf) ~ tier + (1 | utla_name) +
                      retail_and_recreation + grocery_and_pharmacy +
                      parks + transit_stations + workplaces + residential)

    dynamic_models[["interventions_region"]] <-
      as.formula(rt_mean | vreal(prop_sgtf) ~ tier + nhser_name +
                      retail_and_recreation + grocery_and_pharmacy +
                      parks + transit_stations + workplaces + residential)
  }

  dynamic_models[["interventions_random_region"]] <-
    as.formula(rt_mean | vreal(prop_sgtf) ~ tier +  (1 | utla_name) +
                 nhser_name + retail_and_recreation + grocery_and_pharmacy +
                 parks + transit_stations + workplaces + residential)

  if (!main_only) {
  dynamic_models[["interventions_time_region"]] <-
    as.formula(rt_mean | vreal(prop_sgtf) ~ tier + s(time, k = 9) + nhser_name +
                    retail_and_recreation + grocery_and_pharmacy +
                    parks + transit_stations + workplaces + residential)
  }

  dynamic_models[["interventions_time_region_random"]] <-
    as.formula(rt_mean | vreal(prop_sgtf) ~ tier + s(time, k = 9) +
                    (1 | utla_name) + nhser_name +
                    retail_and_recreation + grocery_and_pharmacy +
                    parks + transit_stations + workplaces + residential)

  if (!main_only) {
    dynamic_models[["interventions_time_by_region"]] <-
      as.formula(rt_mean | vreal(prop_sgtf) ~ tier +
                   s(time, k = 9, by = nhser_name) +
                   retail_and_recreation + grocery_and_pharmacy +
                   parks + transit_stations + workplaces + residential)
  }

  dynamic_models[["interventions_time_by_random_region"]] <-
    as.formula(rt_mean | vreal(prop_sgtf) ~ tier +
                 s(time, k = 9, by = nhser_name) +
                 (1 | utla_name) +
                 retail_and_recreation + grocery_and_pharmacy +
                 parks + transit_stations + workplaces + residential)

  if (!main_only) {
    dynamic_models[["interventions_independent_time_region"]] <-
      as.formula(rt_mean | vreal(prop_sgtf) ~ tier + factor(time):nhser_name +
                      retail_and_recreation + grocery_and_pharmacy +
                      parks + transit_stations + workplaces + residential)

    dynamic_models[["interventions_independent_time_random_region"]] <-
      as.formula(rt_mean | vreal(prop_sgtf) ~ tier + factor(time):nhser_name +
                      (1 | utla_name) +
                      retail_and_recreation + grocery_and_pharmacy +
                      parks + transit_stations + workplaces + residential)
  }

  if (parallel) {
    dynamic <- parallel::mclapply(dynamic_models, dynamic_model)
  } else {
    dynamic <- lapply(dynamic_models, dynamic_model)
  }

  return(list(models = list(static = static, dynamic = dynamic),
              data = list(static = static_data, dynamic = dynamic_data)))
}

# fit models
gt <- c("short", "long")
res <- lapply(gt, fit_models,
              data = utla_rt_with_covariates)
names(res) <- gt

# Save fits ---------------------------------------------------------------
output_path <- here("output")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

save_results <- function(name) {
  saveRDS(res[[name]],
          here::here("output", paste0("sgene_fits_", name, "_gt.rds")))
}
lapply(gt, save_results)
