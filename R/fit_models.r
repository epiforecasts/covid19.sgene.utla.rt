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
  mutate(tier = if_else(tier == tier[1], paste0("_", tier), tier)) %>%
  mutate(prop = 1 - prop) %>%
  mutate(tier = if_else(tier == "none", "_none", tier)) %>% 
  mutate(positive_samples = as.integer(prop * samples),
         id = 1:n(), positive_prop = NA_real_)

# model with uncertainty ------------------------------------------------

source(here("R", "variant_rt.r"))

# no uncertainty on sgtf
fit_cert <- variant_rt(
  log_rt = ~ 1,
  data = utla_rt_with_covariates %>%
    rename(rt_mean = rt_mean_long_gt, rt_sd = rt_sd_long_gt)%>%
    mutate(positive_prop = prop) %>%
    mutate(
      positive_prop = ifelse(positive_prop <= 0, 1e-5, ifelse(prop >= 1, 0.9999, prop))
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
fit_models <- function(gt, data, main_only = TRUE, parallel = TRUE,
                       type = c("static", "dynamic"),
                       static_weeks = 1) {
  # filter for target
  dynamic_data <- data %>%
     rename_with(~ sub(paste0("_", gt, "_gt"), "", .x)) %>%
     filter(!is.na(rt_mean))
  static_data <-
    lapply(seq_len(static_weeks),
           function(x) {
             dynamic_data %>%
               filter(week_infection == (max(week_infection) - (x - 1) * 7))
           })

  ##Static model ------------------------------------------------------------
  # set model settings and priors
  static_model <- function(form, i, ...) {
    base_model(form = form, data = static_data[[i]],
               control = list(adapt_delta = 0.95),
               ...)
  }
  # fit models
  static <- list()
  if ("static" %in% type) {
    static[["intercept"]] <-
      lapply(seq_len(static_weeks), function(i) {
        static_model(rt_mean | vreal(prop) ~ 1, i, prior = priors)
      })
    static[["region"]] <-
      lapply(seq_len(static_weeks), function(i) {
        static_model(rt_mean | vreal(prop) ~ nhser_name, i,
                     prior = c(priors,
                               prior(student_t(3, 0, 0.5), class = "b")))
      })
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
  if ("dynamic" %in% type) {
    if (!main_only) {
      dynamic_models[["interventions_only"]] <-
        as.formula(rt_mean | vreal(prop) ~ tier)

      dynamic_models[["interventions"]] <-
        as.formula(rt_mean | vreal(prop) ~ tier +
          retail_and_recreation + transit_stations + workplaces + residential)

      dynamic_models[["interventions_random"]] <-
        as.formula(rt_mean | vreal(prop) ~ tier + (1 | utla_name) +
          retail_and_recreation + transit_stations + workplaces + residential)

      dynamic_models[["interventions_region"]] <-
        as.formula(rt_mean | vreal(prop) ~ tier + nhser_name +
          retail_and_recreation + transit_stations + workplaces + residential)
    }

    dynamic_models[["interventions_random_region"]] <-
      as.formula(rt_mean | vreal(prop) ~ tier + (1 | utla_name) +
        nhser_name + retail_and_recreation + transit_stations + workplaces + residential)

    if (!main_only) {
      dynamic_models[["interventions_time_region"]] <-
        as.formula(rt_mean | vreal(prop) ~ tier + s(time, k = 9) +
                     nhser_name + retail_and_recreation + transit_stations + workplaces + residential)
    }

    dynamic_models[["interventions_time_region_random"]] <-
      as.formula(rt_mean | vreal(prop) ~ tier + s(time, k = 9) +
        (1 | utla_name) + nhser_name +
        retail_and_recreation + transit_stations + workplaces + residential)

    if (!main_only) {
      dynamic_models[["interventions_time_by_region"]] <-
        as.formula(rt_mean | vreal(prop) ~ tier +
          s(time, k = 9, by = nhser_name) +
          retail_and_recreation + transit_stations + workplaces + residential)
    }

    dynamic_models[["interventions_time_by_random_region"]] <-
      as.formula(rt_mean | vreal(prop) ~ tier +
        s(time, k = 9, by = nhser_name) +
        (1 | utla_name) +
        retail_and_recreation + transit_stations + workplaces + residential)

    if (!main_only) {
      dynamic_models[["interventions_independent_time_region"]] <-
        as.formula(rt_mean | vreal(prop) ~ tier + factor(time):nhser_name +
          retail_and_recreation + transit_stations + workplaces + residential)

      dynamic_models[["interventions_independent_time_random_region"]] <-
        as.formula(rt_mean | vreal(prop) ~ tier + factor(time):nhser_name +
          (1 | utla_name) +
          retail_and_recreation + transit_stations + workplaces + residential)
    }

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
              data = utla_rt_with_covariates,
              type = c("static", "dynamic"),
              static_weeks = 6)
names(res) <- gt

# Save fits ---------------------------------------------------------------
output_path <- here("output")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

save_results <- function(name) {
  saveRDS(res[[name]],
          here::here("output", paste0("sgene_fits_", name, "_gt.rds")))
}
lapply(gt, save_results)

