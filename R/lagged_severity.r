# Packages ----------------------------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(vroom)
library(lubridate)
library(brms)
library(future)
library(future.apply)
library(loo)

# set number of parallel cores
no_cores <- availableCores()

# Get functions -----------------------------------------------------------
source(here("R", "lag_model.r"))
source(here("R", "get_infections_data.r"))

# Get data ----------------------------------------------------------------
df <- list()

# all data combinations by UTLA
df[["utla"]][["cfr"]] <- get_infections_data("cases", "deaths", "lagged")
df[["utla"]][["chr"]] <- get_infections_data("cases", "admissions", "lagged")
df[["utla"]][["hfr"]] <- get_infections_data("admissions", "deaths", "lagged",
                                             infections_lag = 7 + unique(df[["utla"]][["chr"]]$lag))

# all data combinations by NHS region
df[["region"]][["cfr"]] <- get_infections_data("cases", "deaths", "lagged", level = "nhs region")
df[["region"]][["chr"]] <- get_infections_data("cases", "admissions", "lagged", level = "nhs region")
df[["region"]][["hfr"]] <- get_infections_data("admissions", "deaths", "lagged", level = "nhs region",
                                               infections_lag = 7 + unique(df[["region"]][["chr"]]$lag))

# Define models -----------------------------------------------------------
models <- list()
## define models to fit
models[["intercept"]] <- as.formula(deaths ~ 1)
models[["primary"]] <- as.formula(deaths ~ s(normalised_cases, k = 5))
models[["loc"]] <- as.formula(deaths ~ (1 | loc))
models[["all"]] <- as.formula(deaths ~ (1 | loc) + s(normalised_cases, k = 5))

## Fit models --------------------------------------------------------------
## core usage
if (no_cores <= 4) {
  stan_cores <- no_cores
  mc_cores <- 1
} else {
  stan_cores <- 4
  mc_cores <- ceiling(no_cores / 4)
}
plan("multisession", workers = mc_cores, earlySignal = TRUE)

# fit model grid in parallel
fit_targets <- expand_grid(loc = c("utla", "region"), 
                           effect_type = c("multiplicative", "additive"), 
                           target = c("cfr", "chr", "hfr"))

fits <- future_lapply(1:nrow(fit_targets), function(i) {
  ft <- fit_targets[i, ]
  ft$data <- list(df[[ft$loc]][[ft$target]])
  message("Fitting ", ft$target, " at the ", ft$loc, " level")
  out <- list()
  fits <- suppressMessages(lapply(models, variant_model,
                                 cores = stan_cores,
                                 data = ft$data[[1]],
                                 additive = ft$effect_type %in% "additive"))
  ft$models <- list(names(models))
  ft$fit <- list(fits)
  ft <- unnest(ft, cols = c("models", "fit"))
  return(ft)},
  future.scheduling = Inf, future.seed = TRUE)

fits <- reduce(fits, bind_rows)

## variant effect ----------------------------------------------------------
fits <- fits %>% 
  mutate(variant_effect_q = map(fit, function(x) {
    samples <- posterior_samples(x, "alpha")
    q <- samples[, "alpha"]
    q <- quantile(q, c(0.025, 0.5, 0.975))
    q <- signif(q, 2)
    return(q)
  })) %>% 
  mutate(variant_effect =  map_chr(variant_effect_q,
                                   ~ paste0(.[2]," (", .[1], ", ", .[3], ")")))

## Compare models ----------------------------------------------------------
## requires custom log_lik functions to be implemented for variant_nb family
expose_functions(fits$fit[[1]], vectorize = TRUE)

fits <- fits %>% 
  mutate(loo = map(fit, loo, save_psis = TRUE))

# Save results ------------------------------------------------------------
saveRDS(fits, here("output", "lagged_severity.rds"))
