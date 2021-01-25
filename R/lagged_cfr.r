# Packages ----------------------------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(vroom)
library(lubridate)
library(brms)
library(parallel)
library(loo)

# set number of parallel cores
no_cores <- detectCores()

# Get functions -----------------------------------------------------------
source(here("R", "variant_model.r"))
source(here("R", "get_infections_data.r"))

# Get data ----------------------------------------------------------------
df <- list()
df[["cfr"]] <- get_infections_data("cases", "deaths", "lagged")
df[["chr"]] <- get_infections_data("cases", "admissions", "lagged")
admissions_lag <- unique(df[["chr"]]$lag)
df[["hfr"]] <- get_infections_data("admissions", "deaths", "lagged",
                        infections_lag = 7 + admissions_lag)

## Fit models --------------------------------------------------------------
models <- list()
## define models to fit
models[["intercept"]] <- as.formula(deaths ~ 1)
models[["time"]] <- as.formula(deaths ~ (1 | time))
models[["utla"]] <- as.formula(deaths ~ (1 | utla))
models[["all"]] <- as.formula(deaths ~ (1 | utla) + (1 | time))

## core usage
if (no_cores <= 4) {
  options(mc.cores = no_cores)
  mc_cores <- 1
} else {
  options(mc.cores = 4)
  mc_cores <- ceiling(no_cores / 4)
}

## fit models
results <- lapply(names(df), function(x) {
  fits <- list()
  cat(x, "multiplicative\n")
  fits[["multiplicative"]] <-
    mclapply(models, variant_model, data = df[[x]], mc.cores = mc_cores)
  cat(x, "additive\n")
  fits[["additive"]] <-
    mclapply(models, variant_model, data = df[[x]], mc.cores = mc_cores,
           additive = TRUE)
  return(fits)
})
names(results) <- names(df)

## variant effect ----------------------------------------------------------
extract_variant_effect <- function(x, additive = FALSE) {
  samples <- posterior_samples(x, "alpha")
  q <- samples[, "alpha"]
  q <- quantile(q, c(0.025, 0.5, 0.975))
  q <- signif(q, 2)
  return(q)
}

var_res <- lapply(names(results), function(x) {
  variant_effect <- list()
  variant_effect[["multiplicative"]] <-
    lapply(results[[x]][["multiplicative"]], extract_variant_effect)
  variant_effect[["additive"]] <-
    lapply(results[[x]][["additive"]], extract_variant_effect, additive = TRUE)
  return(variant_effect)
})
names(var_res) <- names(df)

## Compare models ----------------------------------------------------------
## requires custom log_lik functions to be implemented for variant_nb family
expose_functions(results[[1]][[1]][[1]], vectorize = TRUE)

add_loo <- function(fits) {
  fits %>%
    lapply(add_criterion, "loo") %>%
    lapply(loo, save_psis = TRUE)
}

options(mc.cores = no_cores)
model_loos <- lapply(names(results), function(x) {
  fits <- results[[x]]
  loos <- list()
  loos[["multiplicative"]] <- add_loo(fits[["multiplicative"]])
  loos[["additive"]] <- add_loo(fits[["additive"]])
  lc <- list()
  lc[["multiplicative"]] <- loo_compare(loos[["multiplicative"]])
  lc[["additive"]] <- loo_compare(loos[["additive"]])
  all_loos <- flatten(loos)
  names(all_loos) <- c(paste0(names(all_loos)[1:length(models)], "_multiplictive"),
                       paste0(names(all_loos)[1:length(models)], "_additive"))
  lc[["all"]] <- loo_compare(all_loos)
  return(list(loos = loos, lc = lc))
})
names(model_loos) <- names(df)

# Save results ------------------------------------------------------------
to_save <- lapply(names(df), function(x) {
  output <- list()
  output$data <- df[[x]]
  output$fits <- results[[x]]
  output$effect <- var_res[[x]]
  output$loos <- model_loos[[x]]
  return(output)
})
names(to_save) <- names(df)

saveRDS(to_save, here("output", "associations.rds"))
