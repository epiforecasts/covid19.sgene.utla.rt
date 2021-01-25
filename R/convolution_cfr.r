# Packages ----------------------------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(vroom)
library(lubridate)
library(brms)
library(data.table)
library(loo)
library(future)
library(future.apply)

# set number of parallel cores
no_cores <- availableCores()

# Get functions -----------------------------------------------------------
source(here("R", "brm_convolution.r"))
source(here("R", "get_notifications_data.r"))

# Get data ----------------------------------------------------------------
df <- list()

# all data combinations by UTLA
df[["utla"]][["cfr"]] <- get_notifications_data("cases", "deaths")
df[["utla"]][["chr"]] <- get_notifications_data("cases", "admissions")
df[["utla"]][["hfr"]] <- get_notifications_data("admissions", "deaths")

# all data combinations by NHS region
df[["region"]][["cfr"]] <- get_notifications_data("cases", "deaths", level = "nhs region")
df[["region"]][["chr"]] <- get_notifications_data("cases", "admissions", level = "nhs region")
df[["region"]][["hfr"]] <- get_notifications_data("admissions", "deaths", level = "nhs region")

# Define model ------------------------------------------------------------
models <- list()

models[["intercept"]] <- as.formula(secondary ~ 1 + prop_sgtf)
models[["time"]] <- as.formula(secondary ~ (1 | time) + prop_sgtf)
models[["utla"]] <- as.formula(secondary ~ (1 | loc) + prop_sgtf)
models[["all"]] <- as.formula(secondary ~ (1 | loc) + (1 | time) + prop_sgtf)

# Fit models --------------------------------------------------------------
# set up parallel
## core usage
if (no_cores <= 4) {
  stan_cores <- no_cores
  mc_cores <- 1
} else {
  stan_cores <- 4
  mc_cores <- ceiling(no_cores / 4)
}
plan("multisession", workers = mc_cores, earlySignal = TRUE)

#define context specific args
fit_brm_convolution <- function(formula, ...) {
  brm_convolution(formula, control = list(adapt_delta = 0.99, max_treedepth = 12),
                  iter = 3000, cores = stan_cores, ...)
}

# set context specific priors (based on mean in data)
priors <- list()
priors[["cfr"]] <- c(prior("normal(-4, 0.5)", class = "Intercept"))
priors[["chr"]] <- c(prior("normal(-2.5, 0.5)", class = "Intercept"))
priors[["hfr"]] <- c(prior("normal(-1, 0.5)", class = "Intercept"))

# fit model grid in parallel
fit_targets <- expand_grid(loc = c("region", "utla"), conv = c("fixed", "loc"), 
                           target = c("cfr", "chr", "hfr"))

fits <- future_lapply(1:nrow(fit_targets), function(i) {
  ft <- fit_targets[i, ]
  message("Fitting ", ft$target, " at the ", ft$loc, " level using following convolution: ", ft$conv)
  out <- list()
  fits <- suppressMessages(lapply(models, fit_brm_convolution,
                data = df[[ft$loc]][[ft$target]],
                prior = priors[[ft$target]],
                conv_varying = ft$conv))
  ft$models <- list(names(models))
  ft$fit <- list(fits)
  ft <- unnest(ft, cols = c("models", "fit"))
  return(ft)},
  future.scheduling = Inf, future.seed = TRUE)

fits <- reduce(fits, bind_rows)

# Extract rate summary ----------------------------------------------------
fits <- fits %>% 
  mutate(variant_effect_q = map(fit, function(x) {
    samples <- posterior_samples(x, "b_prop_sgtf")
    q <- samples[, "b_prop_sgtf"]
    q <- exp(q)
    q <- quantile(q, c(0.025, 0.5, 0.975))
    q <- signif(q, 2)
    return(q)
  })) %>% 
  mutate(variant_effect =  map_chr(variant_effect_q,
                                   ~ paste0(.[2]," (", .[1], ", ", .[3], ")")))

saveRDS(fits, here("output", "convolution-associations.rds"))

