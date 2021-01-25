# Packages ----------------------------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(vroom)
library(lubridate)
library(brms)
library(data.table)
library(parallel)
library(loo)

# set number of parallel cores
no_cores <- detectCores()
options(mc.cores = no_cores)

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
## define models to fit
models[["intercept"]] <- as.formula(secondary ~ 1 + prop_sgtf)
models[["time"]] <- as.formula(secondary ~ (1 | time) + prop_sgtf)
models[["utla"]] <- as.formula(secondary ~ (1 | loc) + prop_sgtf)
models[["all"]] <- as.formula(secondary ~ (1 | loc) + (1 | time) + prop_sgtf)


# Fit models --------------------------------------------------------------
#define context specific args
fit_brm_convolution <- function(formula, ...) {
  brm_convolution(formula, control = list(adapt_delta = 0.95, max_treedepth = 12),
                  iter = 3000, ...)
}

mc_cores <- length(models)

# set context specific priors
priors <- c(prior("normal(-4, 0.5)", class = "Intercept"))

# fit model grid in parallel
fit_targets <- expand_grid(loc = c("region", "utla"), conv = c("fixed", "loc"), 
                           target = c("cfr", "chr", "hfr"))

fits <- mclapply(1:nrow(fit_targets), function(i) {
  ft <- fit_targets[i, ]
  out <- list()
  fit <- lapply(models, fit_brm_convolution,
                data = df[[ft$loc]][[ft$target]],
                conv_varying = ft$conv)
  names(fit) <- names(models)
  out[[ft$target]][[ft$loc]][[ft$conv]] <- fit
  return(out)}, mc.cores = mc_cores)

