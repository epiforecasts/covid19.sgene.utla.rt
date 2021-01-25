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

# set context specific priors
priors <- c(prior("normal(-4, 0.5)", class = "Intercept"))

# fit model
fit <- brm_convolution(secondary ~ (1 | loc) + prop_sgtf, data = df[["region"]][["hfr"]], 
                         prior = priors)
