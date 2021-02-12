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
  rename(rt_mean = rt_mean_short_gt)

# get generation time summary parameters
g_mean <- 3.64
g_sd <- 3.08
k <- round((g_sd / g_mean)^2, 2)

# Define priors -----------------------------------------------------------
priors <- prior(constant(0), nlpar = "mG", coef = "Intercept") +
  prior(normal(0, 0.5), nlpar = "mG") +
  prior(normal(0, 0.1), nlpar = "r") +
  set_prior(paste0("constant(", k, ")"), nlpar = "k", coef = "Intercept") +
  set_prior(paste0("constant(", g_mean, ")"), nlpar = "G", coef = "Intercept")

# Define non-linear model -------------------------------------------------
nl_model <- as.formula(
  rt_mean ~ pow(1 + (k * r * (exp(mG) * G)), 1 / k)
  )

# Define model structure
form <- bf(nl_model,
           r ~ 1, 
           k ~ 1,
           G ~ 1,
           mG ~ 1,
           nl = TRUE)

fit <- brm(formula = form, family = gaussian(), 
           data = utla_rt_with_covariates, prior = priors)

print(fit, digits = 3)
