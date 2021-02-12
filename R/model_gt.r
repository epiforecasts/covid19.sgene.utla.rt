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
priors <- prior(constant(0), nlpar = "afk") +
  prior(constant(0), nlpar = "afG") +
  prior(constant(0), nlpar = "afr") +
  prior(student_t(3, 0, 0.25), nlpar = "rp") +
  set_prior(paste0("constant(", k, ")"), nlpar = "kp") +
  set_prior(paste0("constant(", g_mean, ")"), nlpar = "Gp")

# Define non-linear model -------------------------------------------------
nl_model <- as.formula(
  rt_mean ~ (1 + ((1 + afk) * kp) * ((1 + afr) * rp) * ((1 + afG) * Gp))^(1/((1 + afk) * kp))
  )

# Define parameters for model
form <- bf(nl_model,
           rp ~ 1,
           kp ~ 1,
           Gp ~ 1,
           afk ~ 1,
           afG ~ 1,
           afr ~ 1, 
           nl = TRUE)

fit <- brm(formula = form, family = gaussian(), 
           data = utla_rt_with_covariates, prior = priors)