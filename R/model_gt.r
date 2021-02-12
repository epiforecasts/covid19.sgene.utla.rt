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
  rename(rt_mean = rt_mean_long_gt)


# Define non-linear model -------------------------------------------------
nl_model <- as.formula(
  rt_mean ~ (1 + ((1 + afk) * kp) * ((1 + afr) * rp) * ((1 + afG) * Gp))^(1/((1 + afk) * kp))
  )
form <- bf(nl_model,
           rp ~ 1,
           kp ~ 1,
           Gp ~ 1,
           afk ~ 1,
           afG ~ 1,
           afr ~ 1, 
           nl = TRUE)

fit <- make_stancode(formula = form, family = gaussian(), 
           data = utla_rt_with_covariates)
