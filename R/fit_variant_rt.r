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
