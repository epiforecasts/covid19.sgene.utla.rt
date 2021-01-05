# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)

# Load data ---------------------------------------------------------------
rt_weekly <- readRDS(here("data", "rt_weekly.rds"))
sgene_by_utla <- readRDS(here("data", "sgene_by_utla.rds"))
tiers <- readRDS(here("data", "tiers.rds"))
mobility <- readRDS(here("data", "mobility.rds"))

# Combine data sources ----------------------------------------------------
utla_rt <- sgene_by_utla %>%
  inner_join(rt_weekly, by = c("week_infection", "utla_name")) %>%
  select(week_infection, generation_time, nhser_name, utla_name, utla_code,
         prop_variant, prop_variant_sd, samples, cases, rt_mean = mean, 
         rt_sd = sd) %>%
  mutate(time = as.numeric((week_infection - min(week_infection)) / 7))

utla_rt_with_covariates <- utla_rt %>%
  inner_join(tiers, by = c("week_infection", "utla_name")) %>%
  inner_join(mobility, by = c("week_infection", "utla_name"))

saveRDS(utla_rt_with_covariates, here("data", "utla_rt_with_covariates.rds"))
