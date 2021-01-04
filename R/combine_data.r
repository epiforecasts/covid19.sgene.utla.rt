# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)

# Load data ---------------------------------------------------------------
rt_weekly <- readRDS(here("data", "rt_weekly.rds"))
sgene_by_ltla <- readRDS(here("data", "sgene_by_ltla.rds"))
tiers <- readRDS(here("data", "tiers.rds"))
mobility <- readRDS(here("data", "mobility.rds"))

# Combine data sources ----------------------------------------------------
ltla_rt <- sgene_by_ltla %>%
  inner_join(rt_weekly, by = c("week_infection", "ltla_name")) %>%
  select(week_infection, generation_time, nhser_name, ltla_name, ltla_code,
         prop_variant, prop_variant_sd, samples, cases, rt_mean = mean, 
         rt_sd = sd) %>%
  mutate(time = as.numeric((week_infection - min(week_infection)) / 7))

ltla_rt_with_covariates <- ltla_rt %>%
  inner_join(tiers, by = c("week_infection", "ltla_name")) %>%
  inner_join(mobility, by = c("week_infection", "ltla_name")) %>%
  filter(!is.na(prop_variant_sd))

saveRDS(ltla_rt_with_covariates, here("data", "ltla_rt_with_covariates.rds"))