# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)

# Load data ---------------------------------------------------------------
rt <- readRDS(here("data", "rt.rds"))
sgene_by_utla <- readRDS(here("data", "sgene_by_utla.rds"))
tiers <- readRDS(here("data", "tiers.rds"))
mobility <- readRDS(here("data", "mobility.rds"))

# Combine data sources ----------------------------------------------------
week_start <- wday(min(sgene_by_utla$week_infection))

# Make Rt weekly
rt_weekly <- rt %>%
  mutate(week_infection = floor_date(date, "week", week_start = (week_start + 6) %% 7)) %>%
  group_by(utla_name, week_infection, generation_time) %>%
  summarise(mean = mean(mean), sd = mean(sd), n = n(), .groups = "drop") %>%
  filter(n == 7) %>%
  select(-n)

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
