# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)
library(covidregionaldata)
library(magrittr)

week_start <- readRDS(here("data", "sgene_by_utla.rds")) %>%
  .$week_infection %>%
  subtract(7) %>%
  max() %>%
  wday()

# Load data ---------------------------------------------------------------
rt_weekly <- readRDS(here("data", "rt_weekly.rds"))
sgene_by_utla <- readRDS(here("data", "sgene_by_utla.rds"))
tiers <- readRDS(here("data", "tiers.rds"))
mobility <- readRDS(here("data", "mobility.rds"))
utla_cases <- readRDS(here("data", "utla_cases.rds"))

# Combine data sources ----------------------------------------------------
utla_rt <- sgene_by_utla %>%
  inner_join(rt_weekly, by = c("week_infection", "utla_name")) %>%
  select(week_infection, nhser_name, utla_name,
         prop, samples, starts_with("rt_")) %>%
  mutate(time = as.numeric((week_infection - min(week_infection)) / 7))

tiers <- tiers %>%
	complete(week_infection = unique(utla_rt$week_infection),
		       utla_name = unique(utla_name)) %>%
	replace_na(list(tier = "none")) %>%
	arrange(utla_name, week_infection)

cases <- utla_cases %>%
  filter(date >= min(utla_rt$week_infection) + 7) %>%
  mutate(week_infection = floor_date(date, "week", week_start) + 6 - 7) %>%
  group_by(week_infection, utla_name = authority) %>%
  summarise(cases = sum(cases_new, na.rm = TRUE), .groups = "drop")

cases_last_4_weeks <- utla_rt %>%
  inner_join(cases, by = c("utla_name", "week_infection")) %>%
  filter(week_infection > max(week_infection) - 28) %>%
  group_by(utla_name, nhser_name) %>%
  summarise(new_samples = sum(samples),
            new_cases = sum(cases),
            .groups = "drop") %>%
  mutate(sampling = new_samples / new_cases)

utla_rt_with_covariates <- utla_rt %>%
  inner_join(tiers, by = c("week_infection", "utla_name")) %>%
  inner_join(mobility, by = c("week_infection", "utla_name")) %>%
  inner_join(cases_last_4_weeks,
             by = c("utla_name", "nhser_name")) %>%
  mutate(nhser_name = factor(nhser_name)) %>%
  filter(samples >= 20, sampling >= 0.2)

saveRDS(utla_rt_with_covariates, here("data", "utla_rt_with_covariates.rds"))
