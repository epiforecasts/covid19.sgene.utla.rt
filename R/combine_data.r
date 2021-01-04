# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)

# Download Rt estimates ---------------------------------------------------

# extract from epiforecasts.io/covid
rt_estimates <-
  paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
         "master/subnational/united-kingdom-local/cases/summary/rt.csv")
short_rt <- vroom(rt_estimates)
saveRDS(short_rt , here("data-raw", "rt-short-generation-time.csv"))

# load sensitivity analysis with longer generation time
long_rt <- vroom(here("data-raw", "rt-long-generation-time.csv"))

# join and save
rt <- short_rt %>% 
  mutate(generation_time = "short") %>% 
  bind_rows(long_rt %>% 
              mutate(generation_time = "long")) %>% 
  filter(type == "estimate") %>% 
  select(generation_time, ltla_name = region, date, everything(), -strat, -type)

saveRDS(rt, here("data", "rt.rds"))

# Load sgene data ---------------------------------------------------------
sgene_ltla <- readRDS(file.path(processed_path, "sgene_by_ltla.rds"))

week_start <- unique(wday(sgene_ltla$week_infection)) - 1

# Combine data sources ----------------------------------------------------

# make  Rt estimates weekly
rt_weekly <- rt %>%
  mutate(week_infection = floor_date(date, "week", week_start = week_start)) %>%
  group_by(ltla_name, week_infection) %>%
  summarise(mean = mean(mean), sd = mean(sd), n = n(), .groups = "drop") %>%
  filter(n == 7) %>%
  select(-n)




by_ltla_rt <- by_ltla_aggregate %>%
  inner_join(rt_weekly, by = c("week_infection", "ltla_name")) %>%
  select(week_infection, nhser_name, ltla_name, ltla_code, prop_variant,
         prop_variant_sd, samples, cases, rt_mean = mean, rt_sd = sd)
