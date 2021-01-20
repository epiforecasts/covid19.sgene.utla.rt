# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)
library(readr)

# Download Rt estimates ---------------------------------------------------
week_start <- 1

# extract from epiforecasts.io/covid
rt_estimates <-
  paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
         "master/subnational/united-kingdom-local/cases/summary/rt.csv")
short_rt <- vroom(rt_estimates)
vroom_write(short_rt, here("data-raw", "rt-short-generation-time.csv"),
            delim = ",")

# load sensitivity analysis with longer generation time
long_rt <- vroom(here("data-raw", "rt-long-generation-time.csv"))

# join and save
rt <- short_rt %>% 
  mutate(generation_time = "short") %>% 
  bind_rows(long_rt %>% 
              mutate(generation_time = "long")) %>% 
  filter(type == "estimate") %>% 
  select(generation_time, utla_name = region, date, everything(), -strat, -type)

# Make Rt weekly
rt_weekly <- rt %>%
  mutate(week_infection = floor_date(date, "week", week_start = week_start)) %>%
  group_by(utla_name, week_infection, generation_time) %>%
  summarise(mean = mean(mean), sd = mean(sd), n = n(), .groups = "drop") %>%
  filter(n == 7) %>%
  select(-n) %>%
  pivot_longer(c(mean, sd)) %>%
  mutate(gt_rt = paste("rt", name, generation_time, "gt", sep = "_")) %>%
  select(-generation_time, -name) %>%
  pivot_wider(names_from = gt_rt)

saveRDS(rt_weekly, here("data", "rt_weekly.rds"))
