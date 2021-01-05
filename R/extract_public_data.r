# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)
library(readr)

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
  select(generation_time, utla_name = region, date, everything(), -strat, -type)

saveRDS(rt, here("data", "rt.rds"))


# make data weekly
week_start <- wday(max(rt$date))

rt_weekly <- rt %>%
  mutate(week_infection = floor_date(date, "week", week_start = week_start)) %>%
  group_by(utla_name, week_infection, generation_time) %>%
  summarise(mean = mean(mean), sd = mean(sd), n = n(), .groups = "drop") %>%
  filter(n == 7) %>%
  select(-n)

saveRDS(rt_weekly, here("data", "rt_weekly.rds"))

