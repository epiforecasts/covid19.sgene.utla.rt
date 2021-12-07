# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)
library(readr)
library(magrittr)

# set to use either web version or local version of data\
use_web <- TRUE

# Download Rt estimates ---------------------------------------------------
week_start <- readRDS(here("data", "by_utla.rds")) %>%
  .$week_infection %>%
  subtract(7) %>%
  max() %>%
  wday()

# extract from epiforecasts.io/covid or use archived version
if (use_web) {
  rt_estimates <-
  paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
         "master/subnational/united-kingdom-local/cases/summary/rt.csv")
  rt <- vroom(rt_estimates)
  vroom_write(rt, here("data-raw", "rt.csv"), delim = ",")
} else{
  rt <- vroom(here("data-raw", "rt.csv"))
}

# join and save
rt <- rt %>%
  filter(type == "estimate") %>%
  select(utla_name = region, date, everything(), -strat, -type)

# Make Rt weekly
rt_weekly <- rt %>%
  mutate(week_infection =
           floor_date(date, "week", week_start = week_start) + 6) %>%
  group_by(utla_name, week_infection) %>%
  summarise(mean = mean(mean), sd = mean(sd), n = n(), .groups = "drop") %>%
  filter(n == 7) %>%
  select(-n) %>%
  pivot_longer(c(mean, sd)) %>%
  mutate(rt = paste("rt", name, sep = "_")) %>%
  select(-name) %>%
  pivot_wider(names_from = rt)

saveRDS(rt_weekly, here("data", "rt_weekly.rds"))
