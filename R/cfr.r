# Packages ----------------------------------------------------------------
library(dplyr)
library(vroom)
library(brms)

# Data --------------------------------------------------------------------
# get raw cases by data of infection from epiforecasts.io 
case_infs <- vroom(paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
                          "master/subnational/united-kingdom-local/cases/summary/cases_by_infection.csv")) %>% 
  mutate(data = "cases")
# get raw deaths by data of infection from epiforecasts.io 
death_infs <- vroom(paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
                           "master/subnational/united-kingdom-local/deaths/summary/cases_by_infection.csv")) %>% 
  mutate(data = "deaths")

infections <- case_infs %>% 
  bind_rows(death_infs) %>% 
  filter(type == "estimate") %>% 
  select(utla = region, date, data, value = median) %>% 
  pivot_wider(names_from = "data") %>% 
  group_by(utla) %>% 
  complete(date = seq(min(date), max(date), by = "day")) %>% 
  mutate(cases = replace_na(cases, 0)) %>% 
  drop_na(deaths)

# Define models -----------------------------------------------------------


# Fit models --------------------------------------------------------------


# Compare models ----------------------------------------------------------


