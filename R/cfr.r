# Packages ----------------------------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(vroom)
library(lubridate)
library(brms)
library(parallel)
library(loo)

# set number of parallel cores
no_cores <- detectCores()
options(mc.cores = no_cores)
# Data --------------------------------------------------------------------
# get raw cases by data of infection from epiforecasts.io 
cases <- vroom(paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
                          "master/subnational/united-kingdom-local/cases/summary/reported_cases.csv")) %>% 
  mutate(data = "cases")
# get raw deaths by data of infection from epiforecasts.io 
deaths <- vroom(paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
                           "master/subnational/united-kingdom-local/deaths/summary/reported_cases.csv")) %>% 
  mutate(data = "deaths")

# link datasets and pivot wider
reports <- cases %>% 
  bind_rows(deaths) %>% 
  select(utla_name = region, date, data, value = confirm) %>% 
  pivot_wider(names_from = "data")

# Define covariates -------------------------------------------------------
# get variant proportion
sgene_by_utla <- readRDS(here("data", "sgene_by_utla.rds")) %>% 
  drop_na(prop_sgtf) %>% 
  mutate(week_specimen = week_infection + 7) %>% 
  filter(week_specimen > "2020-10-01") %>% 
  group_by(utla_name,week_specimen) %>% 
  slice()

# add week indicator
reports <- reports %>% 
  mutate(week_specimen = floor_date(date, "week", week_start = wday(max(sgene_by_utla$week_specimen)) - 1))

# link with variant proportion
secondary_with_cov <- reports %>% 
  inner_join(sgene_by_utla, by = c("week_specimen", "utla_name")) %>% 
  select(utla = utla_name, date, week_specimen, region = nhser_name, 
         deaths, cases, prop_sgtf, samples)

# make normalised predictors
secondary_with_cov <- secondary_with_cov %>% 
  mutate(normalised_cases = (cases - mean(cases)) / sd(cases),
         time = as.numeric(date),
         time = time - min(time),
         time = (time - mean(time)) / sd(time))

# filter to a reduced set for testing
secondary_with_cov <- secondary_with_cov %>%
  rename(loc = utla, primary = cases) 

# Define model ------------------------------------------------------------
source(here("R/convolution_model.r"))

# set context specific priors
priors <- c(prior("normal(-4, 1)", class = "Intercept"))

# fit model
fit <- convolution_model(deaths ~ (1 | loc) + prop_sgtf, data = secondary_with_cov, prior = priors)

