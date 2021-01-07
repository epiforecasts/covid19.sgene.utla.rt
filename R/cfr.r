# Packages ----------------------------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(vroom)
library(lubridate)
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

# link datasets and pivot wider
infections <- case_infs %>% 
  bind_rows(death_infs) %>% 
  filter(type == "estimate") %>% 
  select(utla_name = region, date, data, value = median) %>% 
  pivot_wider(names_from = "data") %>% 
  group_by(utla_name) %>% 
  complete(date = seq(min(date), max(date), by = "day")) %>% 
  mutate(cases = replace_na(cases, 0)) %>% 
  drop_na(deaths)

# Define covariates -------------------------------------------------------
# get variant proportion
sgene_by_utla <- readRDS(here("data", "sgene_by_utla.rds"))

# make infections weekly summary
weekly_infections <- infections %>% 
  mutate(week_infection = floor_date(date, "week", week_start = wday(max(sgene_by_utla$week_infection)) - 1)) %>% 
  group_by(utla_name, week_infection) %>%
  summarise(cases = sum(cases), deaths = sum(deaths), n = n(), .groups = "drop") %>%
  filter(n == 7) %>%
  select(-n)

# link with variant proportion
deaths_with_cov <- weekly_infections %>% 
  inner_join(sgene_by_utla, by = c("week_infection", "utla_name")) %>% 
  select(utla = utla_name, week_infection, region = nhser_name, 
         deaths, cases, prop_variant, samples)

# make normalised predictors
deaths_with_cov <- deaths_with_cov %>% 
  mutate(normalised_cases = (cases - mean(cases)) / sd(cases),
         time = as.numeric(week_infection),
         time = time - min(time),
         time = (time - mean(time)) / sd(time))

# Define models -----------------------------------------------------------

variant_nb <- custom_family(
  "variant_nb", dpars = c("mu", "phi", "alpha"),
  links = c("logit", "identity", "identity"),
  lb = c(0, 0, 0),
  type = "real",
  vars = "vreal1[n]"
)

stan_funs <- "
real variant_nb_lpdf(int y, real mu, real phi, real alpha,
                     int cases) {
    real scaled_cases = (1 + (alpha - 1)) * mu * cases;
    real sqrt_phi = 1 / sqrt(phi);
    return  neg_binomial_2_lpmf(y | scaled_cases, sqrt_phi);
                            }
real variant_nb_rng(int y, real mu, real phi, real alpha,
                     int cases) {
    real scaled_cases = (1 + (alpha - 1)) * mu * cases;
    real sqrt_phi = 1 / sqrt(phi);
    return  neg_binomial_2_rng(y | scaled_cases, sqrt_phi);
                            }
"
stanvars <- stanvar(block = "functions", scode = stan_funs)

# Set up shared priors ----------------------------------------------------
priors <- c(prior(lognormal(0, 1), class = alpha))

# Fit models --------------------------------------------------------------

base_model <- function(form, iter = 2000, ...) {
  brm(formula = form,
      family = variant_nb,
      stanvars = stanvars, 
      warmup = 500, iter = iter, ...)
}

base_model(deaths | vreal(prop_variant) ~ 1)
# Compare models ----------------------------------------------------------


