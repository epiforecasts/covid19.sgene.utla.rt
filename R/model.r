# Packages ----------------------------------------------------------------
library("brms")
library("here")
library("dplyr")
library("tidyr")
library("readr")
library("lubridate")

# Options -----------------------------------------------------------------
options(mc.cores = 4)

# Get data ----------------------------------------------------------------
processed_path <- here::here("data", "processed")

tiers_file <-
  here::here("data", "raw", "undated",
             "NPI_dataset_full_extract_23_12_2020.csv")

sgene_ltla <- readRDS(file.path(processed_path, "sgene_by_ltla.rds"))

week_start <- unique(wday(sgene_ltla$week_infection)) - 1

rt_estimates <-
  paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
         "master/subnational/united-kingdom-local/cases/summary/rt.csv")
rt <- vroom::vroom(rt_estimates)

rt_by_ltla <- rt %>%
  rename(ltla_name = region) %>%
  filter(type == "estimate")

rt_weekly <- rt_by_ltla %>%
  mutate(week_infection = floor_date(date, "week", week_start = week_start)) %>%
  group_by(ltla_name, week_infection) %>%
  summarise(mean = mean(mean), sd = mean(sd), n = n(), .groups = "drop") %>%
  filter(n == 7) %>%
  select(-n)

ltla_rt <- sgene_ltla %>%
  inner_join(rt_weekly, by = c("week_infection", "ltla_name")) %>%
  select(week_infection, nhser_name, ltla_name, ltla_code, prop_variant,
         prop_variant_sd, samples, cases, rt_mean = mean, rt_sd = sd) %>%
  mutate(time = as.numeric((week_infection - min(week_infection)) / 7))

start_of_week <- unique(wday(ltla_rt$week_infection))

tiers <- read_csv(tiers_file) %>%
  mutate(week_infection = floor_date(date, "week", start_of_week - 1),
         none = if_else(national_lockdown == 1, 0,
                             1 - tier_1 - tier_2 - tier_3 - tier_4)) %>%
  pivot_longer(c(none, starts_with("tier_"), national_lockdown),
               names_to = "tier") %>%
  group_by(ltla_name = ltla, week_infection, tier) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  group_by(ltla_name, week_infection) %>%
  filter(value == max(value)) %>%
  ungroup() %>%
  select(-value)

mobility <- read_csv(here::here("data", "raw", "undated",
                                "google_mobility_2020-12-31.csv")) %>%
  select(week_infection = date, ltla_name = name, region = nhs_nm,
         variable, value) %>%
  inner_join(ltla_rt %>% select(ltla_name, week_infection),
             by = c("ltla_name", "week_infection")) %>%
  group_by(week_infection, region, variable) %>%
  mutate(median = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(value = if_else(is.na(value), median, value)) %>%
  select(-region, -median) %>%
  group_by(variable) %>%
  mutate(value = (value - mean(value)) / sd(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = variable)

ltla_rt_interventions <- ltla_rt %>%
  inner_join(tiers, by = c("week_infection", "ltla_name")) %>%
  inner_join(mobility, by = c("week_infection", "ltla_name")) %>%
  filter(!is.na(prop_variant_sd))

# add small amount of noise to 0 measurement error
ltla_rt_interventions <- ltla_rt_interventions %>%
  mutate(prop_variant_sd = ifelse(prop_variant_sd == 0, 1e-5, prop_variant_sd))

static_rt <- ltla_rt_interventions %>%
  filter(week_infection == max(week_infection))

# Add custom family -------------------------------------------------------
add_var_student <- custom_family(
  "add_var_student", dpars = c("mu", "sigma", "nu", "alpha"),
  links = c("log", "identity", "identity", "identity"),
  lb = c(NA, 0, 1, NA),
  type = "real",
  vars = "vreal1[n]"
)

stan_funs <- "
real add_var_student_lpdf(real y, real mu, real sigma, real nu, real alpha,
                          real f) {
    real combined_mu = (1 + alpha * f) * mu;
    return student_t_lpdf(y | nu, combined_mu, sigma);
                            }
real add_var_student_rng(real mu, real sigma, real nu, real alpha, real f) {
    real combined_mu = (1 + alpha * f) * mu;
    return student_t_rng(nu, combined_mu, sigma);
  }
"
stanvars <- stanvar(block = "functions", scode = stan_funs)

# Set up shared priors ----------------------------------------------------
priors <- c(prior(gamma(2, .1), class = nu),
            prior(student_t(3, 0, 0.5), class = alpha),
            prior(student_t(3, 0, 0.5), class = sigma))

# Set up model ------------------------------------------------------------
base_model <- function(form, iter = 2000, ...) {
  brm(formula = form,
      family = add_var_student,
      stanvars = stanvars, 
      warmup = 500, iter = iter, ...)
}
## Static model ------------------------------------------------------------
# set model settings and priors
static_model <- function(form, ...) {
   base_model(form = form, data = static_rt, 
              control = list(adapt_delta = 0.95), ...)
}

static <- list()
static[["intercept"]] <- static_model(rt_mean | vreal(prop_variant) ~ 1,
                                      prior = priors)
static[["region"]] <- static_model(rt_mean | vreal(prop_variant) ~ nhser_name,
                                   prior = c(priors,
                                             prior(student_t(3, 0, 0.5), class = "b")))

# Dynamic model -----------------------------------------------------------
# set model settings and priors
dynamic_model <- function(form, iter = 2000, ...) {
  base_model(form = form, data = ltla_rt_interventions,
             prior = c(priors,
                       prior(student_t(3, 0, 0.5), class = "b")),
             control = list(adapt_delta = 0.95, max_treedepth = 12),
             iter = iter, ...)
}

# fit models
dynamic <- list()
dynamic[["interventions_only"]] <- 
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier)

dynamic[["interventions"]] <- 
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_random"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + (1 | ltla_name) + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + nhser_name + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_random_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier +  (1 | ltla_name) + nhser_name + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_time_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + s(time, k = 10) + nhser_name + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_time_by_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + s(time, k = 10, by = nhser_name) + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_time_by_random_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + s(time, k = 10, by = nhser_name) + 
                 (1 | ltla_name) + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_independent_time_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + factor(time):nhser_name + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

dynamic[["interventions_independent_time_random_region"]] <-
  dynamic_model(rt_mean | vreal(prop_variant) ~ tier + factor(time):nhser_name +
                  (1 | ltla_name) + 
		retail_and_recreation + grocery_and_pharmacy + 
		parks + transit_stations + workplaces + residential)

# Save fits ---------------------------------------------------------------
output_path <- here::here("output")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

saveRDS(list(models = list(static = static, dynamic = dynamic),
             data = list(static = static_rt, dynamic = ltla_rt_interventions)),
        file.path(output_path, "sgene_fits.rds"))
