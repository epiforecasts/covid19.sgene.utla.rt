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
  filter(week_specimen > "2020-10-01")

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

# Define model ------------------------------------------------------------
# D = c^{+}(1-f)C + c^{-}fC
# where c is the CFR, f is the fraction of cases that are sgtf, C is the number
# of cases by date of infection, and D is the number of deaths by date of infection
# effect of the variant on the CFR is assumed to be:
# c^{-} = \alpha c^{+} or c^{-} = \alpha + c^{+}
# leading to:
# D = (1 - (\alpha - 1)f)c^{+}C or D = (c^{+} + \alpha f)C 
# define custom negative binomial family including variant factor and cases
variant_nb <- function(additive = FALSE) {
  custom_family(
    "variant_nb", dpars = c("mu", "phi", "alpha"),
    links = c("logit", "log", "identity"),
    lb = c(0, 0, ifelse(!additive, 0, NA)),
    type = "int",
    vars = c("f[n]", "cases[n]", "effect[1]")
  )
}

# define stan code to scale cfr by cases and variant fraction
make_stanvars <- function(data, additive = FALSE) {
  stan_funs <- "
real variant_nb_lpmf(int y, real mu, real phi, real alpha,
                     real f, int cases, int effect) {
    real scaled_cases;
    if (effect) {
      scaled_cases = (mu + alpha * f) * cases;
    }else {
      scaled_cases = (1 + (alpha - 1) * f) * mu * cases;
    }
    return  neg_binomial_2_lpmf(y | scaled_cases, phi);
                            }
real variant_nb_rng(int y, real mu, real phi, real alpha,
                    real f, int cases, int effect) {
    real scaled_cases;
    if (effect) {
      scaled_cases = (mu + alpha * f) * cases;
    }else {
      scaled_cases = (1 + (alpha - 1) * f) * mu * cases;
    }
    return  neg_binomial_2_rng(scaled_cases, phi);
                            }
"
  stanvars <- c(stanvar(block = "functions", scode = stan_funs),
                stanvar(block = "data",
                        scode = "  real f[N];",
                        x = data$prop_sgtf,
                        name = "f"),
                stanvar(block = "data",
                        scode = "  int cases[N];",
                        x = data$cases,
                        name = "cases"),
                stanvar(block = "data",
                        scode = "  int effect[1];",
                        x = array(as.numeric(additive)),
                        name = "effect")
  )
  return(stanvars)
}

# define model function
nb_model <- function(form, iter = 2000, data = secondary_with_cov,
                     additive = FALSE, ...) {
  # define priors
  if (additive) {
    priors <- c(prior(normal(0, 0.01), class = alpha))
  }else{
    priors <- c(prior(lognormal(0, 0.5), class = alpha))
  }
  message("Fitting: ", as.character(form))
  make_stancode(formula = form,
      family = variant_nb(additive = additive),
      prior = priors,
      data,
      stanvars = make_stanvars(data, additive = additive),
      control = list(adapt_delta = 0.99),
      warmup = 1000, iter = iter, ...)
}
# Fit models --------------------------------------------------------------
models <- list()
# define models to fit
models[["intercept"]] <- as.formula(deaths ~ 1)
models[["cases"]] <- as.formula(deaths ~ s(normalised_cases, k = 5))
models[["region"]] <- as.formula(deaths ~ region)
models[["time"]] <- as.formula(deaths ~ s(time, k = 9))
models[["utla"]] <- as.formula(deaths ~ (1 | utla))
models[["all"]] <- as.formula(deaths ~ s(normalised_cases, k = 5) + region + (1 | utla))
models[["all_with_residuals"]] <- as.formula(deaths ~ s(normalised_cases, k = 5) + s(time, k = 5) + region + (1 | utla))
models[["all_with_regional_residuals"]] <- as.formula(deaths ~ s(normalised_cases, k = 5) + s(time, k = 5, by = region) + (1 | utla))

# core usage
if (no_cores <= 4) { 
  options(mc.cores = no_cores)
  mc_cores <- 1
}else{
  options(mc.cores = 4)
  mc_cores <- ceiling(no_cores / 4) 
}
# fit models
warning("Fitting models sequentially due to mclapply issues")
fits <- list()
fits[["multiplicative"]] <- lapply(models, nb_model)
fits[["additive"]] <- lapply(models, nb_model, additive = TRUE)

# variant effect ----------------------------------------------------------
extract_variant_effect <- function(x, additive = FALSE) {  
  samples <- posterior_samples(x, "alpha")
  q <- samples[, "alpha"]
  if (!additive) {
    q <- q - 1
  }
  q <- quantile(q, c(0.025, 0.5, 0.975))
  q <- signif(q, 2)
  return(q)
}
variant_effect <- list()
variant_effect[["multiplicative"]] <- lapply(fits[["multiplicative"]], extract_variant_effect) 
variant_effect[["additive"]] <- lapply(fits[["additive"]], extract_variant_effect, additive = TRUE) 

# Compare models ----------------------------------------------------------
# requires custom log_lik functions to be implemented for variant_nb family
expose_functions(fits[[1]][[1]], vectorize = TRUE)

log_lik_variant_nb <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  alpha <- prep$dpars$alpha
  f <- prep$data$f[i]
  y <- prep$data$Y[i]
  cases <- prep$data$cases[i]
  effect <- prep$data$effect[1]
  variant_nb_lpmf(y, mu, phi, alpha, f, cases, effect)
}

posterior_predict_variant_nb <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  alpha <- prep$dpars$alpha
  f <- prep$data$f[i]
  y <- prep$data$Y[i]
  cases <- prep$data$cases[i]
  effect <- prep$data$effect[1]
  variant_nb_rng(mu, phi, alpha, f, cases, effect)
}

options(mc.cores = no_cores)
loos <- list()
add_loo <- function(fits) {
  fits %>% 
    lapply(add_criterion, "loo") %>%
    lapply(loo, save_psis = TRUE)
}
loos[["multiplicative"]] <- add_loo(fits[["multiplicative"]])
loos[["additive"]] <- add_loo(fits[["additive"]])
lc <- list()
lc[["multiplicative"]] <- loo_compare(loos[["multiplicative"]])
lc[["additive"]] <- loo_compare(loos[["additive"]])
all_loos <- flatten(loos)
names(all_loos) <- c(paste0(names(all_loos)[1:length(models)], "_multiplictive"), 
                     paste0(names(all_loos)[1:length(models)], "_additive"))
lc[["all"]] <- loo_compare(all_loos)
  
# Save results ------------------------------------------------------------
output <- list()
output$data <- deaths_with_cov
output$fits <- fits
output$effect <- variant_effect
output$loos <- loos
output$lc <- lc
saveRDS(output, here("output", "cfr.rds"))
