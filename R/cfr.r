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

##' Get a combined data set of two epidemiological indicators
##'
##' Indicators will be matched to the time of the "from" indicator, with a lag
##' given by \code{infections_lag}
##' @param from "cases", "admissions" or "deaths"
##' @param to "cases", "admissions" or "deaths"
##' @param type "deconvoluted" (via the EpiNow2 model) or "lagged" (fixed lag,
##' estimated from data via correlation analysis)
##' @return a data frame of the two indicators matched to the same time
##' @author Sebastian Funk
get_data <- function(from = c("cases", "admissions", "deaths"),
                     to = c("cases", "admissions", "deaths"),
                     type = c("backcalculated", "lagged"),
                     infections_lag = 7) {

  ## Arguments ---------------------------------------------------------------
  from <- match.arg(from, choices = c("cases", "admissions", "deaths"))
  to <- match.arg(to, choices = c("cases", "admissions", "deaths"))
  type <- match.arg(type, choices = c("backcalculated", "lagged"))

  ## Data --------------------------------------------------------------------
  base_url <-
    paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
           "master/subnational/united-kingdom-local")

  if (type == "backcalculated") {
    file_name <- "cases_by_infection.csv"
  } else if (type == "lagged") {
    file_name <- "reported_cases.csv"
  }

  infs <- lapply(c(from, to), function(x) {
    vroom(paste(base_url, x, "summary", file_name, sep = "/")) %>%
      mutate(data = x) %>%
      mutate(date = date - infections_lag)
  })

  ## if type is lagged, add estimate of the fixed lag
  if (type == "lagged") {
    all_dates <-
      as.Date(unique(unlist(lapply(infs, function(x) as.character(x$date)))))
    combined <- tibble(date = all_dates) %>%
      full_join(infs[[1]], by = "date") %>%
      full_join(infs[[2]], by = c("date", "region")) %>%
      filter(!is.na(confirm.x), !is.na(confirm.y)) %>%
      group_by(date, region) %>%
      summarise(cases = sum(confirm.x),
                deaths = sum(confirm.y), .groups = "drop") %>%
      ungroup()

    ## find lag that maximises mean correlation across UTLAs
    max <- 28
    cor <- vapply(seq(0, max), function(lag) {
      combined %>%
        group_by(region) %>%
        mutate(deaths = lead(deaths, n = lag)) %>%
        filter(!is.na(deaths)) %>%
        mutate(id = 1:n()) %>%
        filter(id <= (max(id) - max + lag)) %>%
        summarise(cor = cor(cases, deaths), .groups = "drop") %>%
        filter(!is.na(cor)) %>%
        summarise(mean = mean(cor)) %>%
        .$mean
    }, 0)

    set_lag <- which.max(cor) - 1
    infs[[2]] <- infs[[2]] %>%
      mutate(date = date - set_lag)

    infs <- lapply(infs, function(x) {
      x %>%
        rename(median = confirm) %>%
        mutate(type = "estimate")
    })
  } else {
    set_lag <- NULL
  }

  ## link datasets and pivot wider
  infections <- bind_rows(infs) %>%
    filter(type == "estimate") %>%
    select(utla_name = region, date, data, value = median) %>%
    pivot_wider(names_from = "data") %>%
    group_by(utla_name) %>%
    complete(date = seq(min(date), max(date), by = "day")) %>%
    drop_na()

  ## Define covariates -------------------------------------------------------
  ## get variant proportion
  sgene_by_utla <- readRDS(here("data", "sgene_by_utla.rds")) %>%
    drop_na(prop_sgtf) %>%
    filter(week_infection > "2020-10-01")

  week_start <- lubridate::wday(max(sgene_by_utla$week_infection)) - 1
  ## make infections weekly summary
  weekly_infections <- infections %>%
    mutate(week_infection =
             floor_date(date, "week", week_start = week_start)) %>%
    select(-date) %>%
    group_by(utla_name, week_infection) %>%
    add_count() %>%
    group_by(utla_name, week_infection, n) %>%
    summarise_all(sum) %>%
    ungroup() %>%
    filter(n == 7) %>%
    select(-n)

  ## link with variant proportion
  deaths_with_cov <- weekly_infections %>%
    inner_join(sgene_by_utla, by = c("week_infection", "utla_name")) %>%
    rename(utla = utla_name, region = nhser_name) %>%
    select(-utla_code)

  ## factorise time
  deaths_with_cov <- deaths_with_cov %>%
    mutate(time = factor((week_infection - min(week_infection)) / 7))

  if (!is.null(set_lag)) {
    deaths_with_cov <- deaths_with_cov %>%
      mutate(lag = set_lag)
  }

  ## rename everything to cases/deaths to make things easier later
  deaths_with_cov <- deaths_with_cov %>%
    rename(!!!c(cases = from, deaths = to))

  return(deaths_with_cov)
}

## Define model ------------------------------------------------------------
## D = c^{+}(1-f)C + c^{-}fC
## where c is the CFR, f is the fraction of cases that are sgtf, C is the number
## of cases by date of infection, and D is the number of deaths by date of
## infection. Effect of the variant on the CFR is assumed to be:
## c^{-} = \alpha c^{+} or c^{-} = \alpha + c^{+}
## leading to:
## D = (1 - (\alpha - 1)f)c^{+}C or D = (c^{+} + \alpha f)C

## define custom negative binomial family including variant factor and cases
variant_nb <- function(additive = FALSE) {
  custom_family(
    "variant_nb", dpars = c("mu", "phi", "alpha"),
    links = c("logit", "log", "identity"),
    lb = c(0, 0, ifelse(!additive, 0, NA)),
    type = "int",
    vars = c("f[n]", "cases[n]", "effect[1]")
  )
}

## define stan code to scale cfr by cases and variant fraction
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

## define model function
nb_model <- function(form, iter = 2000, data = deaths_with_cov,
                     additive = FALSE, ...) {
  # define priors
  if (additive) {
    priors <- c(prior(normal(0, 0.01), class = alpha))
  }else{
    priors <- c(prior(lognormal(0, 0.5), class = alpha))
  }
  message("Fitting: ", as.character(form))
  brm(formula = form,
      family = variant_nb(additive = additive),
      prior = priors,
      data,
      stanvars = make_stanvars(data, additive = additive),
      control = list(adapt_delta = 0.99),
      warmup = 1000, iter = iter, ...)
}

df <- list()
df[["cfr"]] <- get_data("cases", "deaths", "lagged")
df[["chr"]] <- get_data("cases", "admissions", "lagged")
admissions_lag <- unique(df[["chr"]]$lag)
df[["hfr"]] <- get_data("admissions", "deaths", "lagged",
                        infections_lag = 7 + admissions_lag)

## Fit models --------------------------------------------------------------
models <- list()
## define models to fit
models[["intercept"]] <- as.formula(deaths ~ 1)
models[["time"]] <- as.formula(deaths ~ time)
models[["utla"]] <- as.formula(deaths ~ (1 | utla))

## core usage
if (no_cores <= 4) {
  options(mc.cores = no_cores)
  mc_cores <- 1
} else {
  options(mc.cores = 4)
  mc_cores <- ceiling(no_cores / 4)
}

## fit models
warning("Fitting models sequentially due to mclapply issues")
results <- lapply(names(df), function(x) {
  fits <- list()
  cat(x, "multiplicative\n")
  fits[["multiplicative"]] <-
    lapply(models, nb_model, data = df[[x]])
  cat(x, "additive\n")
  fits[["additive"]] <-
    lapply(models, nb_model, data = df[[x]],
           additive = TRUE)
  return(fits)
})

## variant effect ----------------------------------------------------------
extract_variant_effect <- function(x, additive = FALSE) {
  samples <- posterior_samples(x, "alpha")
  q <- samples[, "alpha"]
  q <- quantile(q, c(0.025, 0.5, 0.975))
  q <- signif(q, 2)
  return(q)
}

var_res <- lapply(names(results), function(x) {
  variant_effect <- list()
  variant_effect[["multiplicative"]] <-
    lapply(fits[[x]][["multiplicative"]], extract_variant_effect)
  variant_effect[["additive"]] <-
    lapply(fits[[x]][["additive"]], extract_variant_effect, additive = TRUE)
  return(variant_effect)
})

## Compare models ----------------------------------------------------------
## requires custom log_lik functions to be implemented for variant_nb family
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

add_loo <- function(fits) {
  fits %>%
    lapply(add_criterion, "loo") %>%
    lapply(loo, save_psis = TRUE)
}

options(mc.cores = no_cores)
model_loos <- lapply(names(results), function(x)) {
  fits <- results[[x]]
  loos <- list()
  loos[["multiplicative"]] <- add_loo(fits[["multiplicative"]])
  loos[["additive"]] <- add_loo(fits[["additive"]])
  lc <- list()
  lc[["multiplicative"]] <- loo_compare(loos[["multiplicative"]])
  lc[["additive"]] <- loo_compare(loos[["additive"]])
  all_loos <- flatten(loos)
  names(all_loos) <- c(paste0(names(all_loos)[1:length(models)], "_multiplictive"),
                       paste0(names(all_loos)[1:length(models)], "_additive"))
  lc[["all"]] <- loo_compare(all_loos)
  return(list(loos = loos, lc = lc))
}

# Save results ------------------------------------------------------------
output <- list()
output$data <- deaths_with_cov
output$fits <- fits
output$effect <- variant_effect
output$loos <- loos
output$lc <- lc
saveRDS(output, here("output", "cfr.rds"))
