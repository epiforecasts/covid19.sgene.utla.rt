## Packages and options ----------------------------------------------------
library(bayesplot)
library(brms)
library(loo)
library(ggplot2)
library(dplyr)

options(mc.cores = 1)

compare_models <- function(file, type, models = NULL) {

  ## Load fits ---------------------------------------------------------------
  fits <- readRDS(file)
  if (!is.null(models)) {
    fits$models[[type]] <- fits$models[[type]][models]
  }

  filetype <- sub("sgene_fits_", "", basename(file))
  filetype <- paste0(sub("\\.rds", "", filetype), "_", type)

  fit_data <- fits$data[[type]] %>%
    filter(!is.na(prop_variant_sd))

  ## Add custom family functions ---------------------------------------------
  expose_functions(fits$models[[type]][[1]], vectorize = TRUE)

  log_lik_add_var_student <- function(i, prep) {
    mu <- prep$dpars$mu[, i]
    sigma <- prep$dpars$sigma
    nu <- prep$dpars$nu
    alpha <- prep$dpars$alpha
    frac <- prep$data$vreal1[i]
    y <- prep$data$Y[i]
    add_var_student_lpdf(y, mu, sigma, nu, alpha, frac)
  }

  posterior_predict_add_var_student <- function(i, prep, ...) {
    mu <- prep$dpars$mu[, i]
    sigma <- prep$dpars$sigma
    nu <- prep$dpars$nu
    alpha <- prep$dpars$alpha
    frac <- prep$data$vreal1[i]
    add_var_student_rng(mu, sigma, nu, alpha, frac)
  }

  ## Score fits --------------------------------------------------------------
  loos <- fits$models[[type]] %>%
    lapply(add_criterion, "loo") %>%
    lapply(loo, save_psis = TRUE)

  ## Model diagnostics -------------------------------------------------------
  y <- fit_data$rt_mean
  yrep <- lapply(fits$models[[type]], posterior_predict)

  psis <- lapply(fits$models[[type]], function(x) psis(-log_lik(x)))
  lw <- lapply(psis, weights)

  post <- lapply(fits$models[[type]],  function(x) {
    samples <- posterior_samples(x, "alpha")
    q <- quantile(samples[, "alpha"], c(0.025, 0.5, 0.975))
    q <- round(q, 2)
    return(q)
  })

  return(list(post = post,
              loos = loos,
              lc = lc,
              y = y,
              yrep = yrep,
              psis = psis))
}

gt <- c("short", "long")
res <- lapply(gt, function(x) {
  compare_models(here::here("output", paste0("sgene_fits_", x, "_gt.rds")),
                 "dynamic")
})
names(res) <- gt

saveRDS(post, here::here("output", "sgene_model_comparison.rds"))
