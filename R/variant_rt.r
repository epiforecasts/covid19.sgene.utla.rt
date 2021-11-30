#' brms model for estimating variant transmissibility
#' 
#' @description Fit an additive mixture of a baseline and variant reproduction
#' number (Rt) to estimated non-parametric Rt values. The variant Rt is
#' described by a multiplication of the baseline Rt. The baseline Rt is
#' described by any covariate combination supported by `brms`. 
#' 
#' Uncertainty is included in the observed Rt estimates  by assuming that they
#' have student T distribution with standard deviation derived as a combination
#' of a modelled estimate and the observed standard error of the estimates.
#' Uncertainty is also included by modelling variant samples (`positive_prop`) as a
#' latent variable informed by the number of overall samples available and the
#' number of samples which are postive for the variant of interest.
#' 
#' @param log_rt Formula describing the covariates  predicting the log of the
#' baseline variants time-varying reproduction number.
#' @param data Data frame that must contain the responses: `rt_mean`,
#' `rt_sd`, `positive_samples`, `samples` and a dummy variable `positive_prop` which 
#' is `NA_real`. If wanting to not model uncertainty in the variant proportion
#' fill this variable with the point estimate of the proportion.
#' @param brm_fn Function from `brms` to apply. Common options are
#' `brm()` for model fitting,`get_prior()` to see the priors (not may not
#' be accurate in all cases), and `make_stancode()` to see the generated
#' stan code.
#' @param ... Additional arguments to pass to `brm_fn`.
#' @import brms
variant_rt <- function(log_rt = ~ 1, data, brm_fn = brm, ...) {

  # define Rt variant mixture
  rt_model <-
    bf(rt_mean | se(rt_sd, sigma = TRUE) ~ (1 + var) * exp(logRt),
       nl = TRUE) +
    lf(var ~ 0 + mi(positive_prop)) +
    lf(as.formula(paste0("logRt ", paste(log_rt, collapse = " ")))) +
    student()

  # define samples observation model
  samples_obs_model <- bf(positive_samples | trials(samples) ~ 0 + mi(positive_prop)) +
    binomial(link = "identity")

  # define model for latent variant proportion
  variant_model <- bf(positive_prop | mi() ~ 1, family = "beta")

  # set weak priors
  # set positive proportion to be latent and observed directly with no modifier
  priors <- c(
    prior(student_t(3, 0, 0.5), nlpar = "var", resp = "rtmean"),
    prior(student_t(3, 0, 0.5), nlpar = "logRt", resp = "rtmean"),
    prior(constant(1), class = "phi", resp = "positive_samples"),
    prior(constant(1), class = "Intercept", resp = "positive_samples"),
    prior(constant(1), resp = "positivesamples")
  )

# fit model
brm_out <- do.call(brm_fn, c(
  list(
    mvbf(samples_obs_model, variant_model, rt_model, rescor = FALSE),
    data = data,
    prior = priors),
  ...)
)

return(brm_out)
}
