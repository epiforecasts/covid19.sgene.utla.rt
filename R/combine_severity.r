# Packages ----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(here)
library(purrr)
library(brms)

# Read in results ---------------------------------------------------------
convolution_severity <-  readRDS(here("output", "convolution_severity.rds")) %>% 
  mutate(effect_type = "multiplicative",
         convolution = case_when(conv %in% "fixed" ~ "global convolution",
                                 conv %in% "loc" ~ "local convolution")) %>% 
  mutate(variant_effect_q = map(fit, function(x) {
    samples <- posterior_samples(x, "b_prop_sgtf")
    q <- samples[, "b_prop_sgtf"]
    q <- exp(q)
    q <- quantile(q, c(0.025, 0.5, 0.975))
    q <- signif(q, 3)
    return(q)
  }))

lagged_severity <-  readRDS(here("output", "lagged_severity.rds")) %>% 
  mutate(convolution = "global lag", 
         effect_type = ifelse(effect_type %in% "mutliplicative", "multiplicative", effect_type)) %>% 
  mutate(variant_effect_q = map(fit, function(x) {
    samples <- posterior_samples(x, "alpha")
    q <- samples[, "alpha"]
    q <- quantile(q, c(0.025, 0.5, 0.975))
    q <- signif(q, 3)
    return(q)
  }))

# Score lagged models -----------------------------------------------------
lagged_score <- lagged_severity %>% 
  select(loc, effect_type, target, models, variant_effect, loo) %>% 
  group_by(loc, target) %>% 
  summarise(
    models = paste0(effect_type, "-", models),
    ranks = loo_compare(loo),
    .groups = "drop"
  )

saveRDS(lagged_score, here("data", "severity_models_scored.rds"))

# Join results ------------------------------------------------------------
severity <- convolution_severity %>% 
  bind_rows(lagged_severity) %>% 
  select(-conv) %>% 
  select(target, loc, convolution, effect_type, effect = variant_effect, 
         effect_numeric = variant_effect_q, everything()) %>% 
  mutate(convolution = factor(convolution,
                              levels = c("global lag", "global convolution",
                                         "local convolution")))

# Summarise effects -------------------------------------------------------
variant_effects <- severity %>%
  select(-fit, -data, -loo) %>% 
  arrange(loc, effect_type, convolution)

saveRDS(variant_effects, here("data", "severity_sgtf_effect.rds"))

# Save analysis data ------------------------------------------------------
lagged_severity_data <- severity %>% 
  select(target, loc, data) %>% 
  drop_na(data) %>% 
  unique()

saveRDS(lagged_severity_data, here("data", "lagged_severity_data.rds"))

# Save baseline ratio -----------------------------------------------------
baseline_severity_ratio <- severity %>%
  mutate(baseline = map_chr(fit, function(x) {
    samples <- posterior_samples(x, "b_Intercept")
    q <- samples[, "b_Intercept"]
    q <- exp(q)
    q <- quantile(q, c(0.025, 0.5, 0.975))
    q <- signif(q, 2)
    q <- paste0(q[2]," (", q[1], ", ", q[3], ")")
    return(q)
  })) %>% 
  select(-fit, -data, -loo) %>% 
  arrange(loc, effect_type, convolution)

saveRDS(baseline_severity_ratio, here("data", "baseline_severity_ratio.rds"))

