# Packages ----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(here)
library(purrr)

# Read in results ---------------------------------------------------------
convolution_severity <-  readRDS(here("output", "convolution_severity.rds")) %>% 
  mutate(effect_type = "multiplicative",
         convolution = case_when(conv %in% "fixed" ~ "global convolution",
                                 conv %in% "loc" ~ "local convolution"))
lagged_severity <-  readRDS(here("output", "lagged_severity.rds")) %>% 
  mutate(convolution = "global lag", 
         effect_type = ifelse(effect_type %in% "mutliplicative", "multiplicative", effect_type))

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
  select(-conv, -variant_effect_q) %>% 
  select(target, loc, convolution, effect_type, effect = variant_effect, 
         everything()) %>% 
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

