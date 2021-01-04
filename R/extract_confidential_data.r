# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)
library(readr)

# Load Rt data ------------------------------------------------------------
rt <- readRDS(here("data", "rt.rds"))
# extract maximum date
week_start <- wday(max(rt$date))

# Extract tiers data ------------------------------------------------------
tiers_file <- here("data-raw", "NPI_dataset_full_extract_23_12_2020.csv")

tiers <- read_csv(tiers_file) %>%
  mutate(week_infection = floor_date(date, "week", week_start - 1),
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

saveRDS(here("data", "tiers.rds"))

# Extract mobility data ---------------------------------------------------
mobility_file <- here("data-raw", "google_mobility_2020-12-31.csv")

mobility <- read_csv() %>%
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

saveRDS(here("data", "mobility.rds"))

# Extract sgene data from Pillars data ------------------------------------
pillar_file <- here("data-raw", "english_pillars.rds")
english_pillars <- readRDS(pillar_file) %>%
  filter(between(date_specimen,
                 as.Date("2020-09-01"), max(date_specimen) - 4)) %>%
  mutate(week_specimen = floor_date(date_specimen, "week",
                                    week_start = week_start)) %>%
  group_by(week_specimen) %>%
  mutate(ndates = length(unique(date_specimen))) %>%
  ungroup() %>%
  filter(week_specimen < max(week_specimen) | ndates == 7) %>%
  select(-ndates)

by_ltla <- english_pillars %>%
  select(-lower_age_limit, -positive, -total) %>%
  group_by(pillar, week_specimen, nhser_name, ltla_name, ltla_code) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  pivot_longer(starts_with("sgene"), names_to = "sgene_result",
               values_to = "n") %>%
  mutate(sgene_result = sub("^sgene_", "", sgene_result))

by_ltla_aggregate <- by_ltla %>%
  filter(pillar == "Pillar 2") %>%
  mutate(week_infection = week_specimen - 7) %>%
  select(-pillar, -negative) %>%
  pivot_wider(names_from = sgene_result, values_from = n) %>%
  mutate(prop_variant = negative / (positive + negative),
         prop_variant_sd = sqrt(prop_variant * (1 - prop_variant) /
                                  (positive + negative)),
         samples = negative + positive,
         cases = n_a + negative + positive)

saveRDS(by_ltla_aggregate, here("data", "sgene_by_ltla.rds"))
