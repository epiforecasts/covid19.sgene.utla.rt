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
  mutate(week_infection = floor_date(date, "week", week_start),
         none = if_else(national_lockdown == 1, 0,
                        1 - tier_1 - tier_2 - tier_3 - tier_4)) %>%
  pivot_longer(c(none, starts_with("tier_"), national_lockdown),
               names_to = "tier") %>%
  group_by(utla_name = utla20nm, week_infection, tier) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  group_by(utla_name, week_infection) %>%
  filter(value == max(value)) %>%
  ungroup() %>%
  select(-value) %>%
  mutate(utla_name = if_else(utla_name %in% c("Cornwall", "Isles of Scilly"),
                             "Cornwall and Isles of Scilly", utla_name),
         utla_name = if_else(utla_name %in% c("Hackney", "City of London"),
                             "Hackney and City of London", utla_name))

saveRDS(tiers, here("data", "tiers.rds"))

# Extract mobility data ---------------------------------------------------
mobility_file <- here("data-raw", "google_mobility_2020-12-31.csv")

mobility <- read_csv(mobility_file) %>%
  select(week_infection = date, ltla_name = lad_nm, ltla_code = lad_cd,
         variable = variable, value = value)

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
  select(-ndates) %>%
  mutate(utla_name = if_else(utla_name %in% c("Cornwall", "Isles of Scilly"),
                             "Cornwall and Isles of Scilly", utla_name),
         utla_name = if_else(utla_name %in% c("Hackney", "City of London"),
                             "Hackney and City of London", utla_name))

ltla_to_utla <- english_pillars %>%
  select(ltla_code, utla_code) %>%
  distinct()

utla_complete <- english_pillars %>%
  select(utla_code, utla_name, nhser_name) %>%
  distinct() %>%
  expand_grid(mobility %>%
              select(week_infection, variable) %>%
              distinct())

mobility <- mobility %>%
  inner_join(ltla_to_utla, by = "ltla_code") %>%
  right_join(utla_complete, by = c("utla_code", "week_infection", "variable")) %>%
  group_by(week_infection, nhser_name, variable) %>%
  mutate(median = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(value = if_else(is.na(value), median, value)) %>%
  select(-nhser_name, -median) %>%
  group_by(week_infection, utla_name, variable) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  group_by(variable) %>%
  mutate(value = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = variable)

saveRDS(mobility, here("data", "mobility.rds"))

by_utla <- english_pillars %>%
  select(-lower_age_limit, -positive, -total) %>%
  group_by(pillar, week_specimen, nhser_name, utla_name, utla_code) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  pivot_longer(starts_with("sgene"), names_to = "sgene_result",
               values_to = "n") %>%
  mutate(sgene_result = sub("^sgene_", "", sgene_result))

by_utla_aggregate <- by_utla %>%
  filter(pillar == "Pillar 2") %>%
  mutate(week_infection = week_specimen - 7) %>%
  select(-pillar, -negative) %>%
  pivot_wider(names_from = sgene_result, values_from = n) %>%
  mutate(prop_variant = negative / (positive + negative),
         prop_variant_sd = sqrt(prop_variant * (1 - prop_variant) /
                                  (positive + negative)),
         samples = negative + positive,
         cases = n_a + negative + positive) %>%
  
saveRDS(by_utla_aggregate, here("data", "sgene_by_utla.rds"))
