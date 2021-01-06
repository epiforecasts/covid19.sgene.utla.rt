# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)
library(readr)
library(janitor)

# Extract data ---------------------------------------------------
week_start <- 1

old_tier_files <-
  list.files(here("data-raw", "tiers"), "uk_tier_data_parliament")

old_tiers <- lapply(old_tier_files, function(x) {
  df <- read_csv(here("data-raw", "tiers", x)) %>%
    mutate(date = sub("^.*_2020_([0-9]+)_([0-9]+)_.*$", "2020-\\1-\\2", x),
           date = as.Date(date))
  colnames(df) <- sub("^l_", "", colnames(df))
  if (!("tier" %in% colnames(df))) {
    df <- df %>%
      mutate(tier = "n_a")
  } else {
    df <- df %>%
      mutate(tier = as.character(tier))
  }
}) %>%
  bind_rows() %>%
  clean_names() %>%
  mutate(laname = if_else(is.na(laname), category, laname))

ltla_utla <- old_tiers %>%
  filter(!is.na(utla)) %>%
  select(laname, match_utla = utla) %>%
  distinct()

old_tiers <- old_tiers %>%
  filter(country == "E") %>%
  mutate(tier = if_else(tier %in% c("Medium", "1"), "tier_1", tier),
         tier = if_else(tier %in% c("High", "2"), "tier_2", tier),
         tier = if_else(tier %in% c("Very High", "3"), "tier_3", tier),
         tier = if_else(date > "2020-11-04", "national_lockdown", tier)) %>%
  left_join(ltla_utla, by = "laname") %>%
  mutate(utla = if_else(is.na(utla), match_utla, utla)) %>%
  select(-match_utla) %>%
  complete(date = seq(min(date), as.Date("2020-12-01"), by = "day"),
           utla = unique(utla)) %>%
  arrange(utla, date) %>%
  group_by(utla) %>%
  mutate(tier = if_else(is.na(tier) & date == min(date), "tier_1", tier)) %>%
  fill(tier, .direction = "down") %>%
  filter(!is.na(utla)) %>%
  select(date, utla, tier)

# Extract recent tiers ----------------------------------------------------
new_tier_files <-
  list.files(here("data-raw", "tiers"), "^England_LAD")

new_tiers <- lapply(new_tier_files, function(x) {
  df <- read_csv(here("data-raw", "tiers", x)) %>%
    mutate(date = sub("^.*_([0-9]+)_([0-9]+)_(202[0-9]?).*$", "\\3-\\2-\\1", x),
           date = as.Date(date)) %>%
    mutate(Tier = as.character(Tier))
}) %>%
  bind_rows() %>%
  clean_names() %>%
  select(date, tier, utla = ctyua20nm, nhs = rgn19nm)


utla_nhs <- new_tiers %>%
  select(utla, nhs) %>%
  distinct()

new_tiers <- new_tiers %>%
  mutate(tier = if_else(tier == "National Lockdown", "national_lockdown",
                        paste0("tier_", tier))) %>%
  complete(date = seq(min(date), today(), by = "day"), utla = unique(utla)) %>%
  fill(tier, .direction = "down")

# Join into single dataset ------------------------------------------------
tiers <- bind_rows(old_tiers, new_tiers) %>%
  rename(utla_name = utla) %>%
  arrange(utla_name, date) %>%
  mutate(utla_name = if_else(utla_name %in% c("Cornwall", "Isles of Scilly"),
                             "Cornwall and Isles of Scilly", utla_name),
         utla_name = if_else(utla_name %in% c("Hackney", "City of London"),
                             "Hackney and City of London", utla_name)) %>%
  mutate(week_infection = floor_date(date, "week", week_start)) %>%
  group_by(utla_name, week_infection, tier) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(utla_name, week_infection) %>%
  filter(n == max(n)) %>%
  ungroup() %>%
  select(-n)

saveRDS(tiers, here("data", "tiers.rds"))

# Extract mobility data ---------------------------------------------------
mobility_file <- here("data-raw", "google_mobility_2020-12-31.csv")

mobility <- read_csv(mobility_file) %>%
  select(date, laname = lad_nm,
         variable = variable, value = value) %>%
  inner_join(ltla_utla, by = "laname") %>%
  rename(utla = match_utla) %>%
  right_join(utla_nhs, by = c("utla")) %>%
  complete(date = seq(min(date, na.rm = TRUE), max(date, na.rm = TRUE), by = "day"), 
	   utla = unique(utla), 
	   variable = unique(variable)) %>%
  filter(!is.na(date)) %>%
  group_by(date, nhs, variable) %>%
  mutate(median = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(value = if_else(is.na(value), median, value)) %>%
  fill(value, .direction = "downup") %>%
  select(-nhs, -median) %>%
  mutate(utla = if_else(utla %in% c("Cornwall", "Isles of Scilly"),
                        "Cornwall and Isles of Scilly", utla),
         utla = if_else(utla %in% c("Hackney", "City of London"),
                        "Hackney and City of London", utla)) %>%
  mutate(week_infection = floor_date(date, "week", week_start)) %>%
  group_by(week_infection, utla_name = utla, variable) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  group_by(variable) %>%
  mutate(value = (value - mean(value, na.rm = TRUE)) /
           sd(value, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = variable) %>%
  arrange(utla_name, week_infection)

saveRDS(mobility, here("data", "mobility.rds"))
