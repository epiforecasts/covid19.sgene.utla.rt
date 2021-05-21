# Packages ----------------------------------------------------------------
library(vroom)
library(tidyr)
library(dplyr)
library(here)
library(lubridate)
library(readr)
library(janitor)
library(magrittr)

# Extract data ---------------------------------------------------
week_start <- readRDS(here("data", "sgene_by_utla.rds")) %>%
  .$week_infection %>%
  subtract(7) %>%
  max() %>%
  wday()

utlas <- readRDS(here("data", "sgene_by_utla.rds")) %>%
  .$utla_name %>%
  unique()

baseline <- tibble(
  date = c(as.Date("2020-10-05")),
  utla_name = utlas,
  tier = "none"
)

tier_1 <- tibble(
  date = c(as.Date("2020-10-12")),
  utla_name = utlas,
  tier_1 = "tier_1"
)

initial_tiers <- read_csv(here("data-raw", "tiers", "initial_tiers.csv")) %>%
  full_join(tier_1, by = c("date", "utla_name")) %>%
  mutate(tier = if_else(is.na(tier), tier_1, tier)) %>%
  select(-tier_1) %>%
  arrange(utla_name, date)

start_tiers <- baseline %>%
  bind_rows(initial_tiers) %>%
  mutate(date = as.Date(date))

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
  filter(date < "2020-11-05") %>%
  mutate(tier = if_else(tier %in% c("Medium", "1"), "tier_1", tier),
         tier = if_else(tier %in% c("High", "2"), "tier_2", tier),
         tier = if_else(tier %in% c("Very High", "3"), "tier_3", tier)) %>%
  left_join(ltla_utla, by = "laname") %>%
  mutate(utla = if_else(is.na(utla), match_utla, utla)) %>%
  select(-match_utla) %>%
  filter(!is.na(utla)) %>%
  mutate(utla = if_else(utla %in% c("Cornwall", "Isles of Scilly"),
                        "Cornwall and Isles of Scilly", utla),
         utla = if_else(utla %in% c("Hackney", "City of London"),
                        "Hackney and City of London", utla)) %>%
  group_by(date, utla_name = utla, tier) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(date, utla_name) %>%
  filter(n == max(n)) %>%
  ungroup() %>%
  select(-n)

national_lockdown <- tibble(
  date = c(as.Date("2020-11-05")),
  utla_name = utlas,
  tier = "national_lockdown_nov"
)

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
  select(date, tier, utla_name = ctyua20nm, nhs = rgn19nm) %>%
  mutate(utla_name = if_else(utla_name %in% c("Cornwall", "Isles of Scilly"),
                             "Cornwall and Isles of Scilly", utla_name),
         utla_name = if_else(utla_name %in% c("Hackney", "City of London"),
                             "Hackney and City of London", utla_name)) %>%
  group_by(date, nhs, utla_name, tier) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(date, nhs, utla_name) %>%
  filter(n == max(n)) %>%
  ungroup() %>%
  select(-n)

utla_nhs <- new_tiers %>%
  select(utla_name, nhs) %>%
  distinct()

new_tiers <- new_tiers %>%
  select(-nhs) %>%
  mutate(tier = if_else(grepl("^[0-9]+$", tier), paste0("tier_", tier),
                         tolower(sub(" ", "_", tier))))

# Join into single dataset ------------------------------------------------
tiers <- start_tiers %>%
  bind_rows(old_tiers) %>%
  bind_rows(national_lockdown) %>%
  bind_rows(new_tiers) %>%
  complete(date = seq(min(date), today(), by = "day"),
           utla_name = unique(utla_name)) %>%
  arrange(utla_name, date) %>%
  fill(tier, .direction = "down") %>%
  mutate(week_infection = floor_date(date, "week", week_start) - 1) %>%
  group_by(utla_name, week_infection, tier) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(utla_name, week_infection) %>%
  filter(n == max(n)) %>%
  ungroup() %>%
  select(-n)

saveRDS(tiers, here("data", "tiers.rds"))

# Extract mobility data ---------------------------------------------------
gm_lad <- read_csv(here("data-raw", "google_lad.csv"))
mobility_file <- here("data-raw", "2021_GB_Region_Mobility_Report.csv")

mobility <- read_csv(mobility_file) %>%
  filter(!is.na(sub_region_1)) %>%
  pivot_longer(ends_with("_percent_change_from_baseline"),
               names_to = "variable") %>%
  select(name = sub_region_1, date, variable, value) %>%
  mutate(variable = sub("_percent_change_from_baseline", "", variable)) %>%
  inner_join(gm_lad, by = "name") %>%
  select(date, laname = lad_nm,
         variable = variable, value = value) %>%
  inner_join(ltla_utla, by = "laname") %>%
  rename(utla_name = match_utla) %>%
  right_join(utla_nhs, by = c("utla_name")) %>%
  complete(date = seq(min(date, na.rm = TRUE),
                      max(date, na.rm = TRUE), by = "day"),
	   utla_name = unique(utla_name),
	   variable = unique(variable)) %>%
  filter(!is.na(date), !is.na(variable)) %>%
  group_by(date, nhs, variable) %>%
  mutate(median = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(value = if_else(is.na(value), median, value)) %>%
  fill(value, .direction = "downup") %>%
  select(-nhs, -median) %>%
  mutate(utla_name = if_else(utla_name %in% c("Cornwall", "Isles of Scilly"),
                             "Cornwall and Isles of Scilly", utla_name),
         utla_name = if_else(utla_name %in% c("Hackney", "City of London"),
                        "Hackney and City of London", utla_name)) %>%
  mutate(week_infection = floor_date(date, "week", week_start) + 6) %>%
  group_by(week_infection, utla_name, variable) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  group_by(variable) %>%
  mutate(value = (value - mean(value, na.rm = TRUE)) /
           sd(value, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = variable) %>%
  arrange(utla_name, week_infection)

saveRDS(mobility, here("data", "mobility.rds"))
