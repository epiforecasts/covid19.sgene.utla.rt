library("here")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")
library("vroom")
library("lubridate")

# make output directory
fig_path <- here::here("figure")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

rt_estimates <-
  paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
         "master/subnational/united-kingdom-local/cases/summary/rt.csv")
rt <- vroom::vroom(rt_estimates)

rt_by_ltla <- rt %>%
  rename(ltla_name = region) %>%
  filter(type == "estimate")

week_start <- wday(max(rt_by_ltla$date))

rt_weekly <- rt_by_ltla %>%
  mutate(week_infection = floor_date(date, "week", week_start = week_start)) %>%
  group_by(ltla_name, week_infection) %>%
  summarise(mean = mean(mean), sd = mean(sd), n = n(), .groups = "drop") %>%
  filter(n == 7) %>%
  select(-n)

processed_path <- here::here("data", "processed")
english_pillars <- readRDS(file.path(processed_path, "english_pillars.rds")) %>%
  filter(between(date_specimen,
                 as.Date("2020-09-01"), max(date_specimen) - 4)) %>%
  mutate(week_specimen = floor_date(date_specimen, "week",
                                    week_start = week_start)) %>%
  group_by(week_specimen) %>%
  mutate(ndates = length(unique(date_specimen))) %>%
  ungroup() %>%
  filter(week_specimen < max(week_specimen) | ndates == 7) %>%
  select(-ndates)

by_region <- english_pillars %>%
  filter(!is.na(age_group)) %>%
  select(-lower_age_limit, -positive, -total) %>%
  group_by(pillar, week_specimen, nhser_name, age_group) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  pivot_longer(starts_with("sgene"), names_to = "sgene_result",
               values_to = "n") %>%
  mutate(sgene_result = sub("^sgene_", "", sgene_result))

## age and region
p <- ggplot(mapping = aes(x = week_specimen, y = n, fill = sgene_result)) +
  facet_grid(age_group ~ nhser_name) +
  theme_cowplot() +
  scale_fill_brewer(palette = "Set1") +
  xlab("Week")

p_abs <- p +
  geom_bar(data = by_region %>% filter(pillar == "Pillar 2"),
           stat = "identity") +
  ylab("Cases")
p_prop <- p +
  geom_bar(data = by_region %>%
             filter(pillar == "Pillar 2", sgene_result != "n_a"),
           stat = "identity", position = "fill") +
  ylab("Proportion of cases")

ggsave(here::here("figure", "sgene.pdf"), p_abs, height = 20, width = 20)
ggsave(here::here("figure", "sgene_prop.pdf"), p_prop, height = 20, width = 20)

by_ltla <- english_pillars %>%
  select(-lower_age_limit, -positive, -total) %>%
  group_by(pillar, week_specimen, nhser_name, ltla_name, ltla_code) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  pivot_longer(starts_with("sgene"), names_to = "sgene_result",
               values_to = "n") %>%
  mutate(sgene_result = sub("^sgene_", "", sgene_result))

p <- ggplot(mapping = aes(x = week_specimen, y = n, fill = sgene_result)) +
  facet_wrap(. ~ ltla_name) +
  theme_cowplot() +
  scale_fill_brewer("S-Gene", palette = "Set1") +
  xlab("Week")

p_abs <- p +
  geom_bar(data = by_ltla %>% filter(pillar == "Pillar 2"),
           stat = "identity") +
  ylab("Cases")
p_prop <- p +
  geom_bar(data = by_ltla %>%
             filter(pillar == "Pillar 2", sgene_result != "n_a"),
           stat = "identity", position = "fill") +
  ylab("Proportion of cases")

ggsave(here::here("figure", "sgene_ltla.pdf"), p_abs, height = 20, width = 25)
ggsave(here::here("figure", "sgene_ltla_prop.pdf"), p_prop,
       height = 20, width = 25)

by_age <- english_pillars %>%
  select(-lower_age_limit, -positive, -total) %>%
  group_by(pillar, week_specimen, age_group) %>%
  summarise_if(is.numeric, sum) %>%
  pivot_longer(starts_with("sgene"), names_to = "sgene_result",
               values_to = "n") %>%
  mutate(sgene_result = sub("^sgene_", "", sgene_result)) %>%
  filter(!is.na(age_group))

p <- ggplot(by_age %>% filter(pillar == "Pillar 2"),
            aes(x = week_specimen, y = n, fill = sgene_result)) +
  facet_wrap(. ~ age_group) +
  theme_cowplot() +
  scale_fill_brewer("S-Gene", palette = "Set1") +
  xlab("Week")

p_abs <- p +
  geom_bar(data = by_age %>% filter(pillar == "Pillar 2"),
           stat = "identity") +
  ylab("Cases")
p_prop <- p +
  geom_bar(data = by_age %>%
             filter(pillar == "Pillar 2", sgene_result != "n_a"),
           stat = "identity", position = "fill") +
  ylab("Proportion of cases")

ggsave(here::here("figure", "sgene_age.pdf"), p_abs, height = 20, width = 20)
ggsave(here::here("figure", "sgene_age_prop.pdf"), p_prop,
       height = 20, width = 20)

by_ltla_aggregate <- by_ltla %>%
  filter(pillar == "Pillar 2") %>%
  mutate(week_infection = week_specimen - 7) %>%
  select(-pillar, -negative) %>%
  pivot_wider(names_from = sgene_result, values_from = n) %>%
  mutate(prop_variant = negative / (positive + negative),
         prop_variant_sd = sqrt(prop_variant * (1 - prop_variant) /
                                  (positive + negative)),
         samples = negative + positive,
         cases = n_a + negative + positive) %>%
  inner_join(rt_weekly %>% select(week_infection, ltla_name),
             by = c("week_infection", "ltla_name"))

saveRDS(by_ltla_aggregate, here("data", "processed", "sgene_by_ltla.rds"))

by_ltla_rt <- by_ltla_aggregate %>%
  inner_join(rt_weekly, by = c("week_infection", "ltla_name")) %>%
  select(week_infection, nhser_name, ltla_name, ltla_code, prop_variant,
         prop_variant_sd, samples, cases, rt_mean = mean, rt_sd = sd)

for (week in unique(as.character(by_ltla_rt$week_infection))) {
  p <- ggplot(by_ltla_rt %>%
              filter(week_infection == as.Date(week)),
              aes(x = prop_variant, y = rt_mean,
                              fill = nhser_name, size = cases)) +
    geom_jitter(pch = 21) +
    scale_fill_brewer("", palette = "Set1") +
    xlab("Proportion with S gene dropped") +
    ylab("Mean reproduction number") +
    theme_cowplot() +
    labs(size = "Cases since 9 October") +
    ggtitle(week)
  ggsave(here::here("figure", paste0("sgene_rt_", week, ".pdf")),
         width = 10, height = 7)
}

by_week_and_age <- english_pillars %>%
  filter(nhser_name %in% c("London", "South East", "East of England")) %>%
  select(-lower_age_limit, -positive, -total) %>%
  group_by(pillar, week_specimen, age_group) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  pivot_longer(starts_with("sgene"), names_to = "sgene_result",
               values_to = "n") %>%
  mutate(sgene_result = sub("^sgene_", "", sgene_result)) %>%
  filter(!is.na(age_group))

p <- ggplot(mapping = aes(x = age_group, y = n, fill = sgene_result)) +
  facet_wrap(. ~ week_specimen) +
  theme_cowplot() +
  scale_fill_brewer("S-Gene", palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Age group") +
  ylab("Proportion")

p_abs <- p +
  geom_bar(data = by_week_and_age %>% filter(pillar == "Pillar 2"),
           stat = "identity") +
  ylab("Cases")
p_prop <- p +
  geom_bar(data = by_week_and_age %>%
             filter(pillar == "Pillar 2", sgene_result != "n_a"),
           stat = "identity", position = "fill") +
  ylab("Proportion of cases")

ggsave(here::here("figure", "sgene_age_week.pdf"), p_abs,
       height = 10, width = 15)
ggsave(here::here("figure", "sgene_age_week_prop.pdf"), p_prop,
       height = 10, width = 15)

by_week_and_age_dist <- by_week_and_age %>%
  filter(pillar == "Pillar 2", sgene_result != "n_a") %>%
  group_by(week_specimen, sgene_result) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_dist <- ggplot(by_week_and_age_dist,
                 mapping = aes(x = age_group, y = prop, fill = sgene_result)) +
  facet_wrap(. ~ week_specimen) +
  theme_cowplot() +
  scale_fill_brewer("S-Gene", palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Age group") +
  ylab("Proportion (of those tested)") +
  geom_bar(stat = "identity", position = "dodge")

ggsave(here::here("figure", "sgene_age_week_dist.pdf"), p_dist,
       height = 10, width = 15)

by_age_dist <- by_week_and_age %>%
  filter(pillar == "Pillar 2", sgene_result != "n_a") %>%
  group_by(sgene_result, age_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(sgene_result) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_dist <- ggplot(by_age_dist,
                 mapping = aes(x = age_group, y = prop, fill = sgene_result)) +
  theme_cowplot() +
  scale_fill_brewer("S-Gene", palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Age group") +
  ylab("Proportion") +
  geom_bar(stat = "identity", position = "dodge")

ggsave(here::here("figure", "sgene_age_dist.pdf"), p_dist,
       height = 5, width = 7)
