#' Get a combined data set of two epidemiological indicators by date of notification
#'
#' @param from "cases", "admissions" or "deaths"
#' @param to "cases", "admissions" or "deaths"
#' @param specimen_shift Numeric, number of weeks to shift the proportion of SGTF 
#' by. Defaults to 0 which implies specimen date. 
#' @param level Aggregation level of the data
#' @return A data frame containing primary and secondary notifications alongside
#' the proportion of cases that were SGTF positive.
#' @export
#' @importFrom dplyr select mutate filter group_by slice rename summarise
#' @importFrom tidyr inner_join drop_na
#' @importFrom vroom vroom
#' @importFrom lubridate floor_date weeks
#' @author Sam Abbott
get_notifications_data <- function(from = c("cases", "admissions", "deaths"),
                                   to = c("cases", "admissions", "deaths"),
                                   specimen_shift = 0, level = c("utla", "nhs region")){
  
  ## Arguments ---------------------------------------------------------------
  from <- match.arg(from, choices = c("cases", "admissions", "deaths"))
  to <- match.arg(to, choices = c("cases", "admissions", "deaths"))
  level <- match.arg(level, choices = c("utla", "nhs region"))
  
  ## Data --------------------------------------------------------------------
  base_url <-
    paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
           "master/subnational/united-kingdom-local/")
  data_path <- "/summary/reported_cases.csv"
  primary <- vroom(paste0(base_url, from, data_path)) %>% 
    select(loc = region, date, primary = confirm)
  
  secondary <- vroom(paste0(base_url, to, data_path)) %>% 
    select(loc = region, date, secondary = confirm)

  # link datasets 
  reports <- primary %>% 
    inner_join(secondary, by = c("loc", "date"))
  
  # Define covariates -------------------------------------------------------
  # get variant proportion
  sgene_by_utla <- readRDS(here("data", "sgene_by_utla.rds")) %>% 
    drop_na(prop_sgtf) %>% 
    mutate(week_specimen = week_infection + 7) %>% 
    filter(week_specimen > "2020-10-01") %>% 
    group_by(utla_name,week_specimen) %>% 
    # take only the first estimate of the proportion per UTLA (as brm_convolution
    # cannot currently handle repeats)
    slice() %>% 
    rename(loc = utla_name)
  
  # add week indicator
  reports <- reports %>% 
    mutate(week_specimen = floor_date(date, "week", week_start = wday(max(sgene_by_utla$week_specimen)) - 1))
  
  if (specimen_shift != 0) {
    reports %>% 
      mutate(week_specimen = week_specimen + weeks(specimen_shift))
  }
  
  # link with variant proportion
  secondary_with_cov <- reports %>% 
    inner_join(sgene_by_utla, by = c("week_specimen", "loc")) %>% 
    select(loc, date, week_specimen, region = nhser_name, 
           primary, secondary, prop_sgtf, samples)
  
  if (level %in% "nhs region") {
    secondary_with_cov <- secondary_with_cov %>% 
      mutate(positive = samples * prop_sgtf) %>% 
      group_by(date, week_specimen, region) %>% 
      summarise(
        primary = sum(primary),
        secondary = sum(secondary), 
        prop_sgtf = sum(positive) / sum(samples),
        samples = sum(samples),
        .groups = "drop"
      ) %>% 
      select(loc = region, date, week_specimen, primary, secondary, 
             prop_sgtf, samples)
  }
  
  # make normalised predictors
  secondary_with_cov <- secondary_with_cov %>% 
    mutate(normalised_primary = (primary - mean(primary)) / sd(primary),
           time = as.numeric(date),
           time = time - min(time),
           time = (time - mean(time)) / sd(time))
  
  return(secondary_with_cov)
}