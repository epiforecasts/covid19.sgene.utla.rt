##' Get a combined data set of two epidemiological indicators by date of infection
##'
##' Indicators will be matched to the time of the "from" indicator, with a lag
##' given by \code{infections_lag}
##' @param from "cases", "admissions" or "deaths"
##' @param to "cases", "admissions" or "deaths"
##' @param type "backcalculated" (via the EpiNow2 model) or "lagged" (fixed lag,
##' estimated from data via correlation analysis).
##' @param level Aggregation level of the data
##' @return a data frame of the two indicators matched to the same time
##' @author Sebastian Funk
get_infections_data <- function(from = c("cases", "admissions", "deaths"),
                     to = c("cases", "admissions", "deaths"),
                     type = c("lagged", "backcalculated"),
                     level = c("utla", "nhs region"),
                     infections_lag = 7) {
  
  ## Arguments ---------------------------------------------------------------
  from <- match.arg(from, choices = c("cases", "admissions", "deaths"))
  to <- match.arg(to, choices = c("cases", "admissions", "deaths"))
  type <- match.arg(type, choices = c("lagged", "backcalculated"))
  
  ## Data --------------------------------------------------------------------
  base_url <-
    paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
           "master/subnational/united-kingdom-local")
  
  if (type == "backcalculated") {
    file_name <- "cases_by_infection.csv"
  } else if (type == "lagged") {
    file_name <- "reported_cases.csv"
  }
  
  infs <- lapply(c(from, to), function(x) {
    vroom(paste(base_url, x, "summary", file_name, sep = "/")) %>%
      mutate(data = x) %>%
      mutate(date = date - infections_lag)
  })
  
  ## if type is lagged, add estimate of the fixed lag
  if (type == "lagged") {
    all_dates <-
      as.Date(unique(unlist(lapply(infs, function(x) as.character(x$date)))))
    combined <- tibble(date = all_dates) %>%
      full_join(infs[[1]], by = "date") %>%
      full_join(infs[[2]], by = c("date", "region")) %>%
      filter(!is.na(confirm.x), !is.na(confirm.y)) %>%
      group_by(date, region) %>%
      summarise(cases = sum(confirm.x),
                deaths = sum(confirm.y), .groups = "drop") %>%
      ungroup()
    
    ## find lag that maximises mean correlation across UTLAs
    max <- 28
    cor <- vapply(seq(0, max), function(lag) {
      combined %>%
        group_by(region) %>%
        mutate(deaths = lead(deaths, n = lag)) %>%
        filter(!is.na(deaths)) %>%
        mutate(id = 1:n()) %>%
        filter(id <= (max(id) - max + lag)) %>%
        summarise(cor = cor(cases, deaths), .groups = "drop") %>%
        filter(!is.na(cor)) %>%
        summarise(mean = mean(cor)) %>%
        .$mean
    }, 0)
    
    set_lag <- which.max(cor) - 1
    infs[[2]] <- infs[[2]] %>%
      mutate(date = date - set_lag)
    
    infs <- lapply(infs, function(x) {
      x %>%
        rename(median = confirm) %>%
        mutate(type = "estimate")
    })
  } else {
    set_lag <- NULL
  }
  
  ## link datasets and pivot wider
  infections <- bind_rows(infs) %>%
    filter(type == "estimate") %>%
    select(utla_name = region, date, data, value = median) %>%
    pivot_wider(names_from = "data") %>%
    group_by(utla_name) %>%
    complete(date = seq(min(date), max(date), by = "day")) %>%
    drop_na()
  
  ## Define covariates -------------------------------------------------------
  ## get variant proportion
  sgene_by_utla <- readRDS(here("data", "sgene_by_utla.rds")) %>%
    drop_na(prop_sgtf) %>%
    filter(week_infection > "2020-10-01")
  
  week_start <- lubridate::wday(max(sgene_by_utla$week_infection)) - 1
  ## make infections weekly summary
  weekly_infections <- infections %>%
    mutate(week_infection =
             floor_date(date, "week", week_start = week_start)) %>%
    select(-date) %>%
    group_by(utla_name, week_infection) %>%
    add_count() %>%
    group_by(utla_name, week_infection, n) %>%
    summarise_all(sum) %>%
    ungroup() %>%
    filter(n == 7) %>%
    select(-n)
  
  ## link with variant proportion
  deaths_with_cov <- weekly_infections %>%
    inner_join(sgene_by_utla, by = c("week_infection", "utla_name")) %>%
    rename(utla = utla_name, region = nhser_name) %>%
    select(-utla_code)
  
  ## factorise time
  deaths_with_cov <- deaths_with_cov %>%
    mutate(time = factor((week_infection - min(week_infection)) / 7))
  
  ## rename everything to cases/deaths to make things easier later
  deaths_with_cov <- deaths_with_cov %>%
    rename(!!!c(cases = from, deaths = to)) %>% 
    rename(loc = utla)
  
  if (level %in% "nhs region") {
    deaths_with_cov <- deaths_with_cov %>% 
      mutate(positive = samples * prop_sgtf) %>% 
      group_by(week_infection, region) %>% 
      summarise(
        cases = sum(cases),
        deaths = sum(deaths), 
        prop_sgtf = sum(positive) / sum(samples),
        samples = sum(samples),
        .groups = "drop"
      ) %>% 
      select(loc = region, week_infection, cases, deaths, prop_sgtf, samples)
  }
  
  # add normalised cases as a predictor
  deaths_with_cov <- deaths_with_cov %>% 
    mutate(normalised_cases = (cases - mean(cases)) / sd(cases))
  
  if (!is.null(set_lag)) {
    deaths_with_cov <- deaths_with_cov %>%
      mutate(lag = set_lag)
  }
  return(deaths_with_cov)
}