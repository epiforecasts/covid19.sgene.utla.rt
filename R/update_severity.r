
# updated lagged estimates
source(here::here("R", "lagged_severity.r"))
# update convolution estimates
source(here::here("R", "convolution_severity.r"))
# combine and produce results summaries
source(here::here("R", "combine_severity.r"))
# update the report
rmarkdown::render(here::here("severity-report.Rmd"))
