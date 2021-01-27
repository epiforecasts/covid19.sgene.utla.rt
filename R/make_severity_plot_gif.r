# Packages ----------------------------------------------------------------
library(here)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)
library(cowplot)
library(lemon)
library(gganimate) 
#devtools::install_github("thomasp85/transformr")
library(transformr)
library(gifski)

# Get plot data -----------------------------------------------------------
lagged_severity_data <- readRDS(here("data", "lagged_severity_data.rds"))

# Get plot function -------------------------------------------------------
source(here("R", "plot_severity.r"))
plot <- suppressMessages(plot_severity(lagged_severity_data, alpha = 0.6))


gif <- plot +
ggtitle('Date: {closest_state}') +
  transition_states(week_infection)

animate(gif, fps = 5, nframes = 100, renderer = gifski_renderer())

# Save --------------------------------------------------------------------
anim_save(here("output", "severity.gif"))