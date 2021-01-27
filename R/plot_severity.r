library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)
library(cowplot)
library(lemon)

plot_severity <- function(df, alpha = 0.4) {
  
  df <- df %>% 
    mutate(data = map2(data, loc, function(df, loc) {
      if (is.null(df[["region"]])) {
        mutate(df, region =  loc)
      }else{
        df
      }}))
  
  prop_data <- df %>% 
    mutate(data = pmap(list(data, target, loc), 
                       function(df, y, z) {
                         df %>% 
                           mutate(rate = deaths / cases, type = y, loc = z)
                       })) %>% 
    pull(data) %>% 
    bind_rows() %>% 
    mutate(type = case_when(type %in% "cfr" ~ "Case fatality rate",
                            type %in% "chr" ~ "Case hospitalisation rate",
                            type %in% "hfr" ~ "Hospitalisation fatality rate")) %>% 
    mutate(loc = case_when(loc %in% "region" ~ "NHS region",
                           loc %in% "utla" ~ "Upper-tier local authority")) %>% 
    group_by() %>% 
    filter(rate < 1)
  
  suppressWarnings(ggplot(prop_data, aes(x = prop_sgtf, y = rate, size = samples,
                                         fill = region)) +
                     geom_jitter(pch = 21, alpha = alpha) +
                     scale_fill_brewer("", palette = "Set1") +
                     scale_y_continuous(labels = percent) +
                     xlab("Proportion with S gene dropped") +
                     ylab("Rate") +
                     facet_rep_grid(type ~ loc, scales = "free") +
                     theme_cowplot() +
                     labs(size = paste("Samples")) +
                     theme(legend.position = "bottom",
                           legend.direction = "horizontal",
                           legend.box = "vertical"))
}