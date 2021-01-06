#!bin/bash

# update deps
Rscript -e "install.packages('devtools'); devtools::install_deps()"

# update data
Rscript R/extract_public_data.r

# combine data sources
Rscript R/combine_data.r

# fit models
Rscript R/fit_models.r

# compare moels
Rscript R/compare_models.r

# update report
Rscript -e "rmarkdown::render('report.Rmd')"