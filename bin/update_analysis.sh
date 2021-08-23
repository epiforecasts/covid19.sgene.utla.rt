#!bin/bash

# update deps
Rscript -e "install.packages('devtools'); devtools::install_deps()"

# update Rt trajectories
Rscript R/extract_rt.r

# update data
Rscript R/extract_public_data.r

# combine data sources
Rscript R/combine_data.r

# fit models
Rscript R/fit_models.r

# compare moels
Rscript R/compare_models.r

# update report
bash bin/render_report.sh