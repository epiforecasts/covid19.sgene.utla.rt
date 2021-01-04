
# Local area reproduction numbers and S-gene dropouts

This repository contains the data and code for our report exploring the association between lower-tier local area (LTLA) reproduction number estimates in England and the proportion of Covid-19 tests negative for the S-gene. The report can be found [here](https://github.com/epiforecasts/covid19.sgene.ltla.rt/report.pdf/).

## Reproducibility

## Data

All data used in the analysis can be found in the `data` folder in `rds` format. Available data include: 

- `ltla_rt_with_covariates.rds`: LTLA level weekly reproduction number estimates combined with estimates of the proportion of tests that were S-gene negative, normalised Google mobility data, and tier status by local authority over time.
- `rt.rds`: Summarised daily LTLA reproduction number estimates using both a short and a long generation time.
- `sgene_by_ltla.rds`: Weekly test positivity data for the S-gene by LTLA.
- `mobility.rds`: Normalised Google mobility data stratified by context. 
- `tiers.rds`: LTLA level tier level over time.

## Dependencies

The dependencies for this analysis can be installed using (in the working directory of the analysis):

```r
install.packages("devtools")
devtools::install_deps()
```

## Code


All publicly available data sources can be re-extracted using (here and following in the working directory of the analysis):

```r
source(here::here("R/extract_public_data.r"))
```

All data sources can then be combined into the analysis dataset using:

```r
source(here::here("R/combine_data.r"))
```

The statistical models considered can be refit using,

```r
source(here::here("R/fit_models.r"))
```

and compared with

```r
source(here::here("R/compare_models.r"))
```

The report can be regenerated (once the models have been refit and compared) using:

```r
rmarkdown::render("report.Rmd")
```

Alternatively all steps can be reproduced using the following bash script: 

```bash 
bash bin/update_analysis.sh
```
## Installation

You can install the released version of covid19.sgene.ltla.rt from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("covid19.sgene.ltla.rt")
```

