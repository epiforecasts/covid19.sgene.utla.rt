
# Estimating the increase in reproduction number associated with the Delta variant using local area dynamics in England

**This repository is under active development. See the [releases](https://github.com/epiforecasts/covid19.sgene.utla.rt/releases) for stable analyses.**

This repository contains the data and code for analyses exploring the association between upper-tier local area (UTLA) reproduction number estimates in England and the proportion of COVID-19 tests negative/positive for the S-gene as a proxy for the Alpha/Delta variants. A preprint describing the work in more detail can be found on [medRxiv](https://www.medrxiv.org/content/10.1101/2021.11.30.21267056v1.full).

## Reproducibility

### Data

All data used in the analysis can be found in the `data` folder in `rds` format. Available data include: 

- `utla_rt_with_covariates.rds`: UTLA level weekly reproduction number estimates combined with estimates of the proportion of tests that were S-gene negative, normalised Google mobility data, and tier status by local authority over time.
- `rt_weekly.rds`: Summarised weekly UTLA reproduction number estimates using both a short and a long generation time.
- `utla_cases.rds`: UTLA level COVID-19 test positive cases.
- `sgene_by_utla.rds`: Weekly test positivity data for the S-gene by UTLA.
- `mobility.rds`: Normalised Google mobility data stratified by context. 
- `tiers.rds`: UTLA level tier level over time.

### Dependencies

The dependencies for this analysis can be installed using (in the working directory of the analysis):

```r
install.packages("devtools")
devtools::install_deps()
```

### Code

Rt estimates from the [EpiForecasts](http://epiforecasts.io/covid) web site can be updated using (here and following in the working directory of the analysis):

```r
source(here::here("R/extract_rt.r"))
```

All publicly available covariates can be re-extracted using:

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

### Archived results

Instead of rerunning the code from scratch archived results can be retrieved from the [OSF](https://osf.io/h6e8n/) which also contains an archive of the code used for the analysis. Note that this will overwrite any files saved in the `output` folder.

```r
source(here::here("R/download_output.r"))
```
