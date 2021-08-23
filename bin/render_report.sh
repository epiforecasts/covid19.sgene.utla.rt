#!bin/bash

# render pdf report
Rscript -e "rmarkdown::render('report.Rmd')"

# render doc for WOR
Rscript -e "rmarkdown::render('report.Rmd', output_format = 'word_document')"