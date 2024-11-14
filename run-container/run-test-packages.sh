#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

################################################################################################################
# Run module
Rscript -e "rmarkdown::render('run-test-packages.Rmd', clean = FALSE,
                              output_dir = file.path('.'),
                              output_file = c(paste('Report-', 'run-test-packages', '-', Sys.Date(), sep = '')),
                              output_format = 'all')"