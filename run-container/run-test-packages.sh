#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

########################################################################
# Read container path and file
containerdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")/..")
echo "$containerdir"

################################################################################################################
# Run module
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif \
            Rscript -e "rmarkdown::render('run-test-packages.Rmd', clean = FALSE,
                              output_dir = file.path('.'),
                              output_file = c(paste('Report-', 'run-test-packages', '-', Sys.Date(), sep = '')),
                              output_format = 'all')"