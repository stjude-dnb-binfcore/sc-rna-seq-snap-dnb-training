#!/bin/bash

set -e
set -o pipefail

# Set up modules
module load R/4.4.0
module load pandoc/3.2
#module load pandoc/1.19.2.1


# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

################################################################################################################
# Run module
Rscript --vanilla 01-generate-rshiny-app.R
# Rscript --vanilla 02-launch-rshiny-app.R
