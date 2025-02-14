#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# This is NOT needed if we don't use the getGenes function in the `.util/function-calculate-cell-type-signature.R` 
# Download R scripts for dittoSeq library
# This needs to be done only if the user wants to run `run-cell-types-annotation-gene-markers.R` (first time).
# cd util
# git clone https://github.com/dtm2451/dittoSeq.git

################################################################################################################
# Run module
# cd "$(dirname "${BASH_SOURCE[0]}")" 

#Rscript --vanilla run-cell-types-annotation-SingleR.R
Rscript --vanilla run-cell-types-annotation-gene-markers.R
#Rscript --vanilla 04-merge-cell-types-annotations-all.R
################################################################################################################
