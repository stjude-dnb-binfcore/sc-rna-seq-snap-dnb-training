#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

########################################################################
# Read multiple values and assign them to variables by parsing yaml file
fastqc_dir=$(cat ../../project_parameters.Config.yaml | grep 'fastqc_dir:' | awk '{print $2}')
fastqc_dir=${fastqc_dir//\"/}  # Removes all double quotes

# Create results dir
mkdir results
mkdir results/01-fastqc-reports

################################################################################################################
# Run FastQC per library
fastqc -o results/01-fastqc-reports ${fastqc_dir}/*/*R2*.fastq.gz

################################################################################################################
# Run multiqc for all samples in the `results` dir
# to summarize results
 
cd results/01-fastqc-reports
multiqc . 

# rename folder
mv multiqc_data 02-multiqc-reports

# move the files related to multiqc in the main `results` dir
mv 02-multiqc-reports ../
mv multiqc_report.html ../
