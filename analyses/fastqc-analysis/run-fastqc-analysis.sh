#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Read root path
rootdir=$(realpath "./../..")
echo "$rootdir"

########################################################################
# Read metadata_dir from YAML configuration file
metadata_dir=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'metadata_dir:' | awk '{print $2}')
metadata_dir=${metadata_dir//\"/}  # Removes all double quotes
echo "Metadata directory: $metadata_dir"  # Output 

# Define the path to the metadata file (adjust to your actual file)
metadata_file="$metadata_dir/project_metadata.tsv"

# Check if metadata file exists
if [ ! -f "$metadata_file" ]; then
  echo "Error: Metadata file '$metadata_file' does not exist."
  exit 1
fi


########################################################################
# Define the column name for FASTQ (adjust to your actual header name)
fastq_column_name="FASTQ"

# Extract the column number of FASTQ from the header (first line)
column_number=$(head -n 1 "$metadata_file" | tr '\t' '\n' | grep -n "^$fastq_column_name$" | cut -d: -f1)

# Check if the FASTQ column exists
if [ -z "$column_number" ]; then
  echo "Error: Column '$fastq_column_name' not found in the metadata file."
  exit 1
fi

########################################################################
# Declare an array to store all FASTQ paths
declare -a fastqc_dir

# Read the FASTQ column values and split any comma-separated entries
while IFS= read -r fastq_value; do
  # Remove leading/trailing whitespace and carriage returns
  fastq_value=$(echo "$fastq_value" | tr -d '\r' | xargs)

  # Split by comma
  IFS=',' read -ra paths <<< "$fastq_value"
  
  for path in "${paths[@]}"; do
    # Clean individual path values too
    clean_path=$(echo "$path" | tr -d '\r' | xargs)
    [[ -n "$clean_path" ]] && fastqc_dir+=("$clean_path")
  done
done < <(awk -F'\t' -v col="$column_number" 'NR > 1 {print $col}' "$metadata_file")

  
# Check if the array has values
if [ "${#fastqc_dir[@]}" -gt 0 ]; then
  echo "fastqc_dir is an array with ${#fastqc_dir[@]} elements."
else
  echo "fastqc_dir is empty or failed to extract values."
fi

# Optionally print the FASTQ values for verification
echo "FASTQ values extracted:"
for dir in "${fastqc_dir[@]}"; do
  echo "$dir"
done

########################################################################
# Create results dir
mkdir -p results
mkdir -p results/01-fastqc-reports

########################################################################
################################################################################################################
###### STEP 1 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
################################################################################################################
# Run FastQC per library

# Loop through the array and echo each element
#for dir in "${fastqc_dir_array[@]}"; 
for dir in $(echo "${fastqc_dir[@]}" | tr ' ' '\n' | sort); 
  do
  echo "Processing directory: ${dir}"
  fastqc -o results/01-fastqc-reports ${dir}/*R2*.fastq.gz
done

################################################################################################################
###### STEP 2 ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
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

################################################################################################################
###### THE END ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
################################################################################################################
