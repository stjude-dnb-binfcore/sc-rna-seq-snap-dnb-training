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

metadata_file=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'metadata_file:' | awk '{print $2}')
metadata_file=${metadata_file//\"/}  # Removes all double quotes
echo "Metadata file: $metadata_file"  # Output 

# Define the path to the metadata file (adjust to your actual file)
metadata_file="$metadata_dir/$metadata_file"

# Check if metadata file exists
if [ ! -f "$metadata_file" ]; then
  echo "Error: Metadata file '$metadata_file' does not exist."
  exit 1
fi


################################################################################################################
# Extract sample names and fastq paths in parallel
sample_column="ID"
fastq_column="FASTQ"

sample_col_num=$(head -n 1 "$metadata_file" | tr '\t' '\n' | grep -n "^$sample_column$" | cut -d: -f1)
fastq_col_num=$(head -n 1 "$metadata_file" | tr '\t' '\n' | grep -n "^$fastq_column$" | cut -d: -f1)

echo "Sample column: $sample_col_num, FASTQ column: $fastq_col_num"

# Create results dir
mkdir -p results/01-fastqc-reports

# Read metadata file line by line (excluding header)
tail -n +2 "$metadata_file" | while IFS=$'\t' read -r -a fields; do
  sample="${fields[$((sample_col_num - 1))]}"
  fastq_field="${fields[$((fastq_col_num - 1))]}"

  # Split fastq_field by comma
  IFS=',' read -ra paths <<< "$fastq_field"

  rep=1
  for path in "${paths[@]}"; do
    clean_path=$(echo "$path" | tr -d '\r' | xargs)

    echo "Processing sample: $sample, replicate: $rep"
    for file in "$clean_path"/*R2*.fastq.gz; do
      echo "  Running FastQC on: $file"

      base=$(basename "$file" .fastq.gz)
      
      # Run FastQc
      # FastQC doesn't scale linearly with thread count. After about 4–8 threads, the gains start to taper off.
      fastqc -o results/01-fastqc-reports "$file" --threads 8
      
      # Rename files to avoid overwriting if there are multiple techincal replicates
      mv "results/01-fastqc-reports/${base}_fastqc.html" "results/01-fastqc-reports/${base}_rep${rep}_fastqc.html"
      mv "results/01-fastqc-reports/${base}_fastqc.zip" "results/01-fastqc-reports/${base}_rep${rep}_fastqc.zip"
    done

    ((rep++))
  done
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
