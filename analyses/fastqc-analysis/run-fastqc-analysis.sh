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
# Extract column numbers
sample_column="ID"
fastq_column="FASTQ"

sample_col_num=$(head -n 1 "$metadata_file" | tr '\t' '\n' | grep -n "^$sample_column$" | cut -d: -f1)
fastq_col_num=$(head -n 1 "$metadata_file" | tr '\t' '\n' | grep -n "^$fastq_column$" | cut -d: -f1)

echo "Sample column: $sample_col_num, FASTQ column: $fastq_col_num"

mkdir -p results/01-fastqc-reports

# Read each row and process all FASTQ paths per sample
tail -n +2 "$metadata_file" | while IFS=$'\t' read -r -a fields; do
  sample="${fields[$((sample_col_num - 1))]}"
  fastq_field="${fields[$((fastq_col_num - 1))]}"

  # Split on commas to handle technical replicates
  IFS=',' read -ra paths <<< "$fastq_field"

  rep=1
  for raw_path in "${paths[@]}"; do
    clean_path=$(echo "$raw_path" | tr -d '\r' | sed 's/^["'\'']//; s/["'\'']$//')

    echo "Processing sample: $sample, replicate: $rep"

    for file in "$clean_path"/*R2*.fastq.gz; do
      echo "  Running FastQC on: $file"
      base=$(basename "$file" .fastq.gz)

      # Generate unique temp filename (used by FastQC for sample name)
      unique_name="${sample}_rep${rep}"
      temp_fastq="${base}_${unique_name}.fastq.gz"

      # Create symbolic link to avoid modifying original file
      ln -s "$file" "$temp_fastq"

      # Run FastQC on the symlink to embed the correct sample name
      fastqc -o results/01-fastqc-reports "$temp_fastq" --threads 8

      # Clean up the symlink
      rm "$temp_fastq"
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
