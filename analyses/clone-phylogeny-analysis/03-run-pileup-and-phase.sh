#!/bin/bash
# 
#BSUB -P project
#BSUB -J 02-create-pileup-and-phase
#BSUB -oo 02-create-pileup-and-phase-job.out -eo 02-create-pileup-and-phase-job.err
#BSUB -n 24
#BSUB -R "rusage[mem=128GB] span[hosts=1]"
#BSUB -cwd "."

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

#######################################################
# Read multiple values and assign them to variables by parsing yaml file
root_dir=$(cat ../../project_parameters.Config.yaml | grep 'root_dir:' | awk '{print $2}')
root_dir=${root_dir//\"/}  # Removes all double quotes
echo "${root_dir}"  # Output: This is a string with quotes.

module_dir=${root_dir}/analyses/clone-phylogeny-analysis
echo "${module_dir}"

cellranger_data_dir=$(cat ../../project_parameters.Config.yaml | grep 'data_dir:' | awk '{print $2}')
cellranger_data_dir=${cellranger_data_dir//\"/}  # Removes all double quotes
echo "$cellranger_data_dir"  # Output: This is a string with quotes.
########################################################################
# If your `project_metadata` is not in `*.txt` file format
# use the following code line to convert it
# cat "${metadata_dir}"/project_metadata.tsv | sed 's/,/\t/g' > ./input/project_metadata.txt

########################################################################
# Directories and paths to file Inputs/Outputs
mkdir -p ./results/02-create-pileup-and-phase
mkdir -p ${module_dir}/input/sample_barcode/

# Define array with samples
# sample=("SJDSRCT049192_D1" "SJDSRCT049192_R1" "SJRHB063823_D1" "SJRHB030680_R1_sc")
mapfile -t sample < <(grep '^ *-' ../../project_parameters.Config.yaml | sed 's/^ *- *//')

##################
# Check if the variable is an array by testing if it has elements
if [ "${#sample[@]}" -gt 0 ]; then
  echo "sample is an array with ${#sample[@]} elements."
else
  echo "sample is not an array or is empty."
fi
echo "${sample}"  # Output: This is a string with quotes.

# Output the array elements
echo "${sample[@]}"  # Output: sample1 sample2 etc

# Run loop for each sample
for i in "${!sample[@]}"; do
  
  echo "Sample: ${sample[i]}. Processing..."
  
  data_dir=${cellranger_data_dir}/${sample[i]}
  echo "$data_dir"
  
  mkdir -p ${module_dir}/results/02-create-pileup-and-phase/${sample[i]}
  mkdir -p ${module_dir}/input/sample_barcode/${sample[i]}
  
  # Step 2
  # To copy and unzip the barcodes file from the alignment folder for each sample
  cd ..
    
  cp ${data_dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ${module_dir}/input/sample_barcode/${sample[i]}
  gunzip ${module_dir}/input/sample_barcode/${sample[i]}/barcodes.tsv.gz
  
  # Run pileup_and_phase.R (installed in the container)
  Rscript --vanilla /numbat/inst/bin/pileup_and_phase.R \
          --label ${sample[i]} \
          --samples ${sample[i]} \
          --bams ${data_dir}/outs/possorted_genome_bam.bam \
          --barcodes ${module_dir}/input/sample_barcode/${sample[i]}/barcodes.tsv \
          --outdir ${module_dir}/results/02-create-pileup-and-phase/${sample[i]} \
          --gmap ${module_dir}/input/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
          --snpvcf ${module_dir}/input/references/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
          --paneldir ${module_dir}/input/references/1000G_hg38 \
          --ncores 64 \
          --UMItag Auto
            
  echo "Sample: ${sample[i]}. Processing completed succefully."
done
