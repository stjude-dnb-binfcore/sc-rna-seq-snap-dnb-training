#!/bin/bash

#BSUB -P project
#BSUB -J run-test-other-dependecies
#BSUB -q standard
#BSUB -n 1
#BSUB -R "rusage[mem=2500]"
#BSUB -cwd "."

########################################################################
# Set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")"

########################################################################
# Read container path and file
containerdir=$(realpath "$(dirname "${BASH_SOURCE[0]}")/..")
echo "$containerdir"


singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif python3 --version
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif cellranger --version
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif fastqc --version
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif multiqc --version
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif pandoc --version
singularity exec ${containerdir}/rstudio_4.4.0_seurat_4.4.0_latest.sif tex --version
