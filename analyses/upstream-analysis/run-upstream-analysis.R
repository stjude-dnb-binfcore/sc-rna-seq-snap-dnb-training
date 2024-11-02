#################################################################################
# This will run all scripts in the `upstream-analysis`
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)})

#################################################################################
# load config file
configFile <- paste0("../../project_parameters.Config.yaml")
if (!file.exists(configFile)){
  cat("\n Error: configuration file not found:", configFile)
  stop("Exit...")}

# read `yaml` file defining the `params` of the project and strategy analysis
yaml <- read_yaml(configFile)

#################################################################################
# Set up directories and paths to root_dir and analysis_dir
root_dir <- yaml$root_dir
analysis_dir <- file.path(root_dir, "analyses", "upstream-analysis") 

################################################################################################################
# Run Rmd scripts to process data per method
################################################################################################################
# (1) Seurat QC metrics
# Run the seurat_qc script for each sample/library and save html/pdf reports per each
source(paste0(analysis_dir, "/", "01B_run_seurat_qc_multiple_samples.R"))

################################################################################################################
# (2) Estimating and filtering out ambient mRNA (`empty droplets`)

################################################################################################################
# (3) Estimating and filtering out doublets

################################################################################################################
# (4) Merging filtered data


################################################################################################################
# (5) Final QC summary report

################################################################################################################   

