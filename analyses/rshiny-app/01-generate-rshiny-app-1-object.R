#################################################################################
# This will run all scripts in the module
#################################################################################
# Load the Package with a Specific Library Path
# .libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.4")
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  library(reticulate)
  library(gridExtra)
  library(glue)
  library(R.utils)
  library(shiny)
  library(shinyhelper)
  library(DT)
  library(magrittr)
  library(ggdendro)
  library(ShinyCell)
  library(RColorBrewer)
  })

#################################################################################
# load config file
configFile <- paste0("../../project_parameters.Config.yaml")
if (!file.exists(configFile)){
  cat("\n Error: configuration file not found:", configFile)
  stop("Exit...")}

# read `yaml` file defining the `params` of the project and strategy analysis
yaml <- read_yaml(configFile)
#################################################################################
# Parameters
root_dir <- yaml$root_dir
integration_method <- yaml$integration_method
resolution_list <- yaml$resolution_list
PROJECT_NAME <- yaml$PROJECT_NAME
PI_NAME <- yaml$PI_NAME
condition_value <- yaml$condition_value
assay <- yaml$assay_filter_object
  
# Set up directories and paths to root_dir and analysis_dir
analysis_dir <- file.path(root_dir, "analyses") 
module_dir <- file.path(analysis_dir, "rshiny-app") 
annotation_results_dir <- file.path(analysis_dir, "cell-types-annotation", "results") 

# Input files
annotation_file <- file.path(annotation_results_dir, "seurat_obj_annotations.rds")

# Create results_dir
results_dir <- file.path(module_dir, "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)}

################################################################################
### library(ShinyCell) ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
################################################################################
# I modified the `makeShinyFiles.R` that comes from the library(ShinyCell) and rename as `makeShinyFiles_assay.R`
# The script was downloaded from here: https://github.com/SGDDNB/ShinyCell/blob/master/R
# The assay is by default always `RNA` in that script 
# and it doesn't fit the prerequisites for the  `sc-rna-seq-snap` pipeline as we might have a different name for the assays
source(paste0(module_dir, "/util/makeShinyFiles_assay.R"))

################################################################################################################
### Generate R shiny app ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
################################################################################################################

cat("Beginning to process results from", "annotation_file", "\n")
seu1 <- readRDS(annotation_file)
scConf1 <- createConfig(seu1)
makeShinyFiles_assay(seu1, scConf1, shiny.prefix = "sc1", shiny.dir = paste(results_dir, "shinyApp", sep = "/"))


cat("Make R shiny app for all files", "\n")
makeShinyCodesMulti(
  shiny.title = PROJECT_NAME, shiny.footnotes = PI_NAME,
  shiny.prefix = c("sc1"),
  shiny.headers = c("Annotations"),
  shiny.dir = paste(results_dir, "shinyApp", sep = "/")) 
################################################################################################################   
