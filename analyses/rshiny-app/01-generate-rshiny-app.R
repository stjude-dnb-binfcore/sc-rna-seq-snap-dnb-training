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
PROJECT_NAME <- yaml$PROJECT_NAME
PI_NAME <- yaml$PI_NAME
condition_value <- yaml$condition_value
assay <- yaml$assay_filter_object
  
# Set up directories and paths to root_dir and analysis_dir
analysis_dir <- file.path(root_dir, "analyses") 
module_dir <- file.path(analysis_dir, "rshiny-app") 

annotation_results_dir <- file.path(analysis_dir, "cell-types-annotation", "results") 
#broad_SingleR_results_dir <- file.path(annotation_results_dir, "01_cell_types_annotation_SingleR", "01_annotations_broad") 
#fine_SingleR_results_dir <- file.path(annotation_results_dir, "01_cell_types_annotation_SingleR", "02_annotations_fine") 
annotations_all_results_dir <- file.path(annotation_results_dir, "04_cell_types_annotations_all") 

# Input files
#broad_SingleR_file <- file.path(broad_SingleR_results_dir, "seurat_obj_SingleR_broad.rds")
#fine_SingleR_file <- file.path(fine_SingleR_results_dir, "seurat_obj_SingleR_fine.rds")
annotations_all_file <- file.path(annotations_all_results_dir, "seurat_obj_cell_types_annotations_all.rds")

# Create results_dir
results_dir <- file.path(module_dir, "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)}


################################################################################
### library(ShinyCell) ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
################################################################################
# It was downloaded from here: https://github.com/SGDDNB/ShinyCell/blob/master/R
# We modified the function to account for any type of assay
source(paste0(module_dir, "/util/makeShinyFiles_assay.R"))

################################################################################################################
### Generate R shiny app ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
################################################################################################################
#cat("Beginning to process results from", "broad_SingleR_file", "\n")
#seu1 <- readRDS(broad_SingleR_file)
#scConf1 <- createConfig(seu1)
#makeShinyFiles_assay(seu1, scConf1, shiny.prefix = "sc1", shiny.dir = paste(results_dir, "shinyApp", sep = "/"))


#cat("Beginning to process results from", "fine_SingleR_file", "\n")
#seu2 <- readRDS(fine_SingleR_file)
#scConf2 <- createConfig(seu2)
#makeShinyFiles_assay(seu2, scConf2, shiny.prefix = "sc2", shiny.dir = paste(results_dir, "shinyApp", sep = "/"))

cat("Beginning to process results from", "annotations_all_file", "\n")
seu1 <- readRDS(annotations_all_file)
scConf1 <- createConfig(seu1)
makeShinyFiles_assay(seu1, scConf1, shiny.prefix = "sc1", shiny.dir = paste(results_dir, "shinyApp", sep = "/"))

cat("Make R shiny app for all files", "\n")
makeShinyCodesMulti(
  shiny.title = PROJECT_NAME, shiny.footnotes = PI_NAME,
  shiny.prefix = c("sc1"),
  shiny.headers = c("annotations_all"),
  #shiny.prefix = c("sc1", "sc2"),
  #shiny.headers = c("broad_SingleR", "fine_SingleR"),
  shiny.dir = paste(results_dir, "shinyApp", sep = "/")) 
################################################################################################################   
