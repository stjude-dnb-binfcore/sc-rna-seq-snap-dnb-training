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
integration_results_dir <- file.path(analysis_dir, "integrative-analysis", "results") 
clustering_results_dir <- file.path(analysis_dir, "cluster-cell-calling", "results", "01_cluster_cell_calling_custom_multiple") 

annotation_results_dir <- file.path(analysis_dir, "cell-types-annotation", "results") 
broad_SingleR_results_dir <- file.path(annotation_results_dir, "01_cell_types_annotation_SingleR", "01_annotations_broad") 
fine_SingleR_results_dir <- file.path(annotation_results_dir, "01_cell_types_annotation_SingleR", "02_annotations_fine") 

# Input files
integration_file <- file.path(integration_results_dir, glue::glue("seurat_obj_integrated_{integration_method}.rds"))
#clustering_file <- file.path(clustering_results_dir, glue::glue("seurat_obj_integrated_{integration_method}_clusters_{resolution_list}.rds"))
clustering_file <- file.path(clustering_results_dir, glue::glue("seurat_obj_integrated_{integration_method}_clusters_all.rds"))
broad_SingleR_file <- file.path(broad_SingleR_results_dir, "seurat_obj_SingleR_broad.rds")
fine_SingleR_file <- file.path(fine_SingleR_results_dir, "seurat_obj_SingleR_fine.rds")

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

cat("Beginning to process results from", "integration_file", "\n")
seu1 <- readRDS(integration_file)
scConf1 <- createConfig(seu1)
scConf1 <- modColours(scConf1, meta.to.mod = condition_value, new.colours= c("red", "black"))
makeShinyFiles_assay(seu1, scConf1, shiny.prefix = "sc1", shiny.dir = paste(results_dir, "shinyApp", sep = "/"))


cat("Beginning to process results from", "clustering_file", "\n")
seu2 <- readRDS(clustering_file)
scConf2 <- createConfig(seu2)
scConf2 <- modColours(scConf2, meta.to.mod = condition_value, new.colours= c("red", "black"))
makeShinyFiles_assay(seu2, scConf2, shiny.prefix = "sc2", shiny.dir = paste(results_dir, "shinyApp", sep = "/"))


cat("Beginning to process results from", "broad_SingleR_file", "\n")
seu3 <- readRDS(broad_SingleR_file)
scConf3 <- createConfig(seu3)
scConf3 <- modColours(scConf3, meta.to.mod = condition_value, new.colours= c("red", "black"))
makeShinyFiles_assay(seu3, scConf3, shiny.prefix = "sc3", shiny.dir = paste(results_dir, "shinyApp", sep = "/"))


cat("Beginning to process results from", "fine_SingleR_file", "\n")
seu4 <- readRDS(fine_SingleR_file)
scConf4 <- createConfig(seu4)
scConf4 <- modColours(scConf4, meta.to.mod = condition_value, new.colours= c("red", "black"))
makeShinyFiles_assay(seu4, scConf4, shiny.prefix = "sc4", shiny.dir = paste(results_dir, "shinyApp", sep = "/"))


cat("Make R shiny app for all files", "\n")
makeShinyCodesMulti(
  #shiny.title = "snRNA-Seq of DS in Df1 Mice", shiny.footnotes = "Mary Patton's project",
  shiny.title = PROJECT_NAME, shiny.footnotes = PI_NAME,
  shiny.prefix = c("sc1", "sc2", "sc3", "sc4"),
  shiny.headers = c("Integration", "Cluster cell calling", "broad_SingleR", "fine_SingleR"),
  shiny.dir = paste(results_dir, "shinyApp", sep = "/")) 
################################################################################################################   
