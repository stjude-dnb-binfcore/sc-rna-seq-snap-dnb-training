#######################################################################################################
# https://broadinstitute.github.io/KrumlovSingleCellWorkshop2020/data-wrangling-scrnaseq-1.html
# https://www.rdocumentation.org/packages/Seurat/versions/3.1.1/topics/GetAssayData
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(yaml)
  library(tidyverse)
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
# Set up directories and paths to file Inputs/Outputs
root_dir <- yaml$root_dir
#metadata_dir <- yaml$metadata_dir
assay <- yaml$assay_annotation_module
input_dir <- yaml$input_dir_clone_phylogeny_module
analysis_dir <- file.path(root_dir, "analyses", "clone-phylogeny-analysis") 
#input_dir <- file.path(analysis_dir, "input") 
#input_dir <- file.path(analysis_dir, "input", "dyergrp_projects_ALSF_Pediatric_Atlas") 

# List all directories recursively
#all_dirs <- list.dirs(input_dir, recursive = TRUE)

# Filter directories that end with "output"
#sample_input_dir <- all_dirs[grep("/output$", all_dirs)]  # Ensure the directory ends with "output"
#print(sample_input_dir)

# Input files
#seurat_obj_file <- c(dir(path = sample_input_dir, pattern =  "_Seurat.Rds", full.names = TRUE, recursive = TRUE))
#print(dir(input_dir, full.names = TRUE))
seurat_obj_file <- dir(path = input_dir, pattern =  "seurat_obj.*\\.rds", full.names = TRUE, recursive = TRUE)
print(seurat_obj_file)


# Create module_results_dir
module_results_dir <- 
  file.path(analysis_dir, paste0("results"))
if (!dir.exists(module_results_dir)) {
  dir.create(module_results_dir)
}

# Create step_results_dir
step_results_dir <- 
  file.path(module_results_dir, paste0("01-create-count-mat"))
if (!dir.exists(step_results_dir)) {
  dir.create(step_results_dir)
}


#######################################################
# Read metadata file and define `sample_name`
#metadata_file <- file.path(metadata_dir, "project_metadata.tsv") # metadata input file

# Read metadata file and define `sample_name`
#project_metadata <- read.csv(metadata_file, sep = "\t", header = TRUE)
#sample_name <- unique(as.character(project_metadata$ID))
#sample_name <- sort(sample_name, decreasing = FALSE)
#print(sample_name)


#######################################################
# Create list 
seurat_obj_list <- list()
count_mat_list <- list()
sample_name <- c()

for (i in seq_along(seurat_obj_file)) {
  
  # Extract sample name from the file path
  sample_name <- c(sample_name, gsub("Create-", "", str_split_fixed(seurat_obj_file[i], "/", 16)[, 15]))
  print(sample_name[i])  # Print the current sample name for debugging
}

# Sort the sample names outside the loop if needed
sample_name <- sort(sample_name, decreasing = FALSE)
print(sample_name)  # Print all sample names to check if they are stored correctly

for (i in seq_along(sample_name)){

  # Create results_dir per sample
  results_dir <- file.path(step_results_dir, sample_name[i])
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)}
  
  # Read data per sample
  cat("Processing for sample:", sample_name[i], "\n")
  
  # Read the Seurat object
  seurat_obj_list[[i]] <- readRDS(seurat_obj_file[i]) # Use double brackets to assign to the list
  # Save the Seurat object to the results directory
  saveRDS(seurat_obj_list[[i]], file = file.path(results_dir, paste0(sample_name[i], "seurat_obj.rds")))
  cat("Saved Seurat object for sample:", sample_name[i], "\n")
  
  # Extract count matrix
  count_mat_list[[i]] <- as.matrix(seurat_obj_list[[i]]@assays[[assay]]@counts)

  # Save the count matrix
  saveRDS(count_mat_list[[i]], file = paste0(results_dir, "/", "count_mat.rds"))
  
  cat("Number of cells in count_mat_list[[i]]:", ncol(count_mat_list[[i]]), "\n")
  cat("Number of genes in count_mat_list[[i]]:", nrow(count_mat_list[[i]]), "\n")
  
  cat("Processing complete for sample:", sample_name[i], "\n")
  
}


