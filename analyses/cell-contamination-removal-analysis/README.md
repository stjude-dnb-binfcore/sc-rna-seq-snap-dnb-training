# Pipeline for removing cell contamination, and repeating normalization and integration steps of analysis for sc-/sn-RNA-Seq Analysis in 10X Genomics data

## Usage

`run-cell-contamination-removal-analysis.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`.
- `run-cell-contamination-removal-analysis.sh`: comment in/out according to which step user wants to run, i.e., `run-cell-contamination-removal-analysis-steps-1-2-3.R` or `run-cell-contamination-removal-analysis-steps-4.R`.
- `future_globals_value` is hardwired coded in the `run-cell-contamination-removal-analysis-steps-1-2-3.R` and `run-cell-contamination-removal-analysis-steps-4.R`. If necessary, user can increase/decrease resources.


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-cell-contamination-removal-analysis.sh
```

### Run module by using lsf on HPC within the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```


## Folder content
This folder contains a script tasked to remove cells and repeat normalization and integration steps analysis. This object can be used for further analysis, e.g., for the `cell-types-annotation` module.

## Analysis strategy:

After inspection of the two rounds of results in the `cluster-cell-calling` module, user can remove clusters if necessary. This is relevant to projects in which we need to remove cell contamination from the object, e.g., clusters. This is the case, e.g., in PDX projects, there might be both human and mouse clusters identified after the `02-find-markers.Rmd` step of the the `cluster-cell-calling` module. In this case, we reccomened the user to run the `cell-contamination-removal-analysis` module that allows to remove clusters, repeat normalization and integration steps. This object can then be used for cell type annotation or other type of analysis.


## Folder structure 

The structure of this folder is as follows:

```
├── 01-cell-contamination-removal.Rmd
├── 02-integrative-analysis.Rmd
├── 03-cluster-cell-calling.Rmd
├── 04-find-markers.Rmd
├── lsf-script.txt
├── plots
|   ├── 01_cell_contamination_removal
|   ├── 02_integration_{integration_method}
|   ├── 03_cluster_cell_calling_{resolution}
|   ├── Report_cell_contamination_removal_<Sys.Date()>.html
|   ├── Report_integration_{integration_method}_<Sys.Date()>.html
|   ├── Report_cluster_cell_calling_{resolution}_<Sys.Date()>.html
|   ├── Report_cell_contamination_removal_<Sys.Date()>.pdf
|   ├── Report_integration_{integration_method}_<Sys.Date()>.pdf
|   └── Report_cluster_cell_calling_{resolution}_<Sys.Date()>.pdf
├── README.md
├── results
|   ├── 01_cell_contamination_removal
|   ├── 02_integration_{integration_method}
|   └── 03_cluster_cell_calling_{resolution}
├── run-cell-contamination-removal-analysis-steps-1-2-3.R
├── run-cell-contamination-removal-analysis-steps-4.R
└── run-cell-contamination-removal-analysis.sh
```

