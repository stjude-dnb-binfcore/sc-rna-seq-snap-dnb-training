# Pipeline for annotating cell types for sc-/sn-RNA-Seq Analysis in 10X Genomics data

## Usage

`run-cell-types-annotation.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`
- `future_globals_value` is hardwired coded in the `run-cell-types-annotation.R`. If necessary, user can increase/decrease resources.


### Run module on an interactive session on HPC within the container

To run all of the Rscripts in this module sequentially on an interactive session on HPC, please run the following command from an interactive compute node:

```
bash run-cell-types-annotation.sh
```

### Run module by using lsf on HPC within the container

There is also the option to run a lsf job on the HPC cluster by using the following command on an HPC node:

```
bsub < lsf-script.txt
```


## Folder content
This folder contains a script tasked to identify cell types annotation with SingleR by using a broad and fine resolution labeling across the project.

## Folder structure 

The structure of this folder is as follows:

```
├── 01-cell-types-annotation-SingleR-broad.Rmd
├── 02-cell-types-annotation-SingleR-fine.Rmd
├── lsf_script.txt
├── plots
|   ├── 01-cell-types-annotation-SingleR 
|       ├── 01_annotations_broad
|       └── 02_annotations_fine
|   ├── Report_cell_types_annotation_SingleR_broad_<Sys.Date()>.html
|   ├── Report_cell_types_annotation_SingleR_broad_<Sys.Date()>.pdf
|   ├── Report_cell_types_annotation_SingleR_fine_<Sys.Date()>.html
|   └── Report_cell_types_annotation_SingleR_fine_<Sys.Date()>.pdf
├── README.md
├── results
|   ├── 01-cell-types-annotation-SingleR 
|       ├── 01_annotations_broad
|       └── 02_annotations_fine
├── run-cell-types-annotation.R
├── run-cell-types-annotation.sh
└── util
|___└── function-cell-type-fractions.R
```

