# Pipeline for annotating cell types for sc-/sn-RNA-Seq Analysis in 10X Genomics data

## Usage

To run the script in this module from the command line sequentially, use:

```
bash run-cell-types-annotation.sh
```

`run-cell-types-annotation.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `project_parameters.Config.yaml` located at the `root_dir`

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

