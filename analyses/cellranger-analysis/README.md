# Pipeline for running and summarizing Cell Ranger count for single or multiple libraries for sc-/sn-RNA-Seq Analysis in 10X Genomics data

## Usage

To run all the scripts in this module from the command line sequentially, use:

```
bsub < submit-multiple-jobs.sh
```

The `submit-multiple-jobs.sh` script is designed to be run as if it was called from this module directory even when called from outside of this directory. 
   - Step 1: To run the `j1.sh` script to align single or multiple libraries in parallel, i.e., `run-cellranger-analysis`. 
   - Step 2: To run `j2.sh` to summarize alignment results, i.e., `summarize-cellranger-analysis`. The latter script will be on hold and executed once all libraries are aligned and `j1.sh` is complete. This is been taken care of by `waiter.sh` script.

Parameters according to the project and analysis strategy will need to be specified in the following scripts:
- `j1.sh`: 
   - `--file`: define path with `project_metadata` file. This is defined by default in the `./input/` dir within the `/sc-rna-seq-snap/analyses/cellranger-analysis` module. If the file lives in another folder, this path should be modified accordingly within the script. The file needs to be in `*.txt` file format. In addition, there is a code line to convert `*.tsv` file to `*.txt` file, if needed. User can specify the file path in there for the conversion. This file can contain either one or multiple samples as long as the following minimum parameters are provided: `ID`, `SAMPLE`, and `FASTQ` columns.
   - `--transcriptome`: define reference genome to be used for the alignment. This points out to the path with the reference genome. There are 4 options accommodating for human, mouse, and dual index genomes: `GRCh38`, `GRCm39`, `GRCh38ANDGRCm39`, and `GRCh38_GFP_tdTomato`, respectively. If a different reference genome is required (not provided here), this can be added as an argument option in the `./util/run_cellranger.py` along with the path. 
  - `--output_dir`: User can change this accordingly if they want to run and compare alignments with CellRanger, e.g., by using `--force_cells=8000` to constraint the number of cells to be used for alignment. We recommend to run by default, and after careful assessment to edit parameters. We have found that the default parameters set up here work well for most of the cases.

- `./util/run_cellranger.py`: The path to the reference genome is required. Please submit an [issue](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues) to request the path to the path with all reference genomes maintained by our team at the Bioinformatics core at DNB. Otherwise, specify the path to the reference genome of your preference.


- `j2.sh`: No modifications are needed by the user, unless the user modifies the `--output_dir` path in the `j1.sh`, then they will need to adjust accordingly the `--dir` flag. 

- `submit-multiple-jobs.sh`: Here, the user will need to set up the absolute path for the directory in `#BSUB -cwd` and `prefix` parameters. If different alignment strategies are used, user should rename `${prefix}/results/02_cellranger_count/DefaultParameters` accordingly, e.g. `${prefix}/results/02_cellranger_count/ForcedCells8000Parameters`.

- `waiter.sh`: Here, the user will need to set up the absolute path for the directory in `#BSUB -cwd` and `prefix` parameters. Also, user needs to replace `DST` with the Sample ID used for the samples of the project. This requires that all of the Sample IDs are the same, e.g., DST800, DST802, DST811 and so on.


Default memory for running the alignment is set up to 16GB. In case more memory is needed to run a specific sample, this can be modified here:
- `./util/run_cellranger.py`: User can search the following and change accordingly `-n 4 -R "rusage[mem=4000] span[hosts=1]`, if that is necessary. However, we have found that these resources work well for most data. 


## Folder content
This folder contains scripts tasked to run and summarize Cell Ranger count for single or multiple libraries for sc-/sn-RNA-Seq Analysis in 10X Genomics data across the project. For more information and updates, please see [Cell Ranger support page](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-gex-count).

This module uses CellRanger v8.0.1 for the alignment.


## Folder structure 

The structure of this folder is as follows:

```
├── input
|   └── project_metadata.txt
├── j1.sh
├── j2.sh
├── README.md
├── results
|   ├── 01_logs
|   ├── 02_cellranger_count
|   |   └── DefaultParameters
|   └── 03_cellranger_count_summary
├── submit-multiple-jobs.sh
├── util
|   ├── run_cellranger.py
|   └── summarize_cellranger_results.py
└── waiter.sh
```
