# How to run the Singularity container for the Single cell RNA Seq Snap workflow (ScRNASeqSnap)

We have generated definition and singularity files that contain all tools, packages, and dependencies necessary to run the code and analyses modules in the `sc-rna-seq-snap` repository. 


## To use the container in an R interactive session on HPC:

1. Clone the repository
```
git clone https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap.git
```

2. Pull the singularity container and start the singularity container, from the `sc-rna-seq-snap` folder

From an interactive node on HPC, user can open an R interactive session and run:
```
module load singularity/4.1.1
```

🚧🚧🚧
```
singularity pull library://sc-rna-seq-snap-container.sif bash ./run-container/run-container.shRStudio
```
🚧🚧🚧

When RStudio launches, please click "Session" -> "Restart R".


If the user does not have access to the `sc-rna-seq-snap-container.sif`, they can build their own. User can rename the `.sif` file, if they want to (not needed).

```
singularity build ./run-container/sc-rna-seq-snap-container.sif ./run-container/rstudio-v4.4.0-seurat-v4.4.0.def
```

```
build ./run-container/sc-rna-seq-snap-container.sif ./run-container/rstudio-v4.4.0-seurat-v4.4.0.def bash ./run-container/run-container.shRStudio
```


`run-container.shRStudio` is running at `IP_ADDR:PORT`.


## To use the container from the command line on HPC:

🚧🚧🚧





## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
