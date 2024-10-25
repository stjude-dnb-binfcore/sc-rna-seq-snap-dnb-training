# How to run the Singularity container for the Single cell RNA Seq Snap workflow (ScRNASeqSnap)

We have generated definition and singularity files that contain all tools, packages, and dependencies necessary to run the code and analyses modules in the `sc-rna-seq-snap` repository. 


## To use the container in an R interactive session on HPC:

From an interactive node on HPC, the user can open an R interactive session. Please modify memory and resources as needed for the analysis module to run.
```
bsub -P hpcf_interactive -J hpcf_interactive -n 1 -q standard -R "rusage[mem=4G]" -Is "bash"
```

```
module load singularity/4.1.1
```


1. Clone the repository
```
git clone https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap.git
```

and navigate to the directory with the files needed to run the container:
```
cd run-container
```


2. Pull the singularity container 🚧🚧🚧
```
singularity pull library://sc-rna-seq-snap-container.sif
```


3. Start the singularity container
```
bash run-container.sh
```

The `run-container.sh` is running at `IP_ADDR:PORT`. When RStudio launches, please click "Session" -> "Restart R".


4. Build container (if needed)

If the user does not have access to the `sc-rna-seq-snap-container.sif`, they can build their own. 
User can rename the `.sif` file, if they want to (not needed).
```
singularity build sc-rna-seq-snap-container.sif rstudio-v4.4.0-seurat-v4.4.0.def
```

Then, the user can start the container as explained in the step (3).


## To use the container from the command line on HPC: 🚧🚧🚧






## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
