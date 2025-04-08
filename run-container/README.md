# Running the Container for scRNA-Seq Snap Workflow

We provide a Dockerfile and Definition file that include all tools, packages, and dependencies necessary for running the sc-rna-seq-snap analysis modules. These are customized for `Rstudio/R v4.4.0` and `Seurat v4.4.0`.


## Table of Contents

1. [Running the Container on HPC](#running-the-container-on-hpc)
   - [Start an Interactive Session](#start-an-interactive-session)
   - [Load the Singularity Module](#load-the-singularity-module)
   - [Pull the Singularity Container](#pull-the-singularity-container)
   - [Start the Singularity Container](#start-the-singularity-container)
     - [Running Analysis Modules via LSF](#running-analysis-modules-via-lsf)
     - [Running from the Terminal](#running-from-the-terminal)
     - [Running from RStudio](#running-from-rstudio)
     - [Fixing Issues with RStudio Server](#fixing-issues-with-rstudio-server)
   - [Building the Container (If Needed)](#building-the-container-if-needed)
   
2. [Running the Container Outside HPC (Docker)](#running-the-container-outside-hpc-docker)
   - [Pull the Docker Container](#pull-the-docker-container)
   - [Start the Docker Container](#start-the-docker-container)


## Running the Container on HPC

### Start an Interactive Session

Open an interactive node on the HPC and adjust memory/resources as needed:

```
bsub -P hpcf_interactive -J hpcf_interactive -n 2 -q standard -R "rusage[mem=16G]" -Is "bash"
```

### Load the Singularity Module

Please note that a version of Singularity is installed by default on all the cluster nodes at St Jude HPC. Otherwise the user needs to ensure and load Singularity module by running the following on HPC:

```
module load singularity/4.1.1
```

### Pull the Singularity Container

1. Pull the singularity container from the `sc-rna-seq-snap` root_dir

```
singularity pull docker://achronistjude/rstudio_4.4.0_seurat_4.4.0:latest
```


### Start the Singularity Container

a. Running Analysis Modules via LSF

All analysis modules (except for `.analyses/cellranger-analysis`) are designed to be run while executing the container. User only needs to run the lsf script as described in the `README.md` files in each analysis module.


b. Running from the Terminal

User can run analysis module while on interactive node after executing the container:

```
bash run-terminal.sh
```

Then user may navigate to their module of interest, `./sc-rna-seq-snap/analyses/<module_of_interest>`. For example:

```
cd ./sc-rna-seq-snap/analyses/upstream-analysis
bash run-upstream-analysis.sh
```


c. Running from RStudio

User can also run analyses via Rstudio OnDemand after executing the container:

```
bash run-rstudio.sh
```

The `run-rstudio.sh` is running at `IP_ADDR:PORT`. When RStudio launches, please click "Session" -> "Restart R" (at the RStudio web session). For St Jude users, we advice to disconnect from CloudFlare WARP as this might lead to unstable behavior while on VPN.

Again, the user can navigate to their module of interest and explore/run their analyses.


d. Fixing Issues with RStudio Server

If you encounter issues during this step related to RStudio Server and specifically to an invalid secure cookie error. This might be an issue with how the secure cookie is being handled during an HTTP request. In this case, please check if the following directories have been generated and if so, remove them:

```
rm -r .cache/
rm -r .config/
rm -r .local/
rm -r rstudio-container-tmp/
```

These folders cache history and user info. Then, kill the interactive session, start a new one, and hopefully, it works! 🎉


### Building the Container (If Needed)

If the user does not have access to the `rstudio_4.4.0_seurat_4.4.0_latest.sif`, they can build their own. 
User can rename the `.sif` file, if they want to (not needed). Run the following from the `./run-container` dir:
```
singularity build rstudio_4.4.0_seurat_4.4.0_latest.sif rstudio_r_4.4.0_seurat_4.4.0.def
```

Then, the user can start the container as explained in the step (3).


## Running the Container Outside HPC (Docker)

1. Pull the Docker Container from the `sc-rna-seq-snap` root_dir:

```
docker pull docker://achroni/rstudio_4.4.0_seurat_4.4.0:latest
```

2. Start the Docker Container

Run the Docker container from the terminal:

```
docker run --platform linux/amd64 --name review -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/sc-rna-seq-snap docker://achroni/rstudio_4.4.0_seurat_4.4.0:latest
```

Start the container and open a terminal:

```
docker container start review
docker exec -ti review bash
```

Navigate to your module of interest and run the analysis:

```
cd ./sc-rna-seq-snap/analyses/upstream-analysis
bash run-upstream-analysis.sh
```


## Authors

Antonia Chroni, PhD ([@AntoniaChroni](https://github.com/AntoniaChroni)) and 
Walid Abu Al-Afia ([@walidabualafia](https://github.com/walidabualafia)).


## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
