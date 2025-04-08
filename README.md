<p align="center";">
  <img src="figures/img/SCRNA_Logo_Primary.png" alt="ScRNASeqSnap repository logo" width="560px" />
</p>
<p align="center";">
  <a href="https://www.repostatus.org/#active">
    <img src="https://www.repostatus.org/badges/latest/active.svg?style=for-the-badge" alt="The project has reached a stable, usable state and is being actively developed." />
  </a>
  <a href="https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap">
    <img src="https://img.shields.io/badge/version-2.0.1-brightgreen" alt="Version" />
  </a>
</p>


# Single cell RNA Seq Snap workflow (ScRNASeqSnap)

Snap is a comprehensive suite of tools and workflows for analyzing single-cell and single-nucleus RNA (sc/snRNA) data from 10X Genomics sequencing technology supporting human, mouse, and dual genome cohorts.

A pre-built [Docker image](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/blob/main/run-container/README.md) is available, containing all the necessary tools, packages, and dependencies to seamlessly run the code and analysis modules. 

Snap is an initiative of the [Bioinformatics Core](https://www.stjude.org/research/departments/developmental-neurobiology/shared-resources/bioinformatic-core.html) at the Department of Developmental Neurobiology at the St. Jude Children's Research Hospital.


## Tutorial on how to use the single cell RNA seq Snap pipeline

For guidelines and best practices on effectively using the sc-rna-seq-snap pipeline, refer to the [Tutorial and documentation for the snap pipeline](https://github.com/stjude-dnb-binfcore/trainings/blob/main/courses/sc-rna-seq-snap-repo/tutorial/snap-tutorial-docs/Documentation-snap-repo-tutorial-2025-02-18.pdf).

For more detailed information on accessing the code, running it, or requesting memory from the HPCF cluster, check our [wiki page](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/wiki).


### Below is the main directory structure listing the analyses and data files used in this repository

```
├── analyses
|  ├── cell-contamination-removal-analysis
|  ├── cell-types-annotation
|  ├── cellranger-analysis
|  ├── clone-phylogeny-analysis
|  ├── cluster-cell-calling
|  ├── fastqc-analysis
|  ├── integrative-analysis
|  ├── README.md
|  ├── rshiny-app
|  └── upstream-analysis
├── figures
├── LICENSE
├── project_parameters.Config.yaml
├── README.md
├── run-container
├── run-rstudio.sh
├── run-terminal.sh
└── SECURITY.md
```

## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues).

---

*These tools and pipelines have been developed by the Bioinformatics core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
