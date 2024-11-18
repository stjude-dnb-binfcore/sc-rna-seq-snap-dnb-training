# How to use analysis modules in the Single cell RNA Seq Snap workflow (ScRNASeqSnap)

This repository contains tools and workflows for analyzing single cell and single nuclei RNA (sc/snRNA) data from 10X sequencing technology. 

We recommend the following order in running analysis modules which is up to customization based on the user's needs.
1. `fastqc-analysis` module (help=`Required`. Pipeline for FastQC quality control tool for high throughput sequence data analysis.)
2. `cellranger-analysis` module (help=`Required`. Pipeline for running and summarizing Cell Ranger count for single or multiple libraries.)
3. `upstream-analysis` module (help=`Required`. Pipeline for estimating QC metrics and filtering low quality cells.)
4. `integrative-analysis` module (help=`Required`. Pipeline for Integrative analysis.)
5. `cluster-cell-calling` module (help=`Required`. Pipeline for cluster cell calling and gene marker analysis.)
6. `cell-contamination-removal-analysis` module (help=`Optional`. To remove clusters and repeat steps (4) and (5), e.g. for PDX experiments)
7. `cell-types-annotation` module (help=`Required`. Pipeline for annotating cell types.)


## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
