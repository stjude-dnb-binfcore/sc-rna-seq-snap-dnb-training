# rstudio_containers

This is an example startup sequence:
 
```
ssh hpc.stjude.org
hpcf_interactive -n 4
cd ./sc-rna-seq-snap/run-container/working_dir/
bash run-container.sh
```

When RStudio launches, please click "Session" -> "Restart R".

The sif file is required and needs to be added in the same directory where `run-container.sh` lives. Please submit an [issue](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues) to request the file maintained by our team at the Bioinformatics core at DNB. 





## Contact

Contributions, issues, and feature requests are welcome! Please feel free to check [issues](https://github.com/stjude-dnb-binfcore/sc-rna-seq-snap/issues).

---

*These tools and pipelines have been developed by the Bioinformatic core team at the [St. Jude Children's Research Hospital](https://www.stjude.org/). These are open access materials distributed under the terms of the [BSD 2-Clause License](https://opensource.org/license/bsd-2-clause), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
