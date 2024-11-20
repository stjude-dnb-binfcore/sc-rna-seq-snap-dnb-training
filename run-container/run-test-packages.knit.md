---
title: "Data exploratory analysis"
author: "Antonia Chroni for SJCRH DNB_BINF_Core"
papersize: a4
fontsize: 11pt
links-as-notes: true
output:
  html_document:
    toc: TRUE
    toc_float: TRUE
    df_print: paged
    code_folding: hide
    toc_depth: 2
    highlight: tango
    number_sections: TRUE
  pdf_document:
    toc: TRUE
    highlight: tango
    number_sections: TRUE
    latex_engine: lualatex
    keep_tex: FALSE
    fig_caption: yes
    fig_crop: no
    fig_height: 2
    fig_width: 3
    toc_depth: 2
always_allow_html: TRUE
urlcolor: blue
linkcolor: black
citecolor: blue
geometry: margin=1in
header-includes: 
  - \usepackage{titling}
  - \usepackage{fancyhdr}
  - \usepackage{graphicx}
  - \usepackage{float}
---

# Set up






# Session Info


```
## R version 4.4.0 (2024-04-24)
## Platform: x86_64-pc-linux-gnu
## Running under: Ubuntu 22.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Etc/UTC
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] yaml_2.3.10                 tinytex_0.54               
##  [3] lubridate_1.9.3             dplyr_1.1.4                
##  [5] readr_2.1.5                 tidyr_1.3.1                
##  [7] tibble_3.2.1                tidyverse_2.0.0            
##  [9] tidytext_0.4.2              stringr_1.5.1              
## [11] SoupX_1.6.2                 shiny_1.9.1                
## [13] SeuratObject_5.0.2          Seurat_4.4.0               
## [15] R.utils_2.12.3              R.oo_1.27.0                
## [17] R.methodsS3_1.8.2           rlist_0.4.6.2              
## [19] reshape2_1.4.4              RColorBrewer_1.1-3         
## [21] purrr_1.0.2                 patchwork_1.3.0            
## [23] optparse_1.7.5              leiden_0.4.3.1             
## [25] knitr_1.49                  irlba_2.3.5.1              
## [27] igraph_2.1.1                harmony_1.2.1              
## [29] Rcpp_1.0.13-1               ggthemes_5.1.0             
## [31] ggrepel_0.9.6               ggpmisc_0.6.1              
## [33] ggpp_0.5.8-1                ggh4x_0.2.8                
## [35] GGally_2.2.1                future_1.34.0              
## [37] fs_1.6.5                    forcats_1.0.0              
## [39] flextable_0.9.7             flexmix_2.3-19             
## [41] lattice_0.22-6              data.table_1.16.2          
## [43] cowplot_1.1.3               clustree_0.5.1             
## [45] ggraph_2.2.1                numbat_1.4.2               
## [47] Matrix_1.7-0                RcppPlanc_1.0.0            
## [49] scooter_0.0.0.9004          infercnv_1.22.0            
## [51] SingleR_2.8.0               celldex_1.16.0             
## [53] scDblFinder_1.20.0          scater_1.34.0              
## [55] ggplot2_3.5.1               scuttle_1.16.0             
## [57] SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0
## [59] Biobase_2.66.0              GenomicRanges_1.58.0       
## [61] GenomeInfoDb_1.42.0         IRanges_2.40.0             
## [63] S4Vectors_0.44.0            BiocGenerics_0.52.0        
## [65] MatrixGenerics_1.18.0       matrixStats_1.4.1          
## [67] miQC_1.14.0                 remotes_2.5.0              
## [69] devtools_2.4.5              usethis_3.0.0              
## [71] BiocManager_1.30.25        
## 
## loaded via a namespace (and not attached):
##   [1] ica_1.0-3                 plotly_4.10.4            
##   [3] zlibbioc_1.52.0           tidyselect_1.2.1         
##   [5] bit_4.5.0                 doParallel_1.0.17        
##   [7] rjson_0.2.23              blob_1.2.4               
##   [9] urlchecker_1.0.1          S4Arrays_1.6.0           
##  [11] parallel_4.4.0            png_0.1-8                
##  [13] cli_3.6.3                 ggplotify_0.1.2          
##  [15] askpass_1.2.1             openssl_2.2.2            
##  [17] goftest_1.2-3             textshaping_0.4.0        
##  [19] BiocIO_1.16.0             bluster_1.16.0           
##  [21] officer_0.6.7             tokenizers_0.3.0         
##  [23] BiocNeighbors_2.0.0       uwot_0.2.2               
##  [25] curl_6.0.1                mime_0.12                
##  [27] evaluate_1.0.1            tidytree_0.4.6           
##  [29] coin_1.4-3                stringi_1.8.4            
##  [31] rjags_4-16                parallelDist_0.2.6       
##  [33] XML_3.99-0.17             httpuv_1.6.15            
##  [35] AnnotationDbi_1.68.0      magrittr_2.0.3           
##  [37] rappdirs_0.3.3            splines_4.4.0            
##  [39] getopt_1.20.4             logger_0.4.0             
##  [41] sctransform_0.4.1         ggbeeswarm_0.7.2         
##  [43] sessioninfo_1.2.2         DBI_1.2.3                
##  [45] HDF5Array_1.34.0          jquerylib_0.1.4          
##  [47] withr_3.0.2               systemfonts_1.1.0        
##  [49] rprojroot_2.0.4           xgboost_1.7.8.1          
##  [51] lmtest_0.9-40             tidygraph_1.3.1          
##  [53] formatR_1.14              rtracklayer_1.66.0       
##  [55] htmlwidgets_1.6.4         SparseArray_1.6.0        
##  [57] reticulate_1.40.0         zoo_1.8-12               
##  [59] XVector_0.46.0            hahmmr_1.0.0             
##  [61] UCSC.utils_1.2.0          RhpcBLASctl_0.23-42      
##  [63] timechange_0.3.0          foreach_1.5.2            
##  [65] fansi_1.0.6               caTools_1.18.3           
##  [67] ggtree_3.10.1             rhdf5_2.50.0             
##  [69] quantreg_5.99             janeaustenr_1.0.0        
##  [71] alabaster.schemas_1.6.0   gridGraphics_0.5-1       
##  [73] ellipsis_0.3.2            lazyeval_0.2.2           
##  [75] phyclust_0.1-34           survival_3.5-8           
##  [77] scattermore_1.2           BiocVersion_3.20.0       
##  [79] crayon_1.5.3              RcppAnnoy_0.0.22         
##  [81] progressr_0.15.0          tweenr_2.0.3             
##  [83] scistreer_1.2.0           later_1.3.2              
##  [85] ggridges_0.5.6            codetools_0.2-20         
##  [87] profvis_0.4.0             KEGGREST_1.46.0          
##  [89] Rtsne_0.17                limma_3.62.1             
##  [91] gdtools_0.4.1             Rsamtools_2.22.0         
##  [93] filelock_1.0.3            pkgconfig_2.0.3          
##  [95] xml2_1.3.6                spatstat.univar_3.1-1    
##  [97] GenomicAlignments_1.42.0  aplot_0.2.3              
##  [99] spatstat.sparse_3.1-0     alabaster.base_1.6.1     
## [101] ape_5.8                   viridisLite_0.4.2        
## [103] xtable_1.8-4              fastcluster_1.2.6        
## [105] plyr_1.8.9                httr_1.4.7               
## [107] tools_4.4.0               globals_0.16.3           
## [109] pkgbuild_1.4.5            beeswarm_0.4.0           
## [111] nlme_3.1-164              futile.logger_1.4.3      
## [113] lambda.r_1.2.4            dbplyr_2.5.0             
## [115] ExperimentHub_2.14.0      MatrixModels_0.5-3       
## [117] digest_0.6.37             farver_2.1.2             
## [119] tzdb_0.4.0                SnowballC_0.7.1          
## [121] yulab.utils_0.1.8         viridis_0.6.5            
## [123] glue_1.8.0                cachem_1.1.0             
## [125] BiocFileCache_2.14.0      polyclip_1.10-7          
## [127] generics_0.1.3            Biostrings_2.74.0        
## [129] mvtnorm_1.3-2             parallelly_1.39.0        
## [131] pkgload_1.4.0             statmod_1.5.0            
## [133] here_1.0.1                ragg_1.3.3               
## [135] ScaledMatrix_1.14.0       fontBitstreamVera_0.1.1  
## [137] pbapply_1.7-2             httr2_1.0.6              
## [139] spam_2.11-0               dqrng_0.4.1              
## [141] utf8_1.2.4                graphlayouts_1.2.1       
## [143] gtools_3.9.5              alabaster.se_1.6.0       
## [145] gridExtra_2.3             GenomeInfoDbData_1.2.13  
## [147] rhdf5filters_1.18.0       RCurl_1.98-1.16          
## [149] memoise_2.0.1             rmarkdown_2.29           
## [151] scales_1.3.0              gypsum_1.2.0             
## [153] RANN_2.6.2                fontLiberation_0.1.0     
## [155] spatstat.data_3.1-4       cluster_2.1.6            
## [157] spatstat.utils_3.1-1      hms_1.1.3                
## [159] fitdistrplus_1.2-1        munsell_0.5.1            
## [161] colorspace_2.1-1          rlang_1.1.4              
## [163] quadprog_1.5-8            DelayedMatrixStats_1.28.0
## [165] sparseMatrixStats_1.18.0  dotCall64_1.2            
## [167] ggforce_0.4.2             xfun_0.49                
## [169] alabaster.matrix_1.6.0    coda_0.19-4.1            
## [171] TH.data_1.1-2             iterators_1.0.14         
## [173] modeltools_0.2-23         abind_1.4-8              
## [175] libcoin_1.0-10            treeio_1.26.0            
## [177] ggsci_3.2.0               Rhdf5lib_1.28.0          
## [179] futile.options_1.0.1      bitops_1.0-9             
## [181] promises_1.3.0            RSQLite_2.3.8            
## [183] sandwich_3.1-1            DelayedArray_0.32.0      
## [185] compiler_4.4.0            alabaster.ranges_1.6.0   
## [187] beachmat_2.22.0           SparseM_1.84-2           
## [189] polynom_1.4-1             listenv_0.9.1            
## [191] fontquiver_0.2.1          edgeR_4.4.0              
## [193] AnnotationHub_3.14.0      BiocSingular_1.22.0      
## [195] tensor_1.5                MASS_7.3-60.2            
## [197] uuid_1.2-1                BiocParallel_1.40.0      
## [199] spatstat.random_3.3-2     R6_2.5.1                 
## [201] fastmap_1.2.0             multcomp_1.4-26          
## [203] fastmatch_1.1-4           vipor_0.4.7              
## [205] ROCR_1.0-11               ggstats_0.7.0            
## [207] rsvd_1.0.5                nnet_7.3-19              
## [209] gtable_0.3.6              phangorn_2.12.1          
## [211] KernSmooth_2.23-22        miniUI_0.1.1.1           
## [213] deldir_2.0-4              htmltools_0.5.8.1        
## [215] RcppParallel_5.1.9        bit64_4.5.2              
## [217] spatstat.explore_3.3-3    lifecycle_1.0.4          
## [219] zip_2.3.1                 restfulr_0.0.15          
## [221] sass_0.4.9                vctrs_0.6.5              
## [223] spatstat.geom_3.3-4       scran_1.34.0             
## [225] ggfun_0.1.7               sp_2.1-4                 
## [227] future.apply_1.11.3       bslib_0.8.0              
## [229] pillar_1.9.0              gplots_3.2.0             
## [231] metapod_1.14.0            locfit_1.5-9.10          
## [233] jsonlite_1.8.9            argparse_2.2.3
```

