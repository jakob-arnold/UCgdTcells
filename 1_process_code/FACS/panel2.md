FACS (PBMC): Processing
================

``` r
fcs <- read.flowSet(path="../../data_raw/FACS")
markernames(fcs)
```

    ##        FJComp-FL10-A        FJComp-FL11-A        FJComp-FL13-A 
    ##       "GzmB BV510-A"      "CD127 BV605-A"       "CD28 BV786-A" 
    ##         FJComp-FL2-A         FJComp-FL4-A         FJComp-FL5-A 
    ##         "TCF-1 PE-A" "GzmA PC5.5_BB700-A"         "GzmK PC7-A"

``` r
markernames(fcs)[] <- c("GZMB", "CD127", "CD28", "TCF-1", "GZMA", "GZMK")
markernames(fcs)
```

    ## FJComp-FL10-A FJComp-FL11-A FJComp-FL13-A  FJComp-FL2-A  FJComp-FL4-A 
    ##        "GZMB"       "CD127"        "CD28"       "TCF-1"        "GZMA" 
    ##  FJComp-FL5-A 
    ##        "GZMK"

``` r
uc <- c("7516", "6919", "7355", "6344", "6588", "7658", "6593", "6324", "5189")

meta <- data.frame(file_name=list.files("../../data_raw/FACS"))%>%
  separate(file_name, into=c("panel", "TRDV", "id"), sep="_|\\.fcs", remove=F)%>%
   mutate(disease=case_when(
     id %in% uc ~ "UC",
     .default="HD"
   ))%>%
  mutate(ID=paste(disease, id, sep="_"))

vd1 <- prepData(
  fcs,
  md=meta,
  md_cols=list(file = "file_name", id="ID", factors=c("disease", "TRDV")),
  transform=T, cofactor=150, by_time=F, FACS=T
)

vd1@colData@listData[["ID"]] <- vd1@colData@listData[["sample_ID"]]

summary(vd1$disease)
```

    ##    HD    UC 
    ##  6743 15908

``` r
summary(vd1$sample_id)
```

    ## HD_0479 HD_3510 HD_3522 HD_4198 HD_4414 HD_5109 HD_7654 UC_5189 UC_6324 UC_6344 
    ##     948     834    2007    1021     136     197    1600     654    7747     872 
    ## UC_6588 UC_6593 UC_6919 UC_7355 UC_7516 UC_7658 
    ##     849    3528     363     488     385    1022

# Clustering

``` r
vd1 <- runDR(vd1, "UMAP", features=rownames(vd1), cells=1000)
```

``` r
set.seed(1337)
vd1 <- cluster(vd1, features=rownames(vd1), seed=1337)
```

``` r
saveRDS(vd1, "../../data_processed/FACS/vd1.Rds")
```

``` r
sessionInfo()
```

    ## R version 4.4.3 (2025-02-28)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Berlin
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] CATALYST_1.28.0             SingleCellExperiment_1.26.0
    ##  [3] HDCytoData_1.24.0           flowCore_2.16.0            
    ##  [5] SummarizedExperiment_1.34.0 Biobase_2.64.0             
    ##  [7] GenomicRanges_1.56.0        GenomeInfoDb_1.40.0        
    ##  [9] IRanges_2.38.0              S4Vectors_0.42.0           
    ## [11] MatrixGenerics_1.16.0       matrixStats_1.5.0          
    ## [13] ExperimentHub_2.12.0        AnnotationHub_3.12.0       
    ## [15] BiocFileCache_2.12.0        dbplyr_2.5.0               
    ## [17] BiocGenerics_0.50.0         lubridate_1.9.3            
    ## [19] forcats_1.0.0               stringr_1.5.1              
    ## [21] dplyr_1.1.4                 purrr_1.0.2                
    ## [23] readr_2.1.5                 tidyr_1.3.1                
    ## [25] tibble_3.3.0                ggplot2_3.5.1              
    ## [27] tidyverse_2.0.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3          rstudioapi_0.16.0          
    ##   [3] jsonlite_2.0.0              shape_1.4.6.1              
    ##   [5] magrittr_2.0.3              ggbeeswarm_0.7.2           
    ##   [7] TH.data_1.1-2               farver_2.1.2               
    ##   [9] rmarkdown_2.27              GlobalOptions_0.1.2        
    ##  [11] zlibbioc_1.50.0             vctrs_0.6.5                
    ##  [13] DelayedMatrixStats_1.26.0   memoise_2.0.1              
    ##  [15] rstatix_0.7.2               htmltools_0.5.8.1          
    ##  [17] S4Arrays_1.4.1              plotrix_3.8-4              
    ##  [19] curl_6.4.0                  BiocNeighbors_1.22.0       
    ##  [21] broom_1.0.6                 SparseArray_1.4.5          
    ##  [23] plyr_1.8.9                  sandwich_3.1-0             
    ##  [25] zoo_1.8-12                  cachem_1.1.0               
    ##  [27] igraph_2.0.3                lifecycle_1.0.4            
    ##  [29] iterators_1.0.14            pkgconfig_2.0.3            
    ##  [31] rsvd_1.0.5                  Matrix_1.7-3               
    ##  [33] R6_2.6.1                    fastmap_1.2.0              
    ##  [35] GenomeInfoDbData_1.2.12     clue_0.3-65                
    ##  [37] digest_0.6.35               colorspace_2.1-0           
    ##  [39] ggnewscale_0.4.10           AnnotationDbi_1.66.0       
    ##  [41] scater_1.32.0               irlba_2.3.5.1              
    ##  [43] RSQLite_2.3.6               beachmat_2.20.0            
    ##  [45] ggpubr_0.6.0                filelock_1.0.3             
    ##  [47] cytolib_2.16.0              nnls_1.5                   
    ##  [49] colorRamps_2.3.4            timechange_0.3.0           
    ##  [51] httr_1.4.7                  polyclip_1.10-6            
    ##  [53] abind_1.4-8                 compiler_4.4.3             
    ##  [55] bit64_4.0.5                 withr_3.0.0                
    ##  [57] doParallel_1.0.17           ConsensusClusterPlus_1.68.0
    ##  [59] backports_1.4.1             BiocParallel_1.38.0        
    ##  [61] viridis_0.6.5               carData_3.0-5              
    ##  [63] DBI_1.2.2                   ggforce_0.4.2              
    ##  [65] ggsignif_0.6.4              MASS_7.3-64                
    ##  [67] drc_3.0-1                   rappdirs_0.3.3             
    ##  [69] DelayedArray_0.30.1         rjson_0.2.21               
    ##  [71] FlowSOM_2.12.0              gtools_3.9.5               
    ##  [73] tools_4.4.3                 vipor_0.4.7                
    ##  [75] beeswarm_0.4.0              glue_1.8.0                 
    ##  [77] grid_4.4.3                  Rtsne_0.17                 
    ##  [79] reshape2_1.4.4              cluster_2.1.8              
    ##  [81] generics_0.1.3              gtable_0.3.5               
    ##  [83] tzdb_0.4.0                  data.table_1.15.4          
    ##  [85] hms_1.1.3                   ScaledMatrix_1.12.0        
    ##  [87] BiocSingular_1.20.0         car_3.1-2                  
    ##  [89] XVector_0.44.0              RcppAnnoy_0.0.22           
    ##  [91] ggrepel_0.9.5               BiocVersion_3.19.1         
    ##  [93] foreach_1.5.2               pillar_1.11.0              
    ##  [95] splines_4.4.3               circlize_0.4.16            
    ##  [97] tweenr_2.0.3                lattice_0.22-6             
    ##  [99] survival_3.8-3              bit_4.0.5                  
    ## [101] RProtoBufLib_2.16.0         tidyselect_1.2.1           
    ## [103] ComplexHeatmap_2.20.0       scuttle_1.14.0             
    ## [105] Biostrings_2.72.0           knitr_1.46                 
    ## [107] gridExtra_2.3               xfun_0.44                  
    ## [109] stringi_1.8.7               UCSC.utils_1.0.0           
    ## [111] yaml_2.3.8                  evaluate_0.23              
    ## [113] codetools_0.2-20            BiocManager_1.30.23        
    ## [115] cli_3.6.5                   uwot_0.2.2                 
    ## [117] munsell_0.5.1               Rcpp_1.1.0                 
    ## [119] png_0.1-8                   XML_3.99-0.16.1            
    ## [121] parallel_4.4.3              blob_1.2.4                 
    ## [123] sparseMatrixStats_1.16.0    viridisLite_0.4.2          
    ## [125] mvtnorm_1.2-5               ggridges_0.5.6             
    ## [127] scales_1.3.0                crayon_1.5.3               
    ## [129] GetoptLong_1.0.5            rlang_1.1.6                
    ## [131] cowplot_1.1.3               KEGGREST_1.44.0            
    ## [133] multcomp_1.4-25
