Smillie 2019: Processing
================

``` r
library(Seurat)
library(scCustomize)
library(harmony)
library(Matrix)
library(tidyverse)
library(data.table)
```

# Data

``` r
# Meta df
meta <- fread("../../data_raw/Smillie2019/meta.data.txt")
meta <- meta[-1]

# Immune Cells
imm <- readMM("../../data_raw/Smillie2019/gene_sorted-Imm.matrix.mtx")
rownames(imm) <- readLines("../../data_raw/Smillie2019/Imm.genes.tsv")
colnames(imm) <- readLines("../../data_raw/Smillie2019/Imm.barcodes2.tsv")

imm <- CreateSeuratObject(imm)
```

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

    ## Warning: Data is of class dgTMatrix. Coercing to dgCMatrix.

``` r
imm@meta.data <- imm@meta.data%>%
  rownames_to_column("NAME")%>%
  left_join(meta, by="NAME")%>%
  column_to_rownames("NAME")

# Epithelial Cells
ep <- readMM("../../data_raw/Smillie2019/gene_sorted-Epi.matrix.mtx")
rownames(ep) <- readLines("../../data_raw/Smillie2019/Epi.genes.tsv")
colnames(ep) <- readLines("../../data_raw/Smillie2019/Epi.barcodes2.tsv")

ep <- CreateSeuratObject(ep)
```

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')
    ## Warning: Data is of class dgTMatrix. Coercing to dgCMatrix.

``` r
ep@meta.data <- ep@meta.data%>%
  rownames_to_column("NAME")%>%
  left_join(meta, by="NAME")%>%
  column_to_rownames("NAME")

# Stromal Cells
fb <- readMM("../../data_raw/Smillie2019/gene_sorted-Fib.matrix.mtx")
rownames(fb) <- readLines("../../data_raw/Smillie2019/Fib.genes.tsv")
colnames(fb) <- readLines("../../data_raw/Smillie2019/Fib.barcodes2.tsv")

fb <- CreateSeuratObject(fb)
```

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')
    ## Warning: Data is of class dgTMatrix. Coercing to dgCMatrix.

``` r
fb@meta.data <- fb@meta.data%>%
  rownames_to_column("NAME")%>%
  left_join(meta, by="NAME")%>%
  column_to_rownames("NAME")
```

``` r
ep$compartment <- "Epithelial"
imm$compartment <- "Immune"
fb$compartment <- "Stromal"
```

``` r
genes <- Reduce(intersect, list(
  rownames(ep),
  rownames(imm),
  rownames(fb)
  )
  )

all <- merge(
  ep[genes,],y=c(imm[genes,], fb[genes,])
  )
```

``` r
all <- all%>%
  JoinLayers()

rm(imm, ep, fb, meta, genes)
gc()
```

    ##             used   (Mb) gc trigger    (Mb)   max used    (Mb)
    ## Ncells   4086906  218.3    7625631   407.3    7625631   407.3
    ## Vcells 592830614 4523.0 3066845960 23398.2 2857906418 21804.1

# QC

``` r
all$percent.mt <- PercentageFeatureSet(all, pattern="^MT-")
```

``` r
VlnPlot_scCustom(
  all,
  features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by="compartment",
  pt.size=0
)
```

    ## Warning: Default search for "data" layer in "RNA" assay yielded no results;
    ## utilizing "counts" layer instead.

![](Li21_Smillie19_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
all <- subset(
  all, subset=
    (compartment=="Epithelial"&percent.mt<75)|
    (compartment!="Epithelial"&percent.mt<25)&
    nFeature_RNA>500&
    nFeature_RNA<5000&
    nCount_RNA<40000
)

gc()
```

    ##             used   (Mb) gc trigger    (Mb)   max used    (Mb)
    ## Ncells   4319068  230.7    7625632   407.3    7625632   407.3
    ## Vcells 559331466 4267.4 2453476768 18718.6 2857906418 21804.1

# Processing

``` r
all <- all%>%
  NormalizeData()%>%
  FindVariableFeatures()%>%
  ScaleData()%>%
  FindVariableFeatures()%>%
  RunPCA()%>%
  RunUMAP(dims=1:30)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

``` r
pb <- all%>%
  subset(subset=Health!="Non-inflamed")%>%
  AggregateExpression(
    group.by=c("Health", "compartment", "Subject"),
    return.seurat=T
    )
```

# Save Object

``` r
saveRDS(pb, "../../data_processed/Smillie19/PB.Rds")
```
