---
title: "10X 5P scRNA-seq of gd T Cells: Processing"
output: github_document
---



##-----------------------------------------------------------------------------
library(Seurat)
library(harmony)
library(SCpubr)
library(scCustomize)
library(tidyverse)
library(future.apply)
library(scRepertoire)


# Data Processing
## Reading in Data
##-----------------------------------------------------------------------------
sample_dir <- list.dirs("../../data_raw/10X", recursive=T)
sample_dir <- sample_dir[grepl("VD", sample_dir)]

process_sample <- function(sample_dir) {
  
  batch <- str_extract(sample_dir, "R[1-2]")
  sample_name <- str_extract(sample_dir, "VD[^/]+")
  sample_data <- Read10X(data.dir=sample_dir)
  
  rownames(sample_data$'Antibody Capture') <- sub(
    "^Hu\\.Hashtag_", "hash", rownames(sample_data$'Antibody Capture'))
  
  seurat_obj <- CreateSeuratObject(sample_data$'Gene Expression', project=sample_name)
  
  seurat_obj[["HTO"]] <- CreateAssayObject(sample_data$
    'Antibody Capture'[rownames(sample_data$'Antibody Capture') %in% paste0("hash",1:6), ])
  
  seurat_obj$batch <- batch
  
  seurat_obj <- RenameCells(seurat_obj, add.cell.id=paste0(batch, "_",sample_name))
  
  return(seurat_obj)
}

# workers = number of cores available
plan(multisession, workers=8)
all <- future_lapply(sample_dir, process_sample)
plan(sequential)

rm(sample_dir, process_sample)


## Demultiplexing (based on Hashing)
##-----------------------------------------------------------------------------
demux <- function(seurat_obj){
  DefaultAssay(seurat_obj) <- "HTO"
  seurat_obj <- seurat_obj%>%NormalizeData(normalization.method="CLR")%>%
    HTODemux(positive.quantile=.99)
  
  print(summary(as.factor(seurat_obj@meta.data$HTO_classification.global)))
  
  seurat_obj <- subset(seurat_obj, subset=HTO_classification.global=="Singlet")
  
   DefaultAssay(seurat_obj) <- "RNA"
  
  return(seurat_obj)
}

# workers = number of cores available
plan(multisession, workers=8)
all <- future_lapply(all, demux, future.seed=1337)
plan(sequential)

rm(demux)

## Adding Meta Data/Removing Unwanted genes
##-----------------------------------------------------------------------------
meta_function <- function(seurat_obj) {
  seurat_obj@meta.data <- seurat_obj@meta.data%>%
    mutate(
      disease=case_when(
        (batch=="R1" & hash.ID %in% c("hash1", "hash2", "hash3")) |
        (batch=="R2" & hash.ID %in% c("hash1", "hash2")) ~ "HD",
        .default="UC"
      )
    )%>%
    mutate(diseaseID=case_when(
      batch=="R1" & disease=="HD" ~ paste0(disease, str_remove(hash.ID, "hash")),
      batch=="R2" & disease=="HD" ~ paste0(disease, as.numeric(str_remove(hash.ID, "hash"))+3),
      batch=="R1" & disease=="UC" ~ paste0(disease, as.numeric(str_remove(hash.ID, "hash"))-3),
      batch=="R2" & disease=="UC" ~ paste0(disease, as.numeric(str_remove(hash.ID, "hash"))+1)
    ))
    #   unite("diseaseID", disease, batch, hash.ID, sep="_", remove=FALSE)%>%
    # mutate(diseaseID=str_replace(diseaseID, pattern="hash", replacement=""))
    

  seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern="^MT-")
  
  # Remove mitochondrial/ribosomal/un-annoted genes from dataset
  genes_to_remove <- grep(rownames(seurat_obj), pattern="^MT-|^RP[SL]|^ENSG[0-9]", value=T)
  seurat_obj <- subset(seurat_obj, features=setdiff(rownames(seurat_obj), genes_to_remove))
  Idents(seurat_obj) <- "all"

  return(seurat_obj)
}

all <- lapply(all, meta_function)

for (i in seq_along(all)){
  Idents(all[[i]]) <- "all"
}


## Merging objects (split by gd T subsets)
##-----------------------------------------------------------------------------
vd1 <- merge(all[[1]], all[[4]])%>%JoinLayers()
vd2 <- merge(all[[3]], all[[6]])%>%JoinLayers()

rm(i, meta_function, all)


# Vd1
## QC
##-----------------------------------------------------------------------------
p <- VlnPlot_scCustom(vd1, add.noise=F, pt.size=0,
                      features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))&
  theme(axis.title=element_blank())

p[[1]] <- p[[1]]+geom_hline(yintercept=c(500, 3000), linetype="dashed", color="red")+
  ggtitle("Number of Genes")
p[[2]] <- p[[2]]+geom_hline(yintercept=c(10000), linetype="dashed", color="red")+
  ggtitle("UMI counts")
p[[3]] <- p[[3]]+geom_hline(yintercept=c(10), linetype="dashed", color="red")+
  ggtitle("% Mitochondrial\nGenes")
p

# Save QC Plots 
ggsave("../../2_figures_code/figures/10X_Vd1_VlnPlots_QC.pdf", width=6.5, height=3.5)
ggsave("../../2_figures_code/figures/10X_Vd1_VlnPlots_QC.png", dpi=600, width=6.5, height=3.5)

rm(p)

##-----------------------------------------------------------------------------
vd1 <- subset(vd1, subset=nFeature_RNA>500 & nFeature_RNA<3000 &
              nCount_RNA<10000 & percent.mt<10)


## Processing & Clustering
##-----------------------------------------------------------------------------
vd1 <- vd1%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA()%>%
  RunUMAP(dims=1:30)%>%FindNeighbors(dims=1:30)%>%FindClusters()

##-----------------------------------------------------------------------------
do_DimPlot(vd1, group.by="batch")
do_DimPlot(vd1, group.by="seurat_clusters", label=T)&NoLegend()

##-----------------------------------------------------------------------------
FeaturePlot_scCustom(vd1, features=c("MS4A1", "MS4A7", "CD14", "PPBP", "TRDV1","CD3E"), order=F)&
  NoAxes()&NoLegend()

##-----------------------------------------------------------------------------
VlnPlot_scCustom(vd1, features=c("nFeature_RNA", "percent.mt"), pt.size=0)


##-----------------------------------------------------------------------------
vd1 <- subset(vd1, subset=seurat_clusters %in% c("11", "21", "22", "24", "23", 
                                                 "13", "7", "17"), invert=T)


## Batch Correction (Batch + DonorID)
##-----------------------------------------------------------------------------
vd1 <- vd1%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA()%>%
  RunHarmony(c("batch", "diseaseID"))

##-----------------------------------------------------------------------------
ElbowPlot(vd1, reduction="harmony", 50)


##-----------------------------------------------------------------------------
vd1 <- vd1%>%RunUMAP(reduction="harmony", dims=1:20)%>%
  FindNeighbors(reduction="harmony", dims=1:20)%>%FindClusters(resolution=.4)

##-----------------------------------------------------------------------------
do_DimPlot(vd1, group.by="batch")
do_DimPlot(vd1, group.by="disease")
do_DimPlot(vd1, group.by="diseaseID")
do_DimPlot(vd1, group.by="seurat_clusters", label=T)&NoLegend()

## TCR
##-----------------------------------------------------------------------------
contig <- list("R1"=read_tsv( "../../data_raw/10X/R1/VD1/clones.tsv"),
               "R2"=read_tsv("../../data_raw/10X/R2/VD1/clones.tsv"))
contig <- loadContigs(contig, format = "MiXCR")


##-----------------------------------------------------------------------------
# Remove allele variant calls + scoring (*00(###))
contig <- lapply(contig, function(df) {
    df %>%
    mutate_at(paste0(c("v", "d", "j", "c"), "_gene"), ~str_replace(., "\\*.*", ""))
})


##-----------------------------------------------------------------------------
combined <- combineTCR(contig,  removeNA=T, filterMulti=T)
names(combined) <- names(contig)

for (i in seq_along(combined)) {
  combined[[i]]$barcode <- paste0(names(combined[i]), "_VD1_", combined[[i]]$barcode, "-1")
}


##-----------------------------------------------------------------------------
combined <- list(bind_rows(combined))

vd1 <- combineExpression(combined, vd1, proportion=F,cloneSize=c(
  Single=1, Small=5, Medium=25, Large=100, Hyper=100
))


##-----------------------------------------------------------------------------
trgv <- combineTCR(contig,  removeNA=F, filterMulti=T)

names(trgv) <- names(contig)

for (i in seq_along(trgv)) {
  trgv[[i]]$barcode <- paste0(names(trgv[i]), "_VD1_", trgv[[i]]$barcode, "-1")
}

trgv <- bind_rows(trgv)%>%
  separate(TCR1, into=c("TRGV", "TRGJ", "TRGC"), sep="\\.")%>%
  mutate(TRGV=factor(TRGV, levels=c(
    paste0("TRGV", 2:5), "TRGV5P", paste0("TRGV", 7:11)
  )))%>%
  select(barcode, TRGV)


##-----------------------------------------------------------------------------
vd1@meta.data <- vd1@meta.data%>%
  rownames_to_column("barcode")%>%
  left_join(trgv, by="barcode")%>%
  column_to_rownames("barcode")


##-----------------------------------------------------------------------------
rm(i, combined, contig, trgv)


# Vd2
## QC
##-----------------------------------------------------------------------------
p <- VlnPlot_scCustom(vd2, add.noise=F, pt.size=0,
                      features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))&
  theme(axis.title=element_blank())

p[[1]] <- p[[1]]+geom_hline(yintercept=c(200, 3000), linetype="dashed", color="red")+
  ggtitle("Number of Genes")
p[[2]] <- p[[2]]+geom_hline(yintercept=c(10000), linetype="dashed", color="red")+
  ggtitle("UMI counts")
p[[3]] <- p[[3]]+geom_hline(yintercept=c(15), linetype="dashed", color="red")+
  ggtitle("% Mitochondrial\nGenes")
p

# Save QC Plots 
ggsave("../../2_figures_code/figures/10X_Vd2_VlnPlots_QC.pdf", width=6.5, height=3.5)
ggsave("../../2_figures_code/figures/10X_Vd2_VlnPlots_QC.png", dpi=600, width=6.5, height=3.5)

rm(p)


##-----------------------------------------------------------------------------
vd2 <- subset(vd2, subset=nFeature_RNA>200 & nFeature_RNA<3000 & percent.mt<15)


## Processing & Clustering
##-----------------------------------------------------------------------------
vd2 <- vd2%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA()%>%
  RunUMAP(dims=1:30)%>%FindNeighbors(dims=1:30)%>%FindClusters()


##-----------------------------------------------------------------------------
do_DimPlot(vd2, group.by="seurat_clusters", label=T)&NoLegend()
FeaturePlot_scCustom(vd2, features=c("MS4A1", "MS4A7", "CD14", "PPBP", "CD3E", "TRDV2"), order=F)&
  NoLegend()&NoAxes()


##-----------------------------------------------------------------------------
VlnPlot_scCustom(vd2, features=c("nFeature_RNA", "percent.mt", "CD3E"), pt.size=0)&
  geom_boxplot(fill="white", width=.5, coef=0, outlier.shape=NA)


##-----------------------------------------------------------------------------
# Remove Contaminating Cells/low quality Cells
vd2 <- subset(vd2, subset=seurat_clusters %in% c("4", "8", "14", "16"), invert=T)


## Batch Correction (Batch + DonorID)
##-----------------------------------------------------------------------------
vd2 <- vd2%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA()%>%
  RunHarmony(c("batch", "diseaseID"))


##-----------------------------------------------------------------------------
ElbowPlot(vd2, reduction="harmony", 50)


##-----------------------------------------------------------------------------
vd2 <- vd2%>%RunUMAP(reduction="harmony", dims=1:20)%>%
  FindNeighbors(reduction="harmony", dims=1:20)%>%FindClusters()


##-----------------------------------------------------------------------------
do_DimPlot(vd2, group.by="batch")
do_DimPlot(vd2, group.by="disease")
do_DimPlot(vd2, group.by="seurat_clusters", label=T)&NoLegend()


## TCR
##-----------------------------------------------------------------------------
contig <- list("R1"=read_tsv( "../../data_raw/10X/R1/VD2/clones.tsv"),
               "R2"=read_tsv("../../data_raw/10X/R2/VD2/clones.tsv"))
contig <- loadContigs(contig, format = "MiXCR")



##-----------------------------------------------------------------------------
# Remove allele variant calls + scoring (*00(###))
contig <- lapply(contig, function(df) {
    df %>%
    mutate_at(paste0(c("v", "d", "j", "c"), "_gene"), ~str_replace(., "\\*.*", ""))
})



##-----------------------------------------------------------------------------
combined <- combineTCR(contig,  removeNA=T, filterMulti=T)
names(combined) <- names(contig)

for (i in seq_along(combined)) {
  combined[[i]]$barcode <- paste0(names(combined[i]), "_VD2_", combined[[i]]$barcode, "-1")
}

combined <- list(bind_rows(combined))


##-----------------------------------------------------------------------------
vd2 <- combineExpression(combined, vd2, proportion=F,cloneSize=c(
  Single=1, Small=5, Medium=25, Large=100, Hyper=100
))


##-----------------------------------------------------------------------------
trgv <- combineTCR(contig,  removeNA=F, filterMulti=T)

names(trgv) <- names(contig)

for (i in seq_along(trgv)) {
  trgv[[i]]$barcode <- paste0(names(trgv[i]), "_VD2_", trgv[[i]]$barcode, "-1")
}

trgv <- bind_rows(trgv)%>%
  separate(TCR1, into=c("TRGV", "TRGJ", "TRGC"), sep="\\.")%>%
  mutate(TRGV=factor(TRGV, levels=c(
    paste0("TRGV", 2:5), "TRGV5P", paste0("TRGV", 7:11)
  )))%>%
  select(barcode, TRGV)


##-----------------------------------------------------------------------------
vd2@meta.data <- vd2@meta.data%>%
  rownames_to_column("barcode")%>%
  left_join(trgv, by="barcode")%>%
  column_to_rownames("barcode")


##-----------------------------------------------------------------------------
rm(i, combined, contig, trgv)


# Saving Objects
##-----------------------------------------------------------------------------
vd1$cluster <- paste0("C", vd1$seurat_clusters)
vd2$cluster <- paste0("C", vd2$seurat_clusters)

saveRDS(vd1, "../../data_processed/10X/Vd1.Rds")
saveRDS(vd2, "../../data_processed/10X/Vd2.Rds")


# Session Info
##-----------------------------------------------------------------------------
sessionInfo()

