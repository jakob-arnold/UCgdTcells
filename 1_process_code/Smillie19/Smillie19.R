---
title: "Smillie 2019: Processing"
output: github_document
---




##-----------------------------------------------------------------------------
library(Seurat)
library(scCustomize)
library(harmony)
library(Matrix)
library(tidyverse)
library(data.table)


# Data
##-----------------------------------------------------------------------------
# Meta df
meta <- fread("../../data_raw/Smillie2019/meta.data.txt")
meta <- meta[-1]

# Immune Cells
imm <- readMM("../../data_raw/Smillie2019/gene_sorted-Imm.matrix.mtx")
rownames(imm) <- readLines("../../data_raw/Smillie2019/Imm.genes.tsv")
colnames(imm) <- readLines("../../data_raw/Smillie2019/Imm.barcodes2.tsv")

imm <- CreateSeuratObject(imm)

imm@meta.data <- imm@meta.data%>%
  rownames_to_column("NAME")%>%
  left_join(meta, by="NAME")%>%
  column_to_rownames("NAME")

# Epithelial Cells
ep <- readMM("../../data_raw/Smillie2019/gene_sorted-Epi.matrix.mtx")
rownames(ep) <- readLines("../../data_raw/Smillie2019/Epi.genes.tsv")
colnames(ep) <- readLines("../../data_raw/Smillie2019/Epi.barcodes2.tsv")

ep <- CreateSeuratObject(ep)

ep@meta.data <- ep@meta.data%>%
  rownames_to_column("NAME")%>%
  left_join(meta, by="NAME")%>%
  column_to_rownames("NAME")

# Stromal Cells
fb <- readMM("../../data_raw/Smillie2019/gene_sorted-Fib.matrix.mtx")
rownames(fb) <- readLines("../../data_raw/Smillie2019/Fib.genes.tsv")
colnames(fb) <- readLines("../../data_raw/Smillie2019/Fib.barcodes2.tsv")

fb <- CreateSeuratObject(fb)

fb@meta.data <- fb@meta.data%>%
  rownames_to_column("NAME")%>%
  left_join(meta, by="NAME")%>%
  column_to_rownames("NAME")


##-----------------------------------------------------------------------------
ep$compartment <- "Epithelial"
imm$compartment <- "Immune"
fb$compartment <- "Stromal"


##-----------------------------------------------------------------------------
genes <- Reduce(intersect, list(
  rownames(ep),
  rownames(imm),
  rownames(fb)
  )
  )

all <- merge(
  ep[genes,],y=c(imm[genes,], fb[genes,])
  )


##-----------------------------------------------------------------------------
all <- all%>%
  JoinLayers()

rm(imm, ep, fb, meta, genes)
gc()


# QC
##-----------------------------------------------------------------------------
all$percent.mt <- PercentageFeatureSet(all, pattern="^MT-")


##-----------------------------------------------------------------------------
VlnPlot_scCustom(
  all,
  features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by="compartment",
  pt.size=0
)


##-----------------------------------------------------------------------------
all <- subset(
  all, subset=
    (compartment=="Epithelial"&percent.mt<75)|
    (compartment!="Epithelial"&percent.mt<25)&
    nFeature_RNA>500&
    nFeature_RNA<5000&
    nCount_RNA<40000
)

gc()


# Processing
##-----------------------------------------------------------------------------
all <- all%>%
  NormalizeData()%>%
  FindVariableFeatures()%>%
  ScaleData()%>%
  FindVariableFeatures()%>%
  RunPCA()%>%
  RunUMAP(dims=1:30)


##-----------------------------------------------------------------------------
pb <- all%>%
  subset(subset=Health!="Non-inflamed")%>%
  AggregateExpression(
    group.by=c("Health", "compartment", "Subject"),
    return.seurat=T
    )


# Save Object
##-----------------------------------------------------------------------------
saveRDS(pb, "../../data_processed/Smillie19/PB.Rds")





