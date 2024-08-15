#---
#title: "Integrated Datasets (Li 2021, Smillie 2019)"
#output: github_document
#---




##----------------------------------------------------------------
library(Seurat)
library(SCpubr)
library(scCustomize)
library(harmony)
library(Matrix)
library(tidyverse)
library(data.table)


## Li 2021

#Reading in Data
##-----------------------------------------------------------------------------
HD_1 <- Read10X("../../data_raw/Li2021/id6C/")
HD_2 <- Read10X("../../data_raw/Li2021/id7C/")
HD_3 <- Read10X("../../data_raw/Li2021/id8C/")
HD_4 <- Read10X("../../data_raw/Li2021/id9C/")

UC_1 <- Read10X("../../data_raw/Li2021/id1T/")
UC_2 <- Read10X("../../data_raw/Li2021/id2T/")
UC_3 <- Read10X("../../data_raw/Li2021/id3T/")
UC_4 <- Read10X("../../data_raw/Li2021/id4T/")
UC_5 <- Read10X("../../data_raw/Li2021/id5T/")

HD_1 <- CreateSeuratObject(HD_1, project = "HD_1")
HD_2 <- CreateSeuratObject(HD_2, project = "HD_2")
HD_3 <- CreateSeuratObject(HD_3, project = "HD_3")
HD_4 <- CreateSeuratObject(HD_4, project = "HD_4")

UC_1 <- CreateSeuratObject(UC_1, project = "UC_1")
UC_2 <- CreateSeuratObject(UC_3, project = "UC_2")
UC_3 <- CreateSeuratObject(UC_4, project = "UC_3")
UC_4 <- CreateSeuratObject(UC_5, project = "UC_4")
UC_5 <- CreateSeuratObject(UC_5, project = "UC_5")

data_li_raw <- merge(
  HD_1, y = mget(c(paste0("HD_", 2:4), paste0("UC_", 1:5))),
  add.cell.ids = c(paste0("HD_", 1:4), paste0("UC_", 1:5))
  )

rm(list=c(paste0("HD_", 1:4), paste0("UC_", 1:5)))


##-----------------------------------------------------------------------------
data_li_raw <- JoinLayers(data_li_raw)
data_li_raw <- SetIdent(data_li_raw, value="all")


### Meta Data
##-----------------------------------------------------------------------------
data_li_raw@meta.data <- data_li_raw@meta.data %>%
  separate(orig.ident, into = c("disease", "ID"), sep = "_") %>%
  unite(diseaseID ,disease, ID, sep = "", remove = F)


### QC
##----------------------------------------------------------------
data_li_raw[["percent.mt"]] <- PercentageFeatureSet(data_li_raw, pattern = "^MT-")

VlnPlot(data_li_raw, pt.size=0, ncol = 3,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))&NoLegend()


### Subset high quality cells
##-----------------------------------------------------------------------------
data_li <- subset(data_li_raw, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                nCount_RNA < 30000 & percent.mt < 50)


### Remove Genes for Clustering
##-----------------------------------------------------------------------------
non.mt.genes <- rownames(data_li[["RNA"]]$counts[grep("^(MT-)",
                    rownames(data_li[["RNA"]]$counts), invert = T),])

non.ribo.genes <- rownames(data_li[["RNA"]]$counts[grep("^(RPL\\d|RPS\\d|RP\\d|IGH|IGL|IGK)", rownames(data_li[["RNA"]]$counts), invert = T),])

data_li <- subset(data_li, features = non.mt.genes)
data_li <- subset(data_li, features = non.ribo.genes)

rm(list=c("non.mt.genes", "non.ribo.genes"))


### Clustering
##----------------------------------------------------------------
data_li <- data_li%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%
  RunPCA()%>%FindNeighbors(dims=1:30)%>%FindClusters()%>%RunUMAP(dims=1:30)


##----------------------------------------------------------------
do_DimPlot(data_li, group.by = "seurat_clusters", label=T)&NoLegend()
FeaturePlot_scCustom(data_li, features = c("EPCAM", "PTPRC", "CD3E", "TRDC"))&
  NoLegend()&NoAxes()


### Subsetting EPCAM+ & gd Tcells
##-----------------------------------------------------------------------------
ep_li <- subset(data_li, subset=seurat_clusters %in% c(
  "0", "10", "15", "17", "18", "22", "26"
))

gd_li <- subset(data_li, subset=TRDC>1 & seurat_clusters %in% c(
  "3", "5", "7", "8", "9", "12", "14"
))


##-----------------------------------------------------------------------------
rm(list=c("data_li_raw", "data_li"))



## Smillie 2019

### Reading in Data
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
ep_sm <- readMM("../../data_raw/Smillie2019/gene_sorted-Epi.matrix.mtx")
rownames(ep_sm) <- readLines("../../data_raw/Smillie2019/Epi.genes.tsv")
colnames(ep_sm) <- readLines("../../data_raw/Smillie2019/Epi.barcodes2.tsv")

ep_sm <- CreateSeuratObject(ep_sm)

ep_sm@meta.data <- ep_sm@meta.data%>%
  rownames_to_column("NAME")%>%
  left_join(meta, by="NAME")%>%
  column_to_rownames("NAME")

rm("meta")


### QC & Processing
##----------------------------------------------------------------
#Epithelial Cells
ep_sm[["percent.mt"]] <- PercentageFeatureSet(ep_sm, pattern = "^MT-")

VlnPlot(ep_sm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)&NoLegend()


##-----------------------------------------------------------------------------
ep_sm <- subset(ep_sm, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & 
               nCount_RNA < 30000 & percent.mt < 50)


##----------------------------------------------------------------
# Immune Cells
imm[["percent.mt"]] <- PercentageFeatureSet(imm, pattern = "^MT-")

VlnPlot(imm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)&NoLegend()


##-----------------------------------------------------------------------------
imm <- subset(imm, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
               nCount_RNA < 25000 & percent.mt < 15)


##----------------------------------------------------------------
imm <- imm%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%
  RunPCA()%>%FindNeighbors(dims=1:30)%>%
  FindClusters(resolution =1.5)%>%RunUMAP(dims=1:30)


##----------------------------------------------------------------
do_DimPlot(imm, group.by="seurat_clusters", label=T, pt.size=.1)&NoLegend()
FeaturePlot_scCustom(imm, features=c("CD3E", "TRDC"))&NoAxes()&NoLegend()


##-----------------------------------------------------------------------------
gd_sm <- subset(imm, subset=TRDC>1 & seurat_clusters %in% c(
  "2", "18", "37", "44", "23", "36", "0", "46", "16", "5", "6", "7", "30", "29"
))


##-----------------------------------------------------------------------------
rm("imm")



## Combining and Integrating Datasets

### gd T cells
##-----------------------------------------------------------------------------
# changing meta data column names to match between both datasets

gd_li@meta.data <- gd_li@meta.data%>%
  mutate(cohort="Li2021")

gd_sm@meta.data <- gd_sm@meta.data%>%
  rename(disease=Health, ID=Subject)%>%
  mutate(disease=case_when(
    disease=="Healthy" ~ "HD",
    disease=="Inflamed" ~ "UC",
    disease=="Non-inflamed" ~ "SC"
    ),
    cohort="Smillie2019")%>%
  unite(diseaseID, disease, ID, sep="", remove=F)

gd_sm <- subset(gd_sm, subset=disease!="SC")


##-----------------------------------------------------------------------------
# Downsample
gd_sm <- SetIdent(gd_sm, value="all")
gd_sm <- subset(gd_sm, downsample=844)


#### Integtration
##-----------------------------------------------------------------------------
# Only use genes that appear in both datasets
genes <- intersect(rownames(gd_li), rownames(gd_sm))

# Merge
gd <- merge(gd_li[genes,], gd_sm[genes,])%>%JoinLayers()

rm(list=c("genes", "gd_li", "gd_sm"))


##----------------------------------------------------------------
gd <- gd%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%
  RunPCA()%>%RunUMAP(dims=1:30)


##-----------------------------------------------------------------------------
do_DimPlot(gd, group.by = "cohort")


##----------------------------------------------------------------
# Integration
gd <- RunHarmony(gd, "cohort")


##-----------------------------------------------------------------------------
ElbowPlot(gd, 50, reduction = "harmony")


##----------------------------------------------------------------
gd <- gd%>%FindNeighbors(reduction="harmony", dims=1:20)%>%
  FindClusters(resolution = 0.25)%>%
  RunUMAP(reduction="harmony", dims=1:20)


##-----------------------------------------------------------------------------
gd@meta.data <- gd@meta.data%>%
  mutate(cluster=paste0("C", seurat_clusters))


##-----------------------------------------------------------------------------
do_DimPlot(gd, group.by = "cohort")



#### Module Scores
##-----------------------------------------------------------------------------
TRM_genes <- scan("genesets/TRM.txt", character(), quote = "")

gd <- AddModuleScore(gd, features = list(TRM_genes),
                     name = "TRM_signature")

cytokine_genes <- scan("genesets/cytokine.txt", character(),
                       quote = "")
cytotoxic_genes <- scan("genesets/cytotoxic.txt", character(),
                        quote = "")

cytokine_cytotoxic_genes <- c(cytokine_genes, cytotoxic_genes)
cytokine_cytotoxic_genes <- base::unique(cytokine_cytotoxic_genes)

gd <- AddModuleScore(gd, features = list(cytokine_cytotoxic_genes), 
        name = "cytokine_cytotoxic_signature")

rm(list=c("cytokine_genes", "cytotoxic_genes",
          "cytokine_cytotoxic_genes", "TRM_genes"))


### EPCAM+ Cells
##-----------------------------------------------------------------------------
# Change meta data column names to match between datasets
ep_li@meta.data <- ep_li@meta.data%>%
  mutate(cohort="Chinese")

ep_sm@meta.data <- ep_sm@meta.data%>%
  rename(disease=Health, ID=Subject)%>%
  mutate(disease=case_when(
    disease=="Healthy" ~ "HD",
    disease=="Inflamed" ~ "UC",
    disease=="Non-inflamed" ~ "SC"
    ),
    cohort="US")%>%
  unite(diseaseID, disease, ID, sep="", remove=F)

ep_sm <- subset(ep_sm, subset=disease!="SC")


#### Integration
##-----------------------------------------------------------------------------
# Only use genes that appear in both datasets
genes <- intersect(rownames(ep_li), rownames(ep_sm))

ep <- merge(ep_li[genes,], ep_sm[genes,])%>%JoinLayers()

rm(list=c("genes", "ep_li", "ep_sm"))


##----------------------------------------------------------------
ep <- ep%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%
  RunPCA()%>%RunUMAP(dims=1:30)


##-----------------------------------------------------------------------------
do_DimPlot(ep, group.by="cohort", pt.size=.2)


##----------------------------------------------------------------
# Integration
ep <- RunHarmony(ep, "cohort")


##-----------------------------------------------------------------------------
ElbowPlot(ep, reduction = "harmony", 50)


##----------------------------------------------------------------
ep <- ep%>%FindNeighbors(reduction="harmony", dims=1:40)%>%
  FindClusters()%>%RunUMAP(reduction="harmony", dims=1:40)


##-----------------------------------------------------------------------------
# Cell Numbers of CellTypes Per Cluster (annotation derived from
# Smillie datasets) to perform label transfer (enterocytes)
table(ep@meta.data$Cluster, ep@meta.data$seurat_clusters)


##-----------------------------------------------------------------------------
# Subsetting Enterocytes based on Clusters and downsampling
ep <- SetIdent(ep, value="diseaseID")

entero <- subset(ep, subset = seurat_clusters %in% c(
  "0", "1", "7", "11", "16", "17", "21", "25"), downsample=500)


##----------------------------------------------------------------
entero <- entero%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA()%>%
  RunHarmony("cohort")


##-----------------------------------------------------------------------------
ElbowPlot(entero, reduction="harmony", 50)


##----------------------------------------------------------------
entero <- entero%>%FindNeighbors(reduction="harmony",dims=1:20)%>%
  FindClusters()%>%RunUMAP(reduction="harmony", dims=1:20)


##-----------------------------------------------------------------------------
do_DimPlot(entero, group.by="disease")



## Save Seurat Object as Rds files
##-----------------------------------------------------------------------------
saveRDS(gd, "../../data_processed/Li21_Smillie19_integrated_gdTcells.Rds")
saveRDS(entero, "../../data_processed/Li21_Smillie19_integrated_enterocytes.Rds")


## Session Info
##-----------------------------------------------------------------------------
# utils:::print.sessionInfo(sessionInfo()[-10]) 
sessionInfo()



















