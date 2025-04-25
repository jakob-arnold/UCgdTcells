---
title: "Integrated Datasets (Li 2021, Smillie 2019)"
output: github_document
---




##-----------------------------------------------------------------------------
library(Seurat)
library(SCpubr)
library(scCustomize)
library(harmony)
library(Matrix)
library(tidyverse)
library(data.table)
library(readxl)
library(biomaRt)
library(clusterProfiler)
library(ReactomePA)


# Li 2021

Reading in Data
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


## Meta Data
##-----------------------------------------------------------------------------
data_li_raw@meta.data <- data_li_raw@meta.data %>%
  separate(orig.ident, into = c("disease", "ID"), sep = "_") %>%
  unite(diseaseID ,disease, ID, sep = "", remove = F)


## QC
##-----------------------------------------------------------------------------
data_li_raw[["percent.mt"]] <- PercentageFeatureSet(data_li_raw, pattern = "^MT-")

VlnPlot(data_li_raw, pt.size=0, ncol = 3,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))&NoLegend()


## Subset high quality cells
##-----------------------------------------------------------------------------
data_li <- subset(data_li_raw, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                nCount_RNA < 30000 & percent.mt < 50)


## Remove Genes for Clustering
##-----------------------------------------------------------------------------
non.mt.genes <- rownames(data_li[["RNA"]]$counts[grep("^(MT-)",
                    rownames(data_li[["RNA"]]$counts), invert = T),])

non.ribo.genes <- rownames(data_li[["RNA"]]$counts[grep("^(RPL\\d|RPS\\d|RP\\d|IGH|IGL|IGK)", rownames(data_li[["RNA"]]$counts), invert = T),])

data_li <- subset(data_li, features = non.mt.genes)
data_li <- subset(data_li, features = non.ribo.genes)

rm(list=c("non.mt.genes", "non.ribo.genes"))


## Clustering
##-----------------------------------------------------------------------------
data_li <- data_li%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA()%>%
  FindNeighbors(dims=1:30)%>%FindClusters()%>%RunUMAP(dims=1:30)


##-----------------------------------------------------------------------------
do_DimPlot(data_li, group.by = "seurat_clusters", label=T)&NoLegend()
FeaturePlot_scCustom(data_li, features = c("EPCAM", "PTPRC", "CD3E", "TRDC"))&
  NoLegend()&NoAxes()


## Subsetting EPCAM+ & gd Tcells
##-----------------------------------------------------------------------------
ep_li <- subset(data_li, subset=seurat_clusters %in% c(
  "0", "10", "15", "17", "18", "22", "26")
  )

gd_li <- subset(data_li, subset=TRDC>1 & seurat_clusters %in% c(
  "3", "5", "7", "8", "9", "12", "14")
  )


##-----------------------------------------------------------------------------
rm(list=c("data_li_raw", "data_li"))



# Smillie 2019

## Reading in Data
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


## QC & Processing
##-----------------------------------------------------------------------------
#Epithelial Cells
ep_sm[["percent.mt"]] <- PercentageFeatureSet(ep_sm, pattern = "^MT-")

VlnPlot(ep_sm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)&NoLegend()


##-----------------------------------------------------------------------------
ep_sm <- subset(ep_sm, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & 
               nCount_RNA < 30000 & percent.mt < 50)


##-----------------------------------------------------------------------------
# Immune Cells
imm[["percent.mt"]] <- PercentageFeatureSet(imm, pattern = "^MT-")

VlnPlot(imm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0)&NoLegend()


##-----------------------------------------------------------------------------
imm <- subset(imm, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & 
               nCount_RNA < 25000 & percent.mt < 15)


##-----------------------------------------------------------------------------
imm <- imm%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA()%>%
  FindNeighbors(dims=1:30)%>%FindClusters(resolution =1.5)%>%RunUMAP(dims=1:30)


##-----------------------------------------------------------------------------
do_DimPlot(imm, group.by="seurat_clusters", label=T, pt.size=.1)&NoLegend()
FeaturePlot_scCustom(imm, features=c("CD3E", "TRDC"))&NoAxes()&NoLegend()


##-----------------------------------------------------------------------------
gd_sm <- subset(imm, subset=TRDC>1 & seurat_clusters %in% c(
  "2", "18", "37", "44", "23", "36", "0", "46", "16", "5", "6", "7", "30", "29"
))


##-----------------------------------------------------------------------------
rm("imm")



# Combining and Integrating Datasets

## gd T cells
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


### Integtration
##-----------------------------------------------------------------------------
# Downsampling
gd_sm <- SetIdent(gd_sm, value="all")
gd_sm <- subset(gd_sm, downsample=844)

gd_li@assays$SCT <- NULL
DefaultAssay(gd_li) <- "RNA"
gd_sm@assays$SCT <- NULL
DefaultAssay(gd_sm) <- "RNA"


# Only use genes that appear in both datasets
genes <- intersect(rownames(gd_li), rownames(gd_sm))

# Merge
gd <- merge(gd_li[genes,], gd_sm[genes,])%>%JoinLayers()

# rm(list=c("genes", "gd_li", "gd_sm"))


##-----------------------------------------------------------------------------
gd <- gd%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%
  RunPCA()%>%FindNeighbors(dims=1:30)%>%FindClusters()%>%RunUMAP(dims=1:30)


##-----------------------------------------------------------------------------
do_DimPlot(gd, group.by = "cohort")


##-----------------------------------------------------------------------------
# Integration
gd <- RunHarmony(gd, "cohort")


##-----------------------------------------------------------------------------
ElbowPlot(gd, 50, reduction = "harmony")


##-----------------------------------------------------------------------------
gd <- gd%>%FindNeighbors(reduction="harmony", dims=1:20)%>%
  FindClusters(resolution=.25)%>%RunUMAP(reduction="harmony", dims=1:20)


##-----------------------------------------------------------------------------
gd@meta.data <- gd@meta.data%>%
  mutate(cluster=paste0("C", seurat_clusters))


##-----------------------------------------------------------------------------
do_DimPlot(gd, group.by="cohort")


### Module Scores
##-----------------------------------------------------------------------------
# Tissue Residency
TRM_genes <- scan("genesets/TRM.txt", character(), quote="")
gd <- AddModuleScore(gd, features=list(TRM_genes), name="TRM")

# Cytokine/Cytotoxic
cytokine_genes <- scan("genesets/cytokine.txt", character(), quote="")
cytotoxic_genes <- scan("genesets/cytotoxic.txt", character(), quote="")

cytokine_cytotoxic_genes <- c(cytokine_genes, cytotoxic_genes)
cytokine_cytotoxic_genes <- unique(cytokine_cytotoxic_genes)

gd <- AddModuleScore(gd, features=list(cytokine_cytotoxic_genes), name="cyto")

# FLASHseq gd T Subsets
Vg4Vd1 <- read.csv("genesets/Vg4Vd1_genes.csv")%>%
  filter(!gene %in% c("TRGV4", "TRDV1"))%>%pull(gene)

gd <- AddModuleScore(gd, features = list(Vg4Vd1),name="Vg4Vd1_")

Vg9Vd2 <- read.csv("genesets/Vg9Vd2_genes.csv")%>%
  filter(!gene %in% c("TRGV9", "TRDV2"))%>%pull(gene)

gd <- AddModuleScore(gd, features = list(Vg9Vd2), name="Vg9Vd2_")

rm(list=c("cytokine_genes", "cytotoxic_genes","cytokine_cytotoxic_genes",
          "TRM_genes", "Vg9Vd2", "Vg4Vd1", "stem_like_genes")
   )

Stem-Like/Effector Signatures
##-----------------------------------------------------------------------------
# Li 2024: Stem-like T cells are associated with the pathogenesis of ulcerative 
# colitis in humans
li_stem <- read_csv("../../genesets/li24_stem_marker.csv")%>%
  pull(gene)%>%UpdateSymbolList()

# Function for converting mouse to human gene symbols
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl", 
                  host="https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                  host="https://dec2021.archive.ensembl.org/")

mouse_to_human <- function(genes){
  human_genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                   values = genes, mart = mouse, attributesL = c("hgnc_symbol"), 
                   martL = human, uniqueRows=T)
  human_genes <- unique(human_genes[, 2])
  
  return(human_genes)
}

# Siddiqui 2019: Intratumoral Tcf1+PD-1+CD8+ T Cells with Stem-like Properties 
# Promote Tumor Control in Response to Vaccination and Checkpoint Blockade Immunotherapy
siddiqui19 <- read_xlsx("../../genesets/siddiqui19_stem_vs_effector.xlsx")

siddiqui19_stem <- siddiqui19%>%filter(log2_FC >= 1)%>%pull(gene)%>%
  mouse_to_human()

siddiqui19_effector <- siddiqui19%>%filter(log2_FC <= -1)%>%pull(gene)%>%
  mouse_to_human()

# Im 2016: Defining CD8+ T cells that provide the proliferative burst after PD-1 therapy
im16 <- read_tsv("../../genesets/im16_stem_vs_effector.tsv")%>%
  filter(adj.P.Val < 0.05)

im16_stem <- im16%>%filter(logFC >= 1)%>%pull(Gene.symbol)%>%
  mouse_to_human()

im16_effector <- im16%>%filter(logFC <= -1)%>%pull(Gene.symbol)%>%
  mouse_to_human()

# Wu 2016: The TCF1-Bcl6 axis counteracts type I interferon to repress
# exhaustion and maintain T cell stemness
wu16 <- read_tsv("../../genesets/wu16_stem_vs_effector.tsv")%>%
  filter(adj.P.Val < 0.05)%>%
  mutate(gene=str_extract(gene_assignment, "(?<=//\\s)[^/]+(?=\\s//)"))

wu16_stem <- wu16%>%filter(logFC >= 1)%>%pull(gene)%>%
  mouse_to_human()

wu16_effector <- wu16%>%filter(logFC <= -1)%>%pull(gene)%>%
  mouse_to_human()

# Utzschneider 2016: T Cell Factor 1-Expressing Memory-like CD8+ T Cells
# Sustain the Immune Response to Chronic Viral Infections
utz16 <- read_xlsx("../../genesets/utzschneider16_stem_vs_effector.xlsx")

utz16_stem <- utz16%>%filter(logFC >= 1)%>%pull(Gene_Symbol)%>%
  mouse_to_human()

utz16_effector <- utz16%>%filter(logFC <= -1)%>%pull(Gene_Symbol)%>%
  mouse_to_human()

##-----------------------------------------------------------------------------
# Module Score Function
scoring <- function(object, genes, name){
  object <- AddModuleScore(object, features=list(genes), name=name)
  
  return(object)
}

scoring_list <- list(
  siddiqui19_stem = siddiqui19_stem, siddiqui19_effector = siddiqui19_effector,
  im16_stem = im16_stem, im16_effector = im16_effector,
  wu16_stem = wu16_stem, wu16_effector = wu16_effector,
  utz16_stem = utz16_stem, utz16_effector = utz16_effector,
  li_stem = li_stem
)

for (score in names(scoring_list)){
  gd <- scoring(gd, scoring_list[[score]], score)
}

rm(list=c(grep(ls(), pattern="^im16|^li|^sidd|^utz16|^wu16", value=T), "genes", "score",
          "scoring_list" ,"mouse", "human", "mouse_to_human"
          ))


## EPCAM+ Cells
##-----------------------------------------------------------------------------
# Change meta data column names to match between datasets
ep_li@meta.data <- ep_li@meta.data%>%
  mutate(cohort="Li 2021")

ep_sm@meta.data <- ep_sm@meta.data%>%
  rename(disease=Health, ID=Subject)%>%
  mutate(disease=case_when(
    disease=="Healthy" ~ "HD",
    disease=="Inflamed" ~ "UC",
    disease=="Non-inflamed" ~ "SC"
    ),
    cohort="Smillie 2019")%>%
  unite(diseaseID, disease, ID, sep="", remove=F)

ep_sm <- subset(ep_sm, subset=disease!="SC")


### Integration
##-----------------------------------------------------------------------------
# Only use genes that appear in both datasets
genes <- intersect(rownames(ep_li), rownames(ep_sm))

ep <- merge(ep_li[genes,], ep_sm[genes,])%>%JoinLayers()

rm(genes, ep_li, ep_sm)


##-----------------------------------------------------------------------------
ep <- ep%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%
  RunPCA()%>%FindNeighbors(dims=1:30)%>%FindClusters()%>%RunUMAP(dims=1:30)


##-----------------------------------------------------------------------------
do_DimPlot(ep, group.by="cohort", pt.size=.2)


##-----------------------------------------------------------------------------
# Integration
ep <- RunHarmony(ep, "cohort")


##-----------------------------------------------------------------------------
ElbowPlot(ep, reduction = "harmony", 50)


##-----------------------------------------------------------------------------
ep <- ep%>%FindNeighbors(reduction="harmony", dims=1:40)%>%
  FindClusters()%>%RunUMAP(reduction="harmony", dims=1:40)


##-----------------------------------------------------------------------------
do_DimPlot(ep, group.by="cohort", pt.size=.2)



##-----------------------------------------------------------------------------
# Cell Numbers of CellTypes Per Cluster (annotation derived from
# Smillie datasets) to perform label transfer (enterocytes)

ep@meta.data <- ep@meta.data%>%
  mutate(enterocytes=case_when(
    Cluster %in% c("Best4+ Enterocytes", "Enterocyte Progenitors", "Enterocytes",
                   "Immature Enterocytes 1", "Immature Enterocytes 2") ~ "enterocytes",
    is.na(Cluster) ~ "Li21",
    .default = "other"
  ))

ep@meta.data%>%ggplot(aes(x=seurat_clusters, fill=enterocytes))+
  geom_bar(position = "fill")+
  theme_classic()


##-----------------------------------------------------------------------------
# Subsetting Enterocytes based on Clusters and downsampling
entero <- subset(ep, subset = seurat_clusters %in% c(
  "0", "1", "7", "11", "16", "17", "21", "25"
))


##-----------------------------------------------------------------------------
summary(as.factor(entero@meta.data$diseaseID))



##-----------------------------------------------------------------------------
# Downsampling
entero <- SetIdent(entero, value="diseaseID")
entero <- subset(entero, downsample=500)


##-----------------------------------------------------------------------------
entero <- entero%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA()%>%
  RunHarmony("cohort")

##-----------------------------------------------------------------------------
ElbowPlot(entero, reduction="harmony", 50)

##-----------------------------------------------------------------------------
entero <- entero%>%FindNeighbors(reduction="harmony",dims=1:20)%>%FindClusters()%>%
  RunUMAP(reduction="harmony", dims=1:20)


##-----------------------------------------------------------------------------
do_DimPlot(entero, group.by="disease")
do_DimPlot(entero, group.by="cohort")

## Module Scores
##-----------------------------------------------------------------------------
# Cholesterol
chol_genes <- read_tsv("../../genesets/cholesterol.tsv", col_names=F)%>%
  filter(X1=="GENE_SYMBOLS")%>%pull(X2)%>%strsplit(",")

entero <- AddModuleScore(entero, features=chol_genes, name="chol")

# Mevalonate
meva_genes <- read_tsv("../../genesets/mevalonate.tsv", col_names=F)%>%
  filter(X1=="GENE_SYMBOLS")%>%pull(X2)%>%strsplit(",")
entero <- AddModuleScore(entero, features=meva_genes, name="meva")

rm(chol_genes, meva_genes)


## Gene Set Enrichment (GSEA)
##-----------------------------------------------------------------------------
Idents(entero) <- "disease"

diffexpgenes <- FindMarkers(entero, ident.1="UC", ident.2="HD", min.pct=0, logfc.threshold=0)

eg <- bitr(rownames(diffexpgenes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg <- subset(eg,!duplicated(eg$SYMBOL))
row.names(eg) <- eg$SYMBOL
eg$SYMBOL <- NULL
diffexpgenes$SYMBOL <- rownames(diffexpgenes)
diffexpgenes <- subset(diffexpgenes,!duplicated(diffexpgenes$SYMBOL))
row.names(diffexpgenes) <- diffexpgenes$SYMBOL
diffexpgenes <- merge(diffexpgenes,eg, by="row.names")
fc <- diffexpgenes$avg_log2FC
names(fc) <- diffexpgenes$ENTREZID
fc <- sort(fc, decreasing=T)
y <- gsePathway(fc, verbose=F, organism="human", seed=1337, pvalueCutoff=1)
res <- as.data.frame(y)



##-----------------------------------------------------------------------------
write_csv(res, "../../data_processed/Li21_Smillie19_Integrated/GSEA.csv")
saveRDS(y, "../../data_processed/Li21_Smillie19_Integrated/GSEA.Rds")

diffexpgenes%>%
  dplyr::filter(p_val_adj<0.05)%>%
  write.csv("Li21Sm19_Entero_UCvHD_Marker_onlySignificant.csv")

rm(fc, eg, res, y, diffexpgenes)



# Save Seurat Object as Rds files
##-----------------------------------------------------------------------------
saveRDS(gd, "../../data_processed/Li21_Smillie19_Integrated/gdTcells.Rds")
saveRDS(entero, "../../data_processed/Li21_Smillie19_Integrated/enterocytes.Rds")
saveRDS(ep,  "../../data_processed/Li21_Smillie19_Integrated/EPCAM.Rds")


# Session Info
##-----------------------------------------------------------------------------
sessionInfo()





