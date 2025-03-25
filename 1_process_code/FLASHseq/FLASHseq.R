---
title: "FLASH-seq of gd T Cells: Processing"
output: github_document
---



##-----------------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(Matrix)
library(Seurat)
library(scCustomize)
library(SCpubr)
library(ggpubr)
library(biomaRt)
library(readxl)
library(SCENIC)
library(AUCell)


# Processing Data

## Reading Count Matrix

##-----------------------------------------------------------------------------
mtx <- readMM("../../data_raw/FLASHseq/counts.mtx")
cells <- read_tsv("../../data_raw/FLASHseq/cells.tsv", col_names = F)$X1
genes <- read_tsv("../../data_raw/FLASHseq/genes.tsv", col_names = F)$X1

colnames(mtx) <- cells
rownames(mtx) <- genes

gd_raw <- CreateSeuratObject(mtx)

rm(mtx, cells, genes)


## Creating Lookup Table for Patient Data

##-----------------------------------------------------------------------------
# Function to print complete rows (letters A-P = 1-16)
rws <- function(numb) {
  sequence <- paste0(rep(LETTERS[numb], each=24), 1:24)
  return(sequence)
}

lookup_table <- data.frame(
  plate = rep(c("E64","E65", "E66", "E67", "E70"), each = 16*24),
  well = rep(paste0(rep(LETTERS[1:16], each = 24), rep(1:24, times=16)), times=5)
            )%>%
  mutate(ID = case_when(
    plate == "E66" & well %in% c(rws(1:3), paste0("D", 17:24)) ~ "6588",
    plate == "E66" & well %in% c(rws(5:16), paste0("D", 1:16)) ~ "8333",
    plate == "E64" & well %in% c(rws(1:7), paste0("H", 22:24)) ~ "6809",
    plate == "E64" & well %in% c(rws(9:16), paste0("H", 1:21)) ~ "7802",
    plate == "E65" ~ "7130",
    plate == "E67" & well %in% c(rws(1:5), paste0("F", 14:24)) ~ "7802",
    plate == "E67" & well %in% c(rws(7:16), paste0("F", 1:13)) ~ "7130",
    plate == "E70" & well %in% c(rws(1:7), paste0("H", 15:24)) ~ "8626",
    .default = "3600"
        ))%>%
 unite("PlateWell", plate:well)


## Adding Patient Meta Data

##-----------------------------------------------------------------------------
gd_raw@meta.data <- gd_raw@meta.data%>%
  mutate(PlateWell = rownames(gd_raw@meta.data))%>%
  left_join(lookup_table, by="PlateWell")%>%
  column_to_rownames("PlateWell")

gd_raw@meta.data <- gd_raw@meta.data%>%
  mutate(disease = case_when(
    ID %in% c("7130", "7802") ~ "control",
    .default = "UC"))%>%
  mutate(treatment = case_when(
    ID %in% c("8626", "6809", "6588") ~ "Mesalazine",
    .default = "no treatment"))%>%
  unite("disease_ID", c("disease", "treatment", "ID"), sep = "_", remove = F)

rm(rws, lookup_table)


## QC

##-----------------------------------------------------------------------------
gd_raw <- SetIdent(gd_raw, value="all")
gd_raw[["percent.mt"]] <- PercentageFeatureSet(gd_raw, pattern = "^MT-")

p <- VlnPlot_scCustom(gd_raw, add.noise=F, pt.size=0,
                      features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))&
  theme(axis.title=element_blank())

p[[1]] <- p[[1]]+geom_hline(yintercept=c(500, 4000), linetype="dashed", color="red")+
  ggtitle("Number of Genes")
p[[2]] <- p[[2]]+geom_hline(yintercept=c(50000), linetype="dashed", color="red")+
  ggtitle("UMI counts")
p[[3]] <- p[[3]]+geom_hline(yintercept=c(75), linetype="dashed", color="red")+
  ggtitle("% Mitochondrial\nGenes")
p

# Save QC Plots 
ggsave("../../2_figures_code/figures/FLASHseq_VlnPlots_QC.pdf", width=6.5, height=3.5)
ggsave("../../2_figures_code/figures/FLASHseq_VlnPlots_QC.png", dpi=600, width=6.5, height=3.5)

rm(p)


## Subsetting High Quality Cells

##-----------------------------------------------------------------------------
gd <- subset(gd_raw, subset = percent.mt < 75 & nFeature_RNA > 500 &
             nFeature_RNA < 4000 & nCount_RNA > 50000)

# Removing ribosomal and mitochondrial genes
non.R_M.genes <- rownames(gd[["RNA"]]$counts[grep("^(RPL\\d|RPS\\d|MT-)",
                          rownames(gd[["RNA"]]$counts), invert = T),])
gd <- subset(gd, features = non.R_M.genes)

rm(gd_raw, non.R_M.genes)


## Processing

##-----------------------------------------------------------------------------
gd <- gd%>%NormalizeData()%>%FindVariableFeatures()%>%ScaleData()%>%
  RunPCA(verbose=F)%>%RunUMAP(dims=1:30)%>%FindNeighbors(dims=1:30)%>%FindClusters()


##-----------------------------------------------------------------------------
FeaturePlot_scCustom(gd, features = c("MS4A1", "CD19", "MS4A7", "CD14", "CD4",
                                      "TRDC"))&NoAxes()&NoLegend()
do_DimPlot(gd, label=T)&NoAxes()&NoLegend()


## Remove non gd Cells; Recluster

##-----------------------------------------------------------------------------
gd <- subset(gd, idents = c("4", "5"), invert = T)

gd <- gd%>%FindVariableFeatures()%>%ScaleData()%>%RunPCA(verbose=F)%>%
  RunUMAP(dims=1:30)%>%FindNeighbors(dims=1:30)%>%FindClusters()

gd@meta.data <- gd@meta.data%>%
  mutate(cluster = paste0("C", seurat_clusters))


##-----------------------------------------------------------------------------
do_DimPlot(gd, pt.size=2, label=T)&NoAxes()&NoLegend()


# Gene Set Enrichment (GSEA)

##-----------------------------------------------------------------------------
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)

Idents(gd) <- "disease"

diffexpgenes <- FindMarkers(gd, ident.1 = "UC", ident.2 = "control", min.pct	= 0, logfc.threshold=0)

eg <- bitr(rownames(diffexpgenes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg <- subset(eg,!duplicated(eg$SYMBOL))
row.names(eg) <- eg$SYMBOL
eg$SYMBOL <- NULL
diffexpgenes$SYMBOL <- rownames(diffexpgenes)
diffexpgenes <- subset(diffexpgenes,!duplicated(diffexpgenes$SYMBOL))
row.names(diffexpgenes) <- diffexpgenes$SYMBOL
diffexpgenes <- merge(diffexpgenes,eg, by = "row.names")
fc <- diffexpgenes$avg_log2FC
names(fc) <- diffexpgenes$ENTREZID
fc <- sort(fc, decreasing = T)
y <- gsePathway(fc, verbose=T, organism = "human", seed=1337, pvalueCutoff=1)
res <- as.data.frame(y)



##-----------------------------------------------------------------------------
write_csv(res, "../../data_processed/FLASHseq/GSEA.csv")
saveRDS(y, "../../data_processed/FLASHseq/GSEA_res.Rds")

rm(fc, eg, res, y, diffexpgenes)


# TCR Data

##-----------------------------------------------------------------------------
library(scRepertoire)

# Load Contig
contig <- read_tsv("../../data_raw/FLASHseq/clones.tsv")
contig <- contig%>%
  unite(tagValueCELL, tagValueCELL0ROW, tagValueCELL0COL, sep="_")
contig <- list(contig)

contig <- loadContigs(input=contig, format="MiXCR")

# Remove allele variant calls + scoring (*00(###))
contig <- lapply(contig, function(df) {
  df%>%
    mutate_at(paste0(c("v", "d", "j", "c"), "_gene"), ~str_replace(., "\\*.*", ""))
})

contig[[1]] <- contig[[1]]%>%
  filter(chain %in% c("TRG", "TRD"))


##-----------------------------------------------------------------------------
combined <- combineTCR(contig, filterMulti=T, removeNA=F)
combined[[1]]$all <- "all"

# Combine with SeuratObject (obj)
gd <- combineExpression(combined, gd, group.by = "all", prop=F,
                        cloneSize = c(
   Single=1, Rare=5, Medium=10, Large=25, Hyper=25))

# Remove Clones w/ NA values
gd@meta.data <- gd@meta.data%>%
  separate(CTstrict, into=c("trg_strict", "trd_strict"), remove=F,
           sep="_")%>%
   mutate(CTstrict=case_when(
     str_detect(CTstrict, "NA") ~ NA_character_,
     .default=CTstrict),
     trg_strict=case_when(
      str_detect(trg_strict, "NA") ~ NA_character_,
      .default=trg_strict),
     trd_strict=case_when(
       str_detect(trd_strict, "NA") ~ NA_character_,
       .default=trd_strict
     ))


##-----------------------------------------------------------------------------
# Add TRGV/TRDV Data to SeuratObject
combined <- combined[[1]]

combined <- combined%>%
  separate(TCR1, into=c("TRGV", "TRGJ", "TRGC"), sep="\\.", remove=F)%>%
  separate(TCR2, into=c("TRDV", "TRDD", "TRDJ", "TRDC"), sep="\\.", remove=F)

tcr_meta <- combined%>%
  dplyr::select(barcode, TRGV, TRDV)%>%
  mutate(
    TRGV=factor(TRGV, levels=c(paste0("TRGV", 2:5), paste0("TRGV", 7:11))),
    TRDV=factor(TRDV, levels=c("TRDV1", "TRDV2", "TRDV3", "TRAV29/DV5", 
                               "TRAV36/DV7", "TRAV38-2/DV8"))
        )

gd@meta.data <- gd@meta.data%>%
  rownames_to_column("barcode")%>%
  left_join(tcr_meta, by="barcode")%>%
  column_to_rownames("barcode")%>%
  mutate(TRGV_TRDV=case_when(
    (!is.na(TRGV) & !is.na(TRDV)) ~ paste(TRGV, TRDV, sep="_"),
    .default = NA_character_
  ))

rm(tcr_meta)


##-----------------------------------------------------------------------------
# C2/C3 Clonal Overlaps

# TRG_clonotypes <- gd@meta.data %>%
#   filter(!is.na(trg_strict) & cluster %in% c("C2", "C3"))%>%
#   group_by(trg_strict) %>%
#   filter(n_distinct(cluster) > 1) %>%
#   arrange(CTstrict)
# 
# TRG_clonotypes <- unique(TRG_clonotypes$trg_strict)
# TRG_clonotypes
#  
# TRD_clonotypes <- gd@meta.data %>%
#   filter(!is.na(trd_strict) & cluster %in% c("C2", "C3"))%>%
#   group_by(trd_strict) %>%
#   filter(n_distinct(cluster) > 1) %>%
#   arrange(CTstrict)
# 
# TRD_clonotypes <- unique(TRD_clonotypes$trd_strict)
# TRD_clonotypes
# 
# summary(is.na(gd@meta.data$trd_strict)&is.na(gd@meta.data$trg_strict))


# Module Scores

## Cytokine/Cytotoxicity/IFN-response Module Scores

##-----------------------------------------------------------------------------
# Tissue resident memory T cells
TRM_genes <- scan("../../genesets/TRM.txt", character(), quote = "")
gd <- AddModuleScore(gd, features = list(TRM_genes), name = "TRM_signature")

# Cytokine/Cytotoxic scoring
cytokine_genes <- scan("../../genesets/cytokine.txt", character(), quote = "")
cytotoxic_genes <- scan("../../genesets/cytotoxic.txt", character(), quote = "")
cytokine_cytotoxic_genes <- c(cytokine_genes, cytotoxic_genes)
cytokine_cytotoxic_genes <- base::unique(cytokine_cytotoxic_genes)
gd <- AddModuleScore(gd, features = list(cytokine_cytotoxic_genes), 
        name = "cytokine_cytotoxic_signature")

# IFN response
IFN_genes <- scan("../../genesets/IFN_response.txt", character(), quote = "")
gd <- AddModuleScore(gd, features = list(IFN_genes), name = "IFN_signature")

rm(TRM_genes, cytokine_genes, cytotoxic_genes, cytokine_cytotoxic_genes, IFN_genes)


## Stem-Like Signature Genesets

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

im16_effector <- im16%0.1534872>%filter(logFC <= -1)%>%pull(Gene.symbol)%>%
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

rm(list=setdiff(ls(), c("combined", "contig", "gd")))


## Genesets for Fig 3

##-----------------------------------------------------------------------------
gd <- SetIdent(gd, value="TRGV_TRDV")

TRGV4_TRDV1_markers <- FindMarkers(gd, ident.1="TRGV4_TRDV1", only.pos=T,
  min.pct=.2, logfc.threshold=1)%>%filter(p_val_adj<0.05)%>%
  rownames_to_column("gene")
write.csv(TRGV4_TRDV1_markers, "../../genesets/Vg4Vd1_genes.csv")

TRGV9_TRDV2_markers <- FindMarkers(gd, ident.1="TRGV9_TRDV2", only.pos=T,
  min.pct=.2, logfc.threshold=1)%>%
  filter(p_val_adj < 0.05)%>%rownames_to_column("gene")
write.csv(TRGV9_TRDV2_markers, "Vg9Vd2_genes.csv")

rm(TRGV4_TRDV1_markers, TRGV9_TRDV2_markers)


# SCENIC

##-----------------------------------------------------------------------------
scenic <- read_csv("../../data_processed/FLASHseq/pyscenic_output.csv")
scenic <- as.matrix(scenic)
rownames(scenic) <- scenic[,1]
scenic <- scenic[, -1]
scenic <- t(scenic)
rownames(scenic) <- gsub(rownames(scenic), pattern="\\(\\+\\)", replacement="")

gd[["scenic"]] <- CreateAssayObject(scenic)

rm(scenic)


##-----------------------------------------------------------------------------
DefaultAssay(gd) <- "scenic"

gd <- gd%>%ScaleData()%>%RunPCA(features=rownames(gd), reduction.name="scenic_pca")

DefaultAssay(gd) <- "RNA"


##-----------------------------------------------------------------------------
ElbowPlot(gd, reduction="scenic_pca", 50)


##-----------------------------------------------------------------------------
gd <- RunUMAP(gd, reduction.name="scenic_umap", reduction="scenic_pca", dims=1:20)


# Saving Object

##-----------------------------------------------------------------------------
saveRDS(gd, "../../data_processed/FLASHseq/gd.Rds")


# Session Info

##-----------------------------------------------------------------------------
sessionInfo()

