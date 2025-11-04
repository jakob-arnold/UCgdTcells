---
title: "Thomas et al. 2024: Processing scRNA-seq Data"
output: github_document
---



##-----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scCustomize)
library(SCpubr)
library(schard)
library(harmony)
library(ggpubr)


# Data
##-----------------------------------------------------------------------------
cd8 <- h5ad2seurat("../../data_raw/Thomas24/cd8tcells_final.h5ad")
ep <-  h5ad2seurat("../../data_raw/Thomas24/epicolonic_final.h5ad")


##-----------------------------------------------------------------------------
ep
cd8

# EpCAM/Enterocytes
##-----------------------------------------------------------------------------
table(ep$Disease,
      ep$Treatment)


##-----------------------------------------------------------------------------
table(ep$Disease,
      ep$Remission_status)


##-----------------------------------------------------------------------------
table(ep$minor)


##-----------------------------------------------------------------------------
ep@meta.data <- ep@meta.data%>%
  unite("resp_tp", Remission_status, Treatment, remove=F)%>%
  unite("resp_tp_infl", Remission_status, Treatment, Inflammation, remove=F)


##-----------------------------------------------------------------------------
table(ep$resp_tp,
      ep$Inflammation)


##-----------------------------------------------------------------------------
table(ep$resp_tp_infl,
      ep$Disease)


##-----------------------------------------------------------------------------
entero <- subset(
  ep,
  subset=Disease=="UC" &
    minor=="Non_ileal_enterocyte" &
    resp_tp_infl %in% c("Non_Remission_Post_Inflamed",
                        "Non_Remission_Pre_Inflamed",
                        "Remission_Post_Non_Inflamed",
                        "Remission_Pre_Inflamed")
               )

gc()


##-----------------------------------------------------------------------------
table(entero$resp_tp_infl,
      entero$Disease)

## QC
##-----------------------------------------------------------------------------
entero$percent.mt <- PercentageFeatureSet(entero, pattern="^MT-")


##-----------------------------------------------------------------------------
genes_to_remove <- grep(rownames(entero), pattern="^(RPL\\d|RPS\\d|MT-)")

entero <- subset(entero, features=rownames(entero)[-genes_to_remove])


##-----------------------------------------------------------------------------
VlnPlot_scCustom(entero,
                 features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by="orig.ident",
                 pt.size=0)&
  geom_boxplot(width=.3, fill="white", outlier.shape=NA, coef=0)

##-----------------------------------------------------------------------------
entero <- subset(entero, subset=
                   nFeature_RNA>800 &
                   nCount_RNA<50000)


## Pseudo-Bulk
##-----------------------------------------------------------------------------
pb <- AggregateExpression(
  entero,
  group.by=c("Remission_status", "Treatment", "Patient"),
  return.seurat=T
)

##-----------------------------------------------------------------------------
table(pb$Remission_status,
      pb$Treatment)


##-----------------------------------------------------------------------------
table(pb$Patient)


##-----------------------------------------------------------------------------
# Remove Patients without paired data
pb <- subset(pb, Patient %in% c("UC16", "UC20"), invert=T)


##-----------------------------------------------------------------------------
table(pb$Patient)

##-----------------------------------------------------------------------------
pb@meta.data <- pb@meta.data%>%
  unite("resp_tp", Remission_status, Treatment, remove=F)%>%
  mutate(resp_tp=factor(resp_tp, levels=c("Remission_Pre",
                                          "Remission_Post",
                                          "Non-Remission_Pre",
                                          "Non-Remission_Post")))


# gd T
##-----------------------------------------------------------------------------
gd <- subset(cd8, TRDC>0 & Disease %in% c("Healthy", "UC"))


##-----------------------------------------------------------------------------
gd@meta.data <- gd@meta.data%>%
  unite("tp_responder", Treatment, Remission_status, remove=F)%>%
  unite("tp_responder_inflammation", tp_responder, Inflammation,remove=F)


##-----------------------------------------------------------------------------
table(gd$tp_responder,
      gd$Inflammation)


##-----------------------------------------------------------------------------
gd <- gd%>%
  subset(subset=tp_responder_inflammation %in% c(
    "Post_Non_Remission_Inflamed", "Post_Remission_Non_Inflamed",
    "Pre_Non_Remission_Inflamed", "Pre_Remission_Inflamed"
    
  )|Disease=="Healthy"
  )


##-----------------------------------------------------------------------------
gd@meta.data <- gd@meta.data%>%
  unite("res_tp", Remission_status, Treatment, remove=F)%>%
  mutate(res_tp=factor(res_tp, levels=c(
    "Non_Remission_Pre", "Non_Remission_Post",
    "Remission_Pre", "Remission_Post", "NA_None"
  )))


##-----------------------------------------------------------------------------
table(gd$res_tp)


##-----------------------------------------------------------------------------
gd@meta.data <- gd@meta.data%>%
  mutate(Remission_status=str_replace(Remission_status, "_", "-"))%>%
  unite("resp_tp", Remission_status, Treatment, remove=F)%>%
  unite("resp_tp_infl", Remission_status, Treatment, Inflammation, remove=F)


##-----------------------------------------------------------------------------
table(gd$resp_tp_infl,
      gd$Disease)


##-----------------------------------------------------------------------------
gd <- subset(gd, subset=Disease!="Healthy")


##-----------------------------------------------------------------------------
gd_pb <- AggregateExpression(
  gd,
  group.by=c("Remission_status", "Treatment", "Patient"),
  return.seurat=T
)


##-----------------------------------------------------------------------------
table(gd_pb$Patient)

##-----------------------------------------------------------------------------
# Remove Patients without paired data
gd_pb <- subset(gd_pb, subset=Patient %in% c("UC16", "UC20"), invert=T)


##-----------------------------------------------------------------------------
table(gd_pb$Patient)

## gd FLASHseq for Scoring
##-----------------------------------------------------------------------------
gd_flash <- readRDS("../../data_processed/FLASHseq/gd.Rds")


##-----------------------------------------------------------------------------
Idents(gd_flash) <- "celltype"

cell_marker <- FindAllMarkers(gd_flash, only.pos=T, min.pct=0.05, logfc.threshold=0.5)%>%
  filter(p_val_adj<0.05)


##-----------------------------------------------------------------------------
gene_list <- split(cell_marker$gene, cell_marker$cluster)

for (celltype in names(gene_list)) {
  gd_pb <- AddModuleScore(gd_pb, features=list(gene_list[[celltype]]), name=celltype)
}


# Saving Objects
##-----------------------------------------------------------------------------
# saveRDS(pb, "../../data_processed/Thomas_24/thomas24_entero_pb.Rds")
# saveRDS(gd_pb, "../../data_processed/Thomas_24/thomas24_gd_pb.Rds")


# Session Info
##-----------------------------------------------------------------------------
sessionInfo()





















