---
title: "FLASH-seq: Code for Figures"
output: github_document
---



##-----------------------------------------------------------------------------
library(Seurat)
library(SCpubr)
library(scCustomize)
library(tidyverse)
library(scRepertoire)
library(ReactomePA)
library(enrichplot)
library(circlize)
library(patchwork)
library(ggpubr)
library(rstatix)
library(SCENIC)
library(AUCell)
library(paletteer)
library(scales)


# Object
##-----------------------------------------------------------------------------
gd <- readRDS("../../data_processed/FLASHseq/gd.Rds")


# Marker Lists
##-----------------------------------------------------------------------------
Idents(gd) <- "seurat_clusters"

marker <- FindAllMarkers(gd, only.pos=T)%>%
  filter(p_val_adj<0.05)

##-----------------------------------------------------------------------------
write.csv(marker, "FLASHseq_Marker_AllClusters.csv")


# Palettes
##-----------------------------------------------------------------------------
pal_seur <- c("#4682B4", "#AF46B4", "gray40", "#B47846", "#4BB446", "red3")
pal_seur <- setNames(pal_seur, 0:5)
pal_clus <- setNames(pal_seur, paste0("C", 0:5))

pal_cell <- c("IEL"="#4682B4", "Memory-Like"="#B47846", "Cytotoxic"="#4BB446",
              "Cycling"="red3")

pal_dis <- c("HD"="#008B8B","UC"="#8B0000")
names_dis <- c("UC"="UC (A)")

pal_trdv <- setNames(c("#8C564BFF", "#9467BDFF", "#D62728FF", "#2CA02CFF", 
                       "#FF7F0EFF", "#1F77B4FF"), levels(gd$TRDV))

pal_trgv <- setNames(c("#BCBD22FF", "#7F7F7FFF", "#E377C2FF", "#8C564BFF", "pink",
                       "#9467BDFF","#D62728FF", "#2CA02CFF", "#FF7F0EFF", 
                       "#1F77B4FF"), levels(gd$TRGV))


# Figure 1
##-----------------------------------------------------------------------------
# 1b
p <- do_DimPlot(gd, label=T, group.by="seurat_clusters", legend.position="none",
           colors.use=pal_seur, pt.size=.8)

p[[1]][["layers"]][[3]][["geom"]][["default_aes"]][["alpha"]] <- .7
p

ggsave("../figures/FLASHseq_UMAP_by_Cluster.pdf", width=2.7, height=3)
ggsave("../figures/FLASHseq_UMAP_by_Cluster.png", width=2.7, height=3, dpi=600)
rm(p)


##-----------------------------------------------------------------------------
# 1c
do_DimPlot(gd, group.by="disease", legend.position="left", pt.size=.8)+
  scale_color_manual(values=pal_dis ,labels=names_dis)+
  theme(legend.position=c(.82, .95), legend.background=element_blank())

ggsave("../figures/FLASHseq_UMAP_by_disease.pdf", width=2.7, height=3)
ggsave("../figures/FLASHseq_UMAP_by_disease.png", width=2.7, height=3, dpi=600)


##-----------------------------------------------------------------------------
# 1c
do_DimPlot(gd, group.by="celltype", pt.size=.8,colors.use=pal_cell,
           legend.icon.size=3, legend.ncol=2, legend.byrow=F)+
  theme(legend.position=c(.45, 1.05), legend.background=element_blank(),
        legend.key.height=unit(1, "mm"), legend.key.spacing.y=unit(1, "mm"))

  
ggsave("../figures/FLASHseq_UMAP_by_CellType.pdf", width=2.7, height=3)
ggsave("../figures/FLASHseq_UMAP_by_CellType.png", width=2.7, height=3, dpi=600)



##-----------------------------------------------------------------------------
# 1d
Idents(gd) <- "cluster"

genes <- c("ITGAE","ITGA1","CD160","GZMA","ID3","KIT","IRF8","TCF7", "IL7R", 
           "PDCD1", "GZMK","CD5","CD27","CCR7","ZNF683","TNF","IFNG", "EOMES",
           "S1PR1","SELL","FCGR3A", "GNLY", "NKG7","GZMB", "PRF1","TBX21",
           "ZEB2","KLRG1", "ITGB2", "MKI67")

p <- DotPlot(gd, group.by="cluster", features=genes, scale.by="radius")+
  scale_color_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  guides(color=guide_colorbar(frame.colour="black", frame.linewidth=.5, order=1,
          ticks.colour="black", ticks.linewidth=.5, title="Avg.\nScal.\nExpr."),
         size=guide_legend(override.aes=list(fill="white", shape=21),title="% Expr.",
                           order=2)
        )+
  geom_point(aes(size=(pct.exp+6)), shape=21, fill=NA, colour="black", stroke=.25)+
  theme(axis.title=element_blank(),axis.ticks.length.x=unit(1.5, "mm"),
        panel.grid.major=element_line(color="gray95"), legend.justification="top")+
  coord_flip()+RotatedAxis()


p_1 <- gd@meta.data%>%mutate(cluster=factor(cluster, levels=levels(p$data$id)))

p_1 <- ggplot(p_1, aes(x=cluster, fill=disease))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_dis, labels=names_dis)+
  theme_void()+
  theme(text=element_text(size=18),legend.title=element_blank(),
        legend.position=c(-.4, .5), legend.key.size=unit(4, "mm"), 
        legend.spacing.y=unit(0.5, "mm"), legend.margin=margin(t=-8)
        )+
  guides(fill=guide_legend(ncol=1, byrow=T))+
  geom_hline(yintercept=.5, linetype="dashed", linewidth=.75)

p/p_1+plot_layout(heights=c(1, 0.05))

ggsave("../figures/FLASHseq_DotPlot_1.pdf", width=3.2, height=6.75)
ggsave("../figures/FLASHseq_DotPlot_1.png", width=3.2, height=6.75, dpi=600)

rm(p, p_1, genes)


##-----------------------------------------------------------------------------
# 1e+1f
p <- VlnPlot_scCustom(gd, colors_use=pal_clus, pt.size=0, sort=T, group.by="cluster",
        features=c("cytokine_cytotoxic_signature1", "TRM_signature1"))&
  geom_boxplot(width=.3, fill="white", outlier.shape=NA, coef=0, color="black")&
  NoLegend()&
  scale_y_continuous(expand=expansion(c(0.05, 0.05)))&
  labs(y="Module Score", x=NULL)&
  theme(plot.title=element_text(hjust=0))

for (i in seq_along(p)){
  p[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- F
}

p[[1]]+ggtitle(label="Cytotoxicity and Cytokine\nProduction Signature",
               subtitle=expression(italic("Szabo et al., 2019")))

ggsave("../figures/FLASHseq_VlnPlot_Cyto.pdf", width=4, height=3.5)
ggsave("../figures/FLASHseq_VlnPlot_Cyto.png", width=4, height=3.5, dpi=600)

p[[2]]+ggtitle(label="Tissue Residency\nSignature",
               subtitle=expression(italic("Kumar et al., 2017")))

ggsave("../figures/FLASHseq_VlnPlot_TRM.pdf", width=4, height=3.5)
ggsave("../figures/FLASHseq_VlnPlot_TRM.png", width=4, height=3.5, dpi=600)

rm(i, p)


##-----------------------------------------------------------------------------
# 1g+1h
do_DimPlot(gd, group.by="celltype", pt.size=2, font.size=24,
           reduction="scenic_umap", colors.use=pal_cell, legend.icon.size=5)+
  theme(legend.position=c(0.2, 0.97), legend.background=element_blank(),
        legend.key.spacing.y=unit(2, "mm"))

ggsave("../figures/FLASHseq_SCENIC_UMAP_by_CellType.pdf", width=6, height=5)
ggsave("../figures/FLASHseq_SCENIC_UMAP_by_CellType.png", dpi=600, width=6, height=5)

do_DimPlot(gd, group.by="disease",pt.size=2, font.size=24, legend.icon.size=5,
           reduction="scenic_umap")+
  scale_color_manual(values=pal_dis, labels=names_dis)+
  theme(legend.position=c(0.15, 0.95), legend.key.spacing.y=unit(2, "mm"))

ggsave("../figures/FLASHseq_SCENIC_UMAP_by_Disease.pdf", width=6, height=5)
ggsave("../figures/FLASHseq_SCENIC_UMAP_by_Disease.png", dpi=600, width=6, height=5)

##-----------------------------------------------------------------------------
# Supps

do_DimPlot(gd, group.by="cluster", pt.size=2, font.size=24, legend.ncol=2,
           reduction="scenic_umap", colors.use=pal_clus, legend.icon.size=5)+
  theme(legend.position=c(0.2, 0.97), legend.background=element_blank(),
        legend.key.spacing.y=unit(2, "mm"))

ggsave("../figures/FLASHseq_SCENIC_UMAP_by_Cluster.pdf", width=6, height=5)
ggsave("../figures/FLASHseq_SCENIC_UMAP_by_Cluster.png", dpi=600, width=6, height=5)



##-----------------------------------------------------------------------------
# 1i
# Import AUC scores
rss <- importAUCfromText(
  "../../data_processed/FLASHseq/pySCENIC/auc.csv")
rss@NAMES <- gsub(rss@NAMES, pattern="\\(\\+\\)", replacement="")

rss <- calcRSS(AUC=getAUC(rss),cellAnnotation=gd@meta.data$celltype)%>%
  as.data.frame()%>%
  rownames_to_column("features.plot")%>%
  pivot_longer(cols=levels(gd$celltype), names_to="id", values_to="RSS")

df <- DotPlot(gd, assay="scenic", group.by="celltype",
              features=unique(rss$features.plot))$data%>%
  inner_join(rss, by=c("features.plot", "id"))



##-----------------------------------------------------------------------------
# DotPlot
regulons <- c("GATA3", "SOX10", "AHR", "IRF4", "IRF5", "IRF8", "IRF9", "STAT3",
              "STAT6", "ATF3", "TCF7", "LEF1", "KLF2", "EOMES", "TBX21", "ETS1",
              "FOXM1", "MYB", "E2F1")

df%>%
  filter(features.plot %in% regulons)%>%
  mutate(features.plot=factor(features.plot, levels=rev(regulons)),
         id=factor(id, levels=levels(gd$celltype))
         )%>%
  
ggplot(aes(x=id, y=features.plot, size=RSS, fill=avg.exp.scaled))+
  geom_point(color="black", shape=21)+
  theme_minimal()+
  scale_fill_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  guides(fill=guide_colorbar(
    frame.colour="black", frame.linewidth=.5,
    ticks.colour="black", ticks.linewidth=.5,
    order=1, title="Avg.\nScaled\nAUC"
  ))+
  theme(axis.text=element_text(color="black"),axis.title=element_blank())+
  RotatedAxis()

ggsave("../figures/FLASHseq_SCENIC_RSS_by_CellType.pdf", width=2.5, height=5)
ggsave("../figures/FLASHseq_SCENIC_RSS_by_CellType.png", dpi=600 ,width=2.5, height=5)

rm(rss, regulons, df)


# Figure 3

##-----------------------------------------------------------------------------
# 3a
do_DimPlot(gd, group.by="TRDV", pt.size=.8, font.size=12, legend.icon.size=3,
           legend.ncol=3,
           idents.keep=unique(gd$TRDV[!is.na(gd$TRDV)]),
           colors.use=pal_trdv, na.value="gray80")+

  theme(legend.key.spacing.y=unit(1.5, "mm"), legend.key.height=unit(1.5, "mm"),
        legend.margin=margin(t=-15))

ggsave("../figures/FLASHseq_UMAP_by_TRDV.pdf", width=3.5, height=3)
ggsave("../figures/FLASHseq_UMAP_by_TRDV.png", width=3.5, height=3, dpi=600)



##-----------------------------------------------------------------------------
# 3d
do_DimPlot(gd, group.by="TRGV", pt.size=.8, font.size=12, legend.icon.size=3,
           legend.ncol=3,
           idents.keep=unique(gd$TRGV[!is.na(gd$TRGV)]),
           colors.use=pal_trgv, na.value="gray80")+

  theme(legend.key.spacing.y=unit(1.5, "mm"), legend.key.height=unit(1.5, "mm"),
        legend.margin=margin(t=-15))


ggsave("../figures/FLASHseq_UMAP_by_TRGV.pdf", width=3.5, height=3)
ggsave("../figures/FLASHseq_UMAP_by_TRGV.png", width=3.5, height=3, dpi=600)



##-----------------------------------------------------------------------------
arranged <- gd@meta.data%>%count(TRGV_TRDV)%>%arrange(desc(n))
arranged$TRGV_TRDV <- factor(arranged$TRGV_TRDV, levels=unique(arranged$TRGV_TRDV)[-1])

pal_gvdv <- setNames(c("#1F77B4FF", "#8C564BFF", "darkgreen", "darkred", "darkorange3", "gray30"),
                     c(levels(arranged$TRGV_TRDV)[1:5], "other"))
names_gvdv <- setNames(str_replace(arranged$TRGV_TRDV[2:6], "_", "-"),arranged$TRGV_TRDV[2:6])



##-----------------------------------------------------------------------------
# 3g
p <- do_DimPlot(gd, group.by="TRGV_TRDV", pt.size=.8, font.size=12, legend.icon.size=3,
           legend.ncol=2,
           idents.keep=arranged$TRGV_TRDV[2:6],
           na.value="gray80")+
  scale_color_manual(values=pal_gvdv, labels=names_gvdv)+
  theme(legend.key.spacing.y=unit(1.5, "mm"), legend.key.height=unit(1.5, "mm"),
        legend.margin=margin(t=-15))

p[[1]]$data$TRGV_TRDV <- factor(p[[1]]$data$TRGV_TRDV, levels=arranged$TRGV_TRDV[2:6])
p

ggsave("../figures/FLASHseq_UMAP_by_TRGV-TRDV.pdf", width=3.5, height=3)
ggsave("../figures/FLASHseq_UMAP_by_TRGV-TRDV.png", dpi=600, width=3.5, height=3)

rm(p)


##-----------------------------------------------------------------------------
# 3b
ggplot(gd@meta.data, aes(x=disease, fill=fct_rev(TRDV)))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_trdv, na.value="gray80")+
  scale_x_discrete(labels=names_dis)+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black"), legend.margin=margin(l=-7),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"))+
  RotatedAxis()
ggsave("../figures/FLASHseq_BarPlot_TRDV_by_Disease.pdf", width=2.35, height=3)
ggsave("../figures/FLASHseq_BarPlot_TRDV_by_Disease.png", width=2.35, height=3, dpi=600)



##-----------------------------------------------------------------------------
gd@meta.data%>%filter(celltype!="Cycling")%>%
ggplot(aes(x=celltype, fill=fct_rev(TRDV)))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_trdv, na.value="gray80")+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black"),
        legend.margin=margin(l=-7),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"))+
  RotatedAxis()
ggsave("../figures/FLASHseq_BarPlot_TRDV_by_CellType.pdf", width=2.75, height=3)
ggsave("../figures/FLASHseq_BarPlot_TRDV_by_CellType.png", width=2.75, height=3, dpi=600)


##-----------------------------------------------------------------------------
# 3e
ggplot(gd@meta.data, aes(x=disease, fill=fct_rev(TRGV)))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_trgv, na.value="gray80")+
  scale_x_discrete(labels=names_dis)+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), 
        axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"),
        legend.margin=margin(l=-7))+
  RotatedAxis()
ggsave("../figures/FLASHseq_BarPlot_TRGV_by_Disease.pdf", width=2, height=3)
ggsave("../figures/FLASHseq_BarPlot_TRGV_by_Disease.png", width=2, height=3, dpi=600)



##-----------------------------------------------------------------------------
gd@meta.data%>%filter(celltype!="Cycling")%>%
ggplot(aes(x=celltype, fill=fct_rev(TRGV)))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_trgv, na.value="gray80")+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), 
        axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"),
        legend.margin=margin(l=-7))+
  RotatedAxis()
ggsave("../figures/FLASHseq_BarPlot_TRGV_by_CellType.pdf", width=2.25, height=3.25)
ggsave("../figures/FLASHseq_BarPlot_TRGV_by_CellType.png", width=2.25, height=3.25, dpi=600)


##-----------------------------------------------------------------------------
# 3h
gd@meta.data%>%
  mutate(TRGV_TRDV=case_when(
    TRGV_TRDV %in% levels(arranged$TRGV_TRDV)[1:5] ~ TRGV_TRDV,
    TRGV_TRDV %in% levels(arranged$TRGV_TRDV[-1:-5]) ~ "other"
  ),
          TRGV_TRDV=factor(TRGV_TRDV, levels=c(levels(arranged$TRGV_TRDV[1:5]), "other"))        
  )%>%
  filter(celltype!="Cycling")%>%
ggplot(aes(x=celltype, fill=TRGV_TRDV))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_gvdv, labels=names_gvdv, na.value="gray80")+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=12), legend.margin=margin(l=-7),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"))+
  RotatedAxis()

ggsave("../figures/FLASHseq_BarPlot_TRDV-TRGV_by_CellType.pdf", width=3, height=4)
ggsave("../figures/FLASHseq_BarPlot_TRDV-TRGV_by_CellType.png", width=3, height=4, dpi=600)


##-----------------------------------------------------------------------------
# # Test Function
# test <-   function(y.pos=NULL, x.pos=NULL){
#   stat_compare_means(bracket.size=0, method="t.test", label.y=y.pos, label.x=x.pos,
#   aes(label=ifelse(..p.. < 0.0001, "p<0.0001", paste0("p=", ..p.format..)))
#                     )
# }


##-----------------------------------------------------------------------------
# 3c
p <- VlnPlot_scCustom(gd, group.by="disease", features=c("TRDV2", "TRDV1"),
                 num_columns=1,colors_use=pal_dis, add.noise=F, pt.size=0)&
  geom_boxplot(fill="white", color="black", width=.25, coef=0, outlier.shape=NA)&
  scale_y_continuous(expand=expansion(c(0.2, 0.05)))&
  scale_x_discrete(labels=names_dis)&
  labs(x=NULL)&
  theme(axis.text.x=element_text(angle=0, hjust=.5))&
  stat_compare_means(bracket.size=0, method="t.test", label.x=1.27,
  aes(label=ifelse(..p.. < 0.0001, "p<0.0001", paste0("p=", ..p.format..)))
                    )

for (i in seq_along(p)) {
  p[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- F
}
p
ggsave("../figures/FLASHseq_VlnPlots_TRDVs_by_Disease.pdf", width=2.5, height=5)
ggsave("../figures/FLASHseq_VlnPlots_TRDVs_by_Disease.png", width=2.5, height=5, dpi=600)


##-----------------------------------------------------------------------------
# 3f
p <- VlnPlot_scCustom(gd, group.by="disease", features=c("TRGV4", "TRGV9"),
                 num_columns=1,colors_use=pal_dis, add.noise=F, pt.size=0)&
  geom_boxplot(fill="white", color="black", width=.25, coef=0, outlier.shape=NA)&
  scale_y_continuous(expand=expansion(c(0.2, 0.05)))&
  scale_x_discrete(labels=names_dis)&
  labs(x=NULL)&
  theme(axis.text.x=element_text(angle=0, hjust=.5))&
  stat_compare_means(bracket.size=0, method="t.test", label.x=1.27,
  aes(label=ifelse(..p.. < 0.0001, "p<0.0001", paste0("p=", ..p.format..)))
                    )

for (i in seq_along(p)) {
  p[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- F
}
p
ggsave("../figures/FLASHseq_VlnPlots_TRGVs_by_Disease.pdf", width=2.5, height=5)
ggsave("../figures/FLASHseq_VlnPlots_TRGVs_by_Disease.png", width=2.5, height=5, dpi=600)


##-----------------------------------------------------------------------------
# 3i
gd@meta.data %>%
  group_by(ID, disease)%>%
  summarise(
    total_valid=sum(!is.na(TRGV_TRDV)),
    TRGV9_TRDV2=sum(TRGV_TRDV == "TRGV9_TRDV2", na.rm=TRUE),
    TRGV4_TRDV1=sum(TRGV_TRDV == "TRGV4_TRDV1", na.rm=TRUE),
    .groups="drop"
  )%>%
  pivot_longer(cols=c(TRGV9_TRDV2, TRGV4_TRDV1), names_to="TRGV_TRDV", values_to="n")%>%
  mutate(percent=(n / total_valid) * 100)%>%

ggplot(aes(x=disease, y=percent, fill=disease))+
  geom_boxplot(color="black", outlier.shape=NA)+
  geom_point(position=position_jitter(seed=1337, height=0, width=0.3), shape=21, fill="white")+
  scale_fill_manual(values=pal_dis)+
  scale_x_discrete(labels=names_dis)+
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.1)))+
  theme_minimal()+
  labs(x=NULL, y="% of all subsets")+
  theme(axis.line.y=element_line(), axis.ticks.y=element_line(),
        panel.grid=element_blank(), axis.text=element_text(color="black"),
        legend.position="none", strip.text=element_text(face="bold", size=12))+
  facet_wrap(~TRGV_TRDV, scales="free_y")+
  stat_compare_means(bracket.size=0, method="t.test", label.x=1.6,
                     aes(label=ifelse(..p.. < 0.0001, "p<0.0001", paste0("p=", ..p.format..)))
                     )

ggsave("../figures/FLASHseq_BarPlot_TRGVs-TRDVs_by_Disease.pdf", width=3.5, height=3)
ggsave("../figures/FLASHseq_BarPlot_TRGVs-TRDVs_by_Disease.png", width=3.5, height=3, dpi=600)


##-----------------------------------------------------------------------------
# 3j
chord.df <- getCirclize(gd, "CTstrict_filled", group.by="cluster",proportion=F)%>%
              arrange(from)

pdf("../figures/FLASHseq_CloneOverlap_byCluster.pdf")

par(cex=2.5, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack=c("name", "grid"),
             self.link=1, grid.col=pal_clus)
title("Clonal Overlap")

dev.off()
 #
png("../figures/FLASHseq_CloneOverlap_byCluster.png", res=600, width=4200, height=4200)

par(cex=2.5, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack=c("name", "grid"),
             self.link=1, grid.col=pal_clus)
title("Clonal Overlap")

dev.off()

rm(chord.df)


##-----------------------------------------------------------------------------
# 3j
chord.df <- getCirclize(gd, "CTstrict_filled", group.by="celltype")

pdf("../figures/FLASHseq_CloneOverlap_byCellType.pdf")

par(cex=2.2, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack=c("name", "grid"),
             self.link=1, grid.col=pal_cell)
title("Clonal Overlap")

dev.off()
 #
png("../figures/FLASHseq_CloneOverlap_by_CellType.png", res=600, width=4200, height=4200)

par(cex=2.2, mar=c(0 , 0 , 2, 0))
chordDiagram(chord.df, annotationTrack=c("name", "grid"),
             self.link=1, grid.col=pal_cell)
title("Clonal Overlap")

dev.off()

rm(chord.df)



##-----------------------------------------------------------------------------
# 3k
shared_clones <- intersect(
  gd$CTstrict_filled[gd$seurat_clusters=="3"],
  gd$CTstrict_filled[gd$seurat_clusters=="4"]
)[-1]


##-----------------------------------------------------------------------------
process_clone <- function(clone) {
  # Split the string into gamma and delta chains
  parts <- unlist(strsplit(clone, "_"))
  
  # Handle case where one of the chains is missing (NA)
  if (length(parts) == 1) {
    # Only one chain exists (either gamma or delta is NA)
    gamma_chain <- "NA"
    delta_chain <- parts[1]
  } else {
    # Two chains exist, gamma and delta
    gamma_chain <- parts[2]
    delta_chain <- parts[1]
  }
  
  # Split the chains into their respective segments
  gamma_segments <- strsplit(gamma_chain, ";")[[1]]
  delta_segments <- strsplit(delta_chain, ";")[[1]]
  
  # Extract CDR3 region and calculate its length
  gamma_cdr3 <- if (length(gamma_segments) > 1) gamma_segments[2] else ""
  delta_cdr3 <- if (length(delta_segments) > 1) delta_segments[2] else ""
  
  gamma_cdr3_length <- nchar(gamma_cdr3)
  delta_cdr3_length <- nchar(delta_cdr3)
  
  # Check if any chain is missing, handle accordingly
  if (gamma_chain == "NA") {
    gamma_label <- ""
    delta_label <- paste(delta_segments[1], "_CDR3(", delta_cdr3_length, " nt)", sep="")
    result <- paste(delta_label, "(TRG NA)", sep=" ")
  } else if (delta_chain == "NA") {
    delta_label <- ""
    gamma_label <- paste(gamma_segments[1], "_CDR3(", gamma_cdr3_length, " nt)", sep="")
    result <- paste(gamma_label, "(TRD NA)", sep=" ")
  } else {
    # Both chains exist
    gamma_label <- paste(gamma_segments[1], "_CDR3(", gamma_cdr3_length, " nt)", sep="")
    delta_label <- paste(delta_segments[1], "_CDR3(", delta_cdr3_length, " nt)", sep="")
    result <- paste(gamma_label, delta_label, sep=" ")
  }
  

  return(result)
}


clones_named <- sapply(shared_clones, process_clone)

names(clones_named) <- shared_clones



##-----------------------------------------------------------------------------
colors <- setNames(c("red", "blue", "green", "purple", "yellow"), shared_clones)

gd <- SetIdent(gd, value="CTstrict_filled")

do_DimPlot(gd, idents.keep=shared_clones, pt.size=3, na.value="gray80",
           legend.ncol=1, border.size=1.5)+
  scale_color_manual(values=colors, labels=clones_named)+
  theme(legend.key.spacing.y=unit(0, "mm"))

ggsave("../figures/FLASHseq_UMAP_by_C3C4_Clones.pdf", width=7, height=7)
ggsave("../figures/FLASHseq_UMAP_by_C3C4_Clones.png", dpi=600, width=7, height=7)

rm(shared_clones, clones_named, colors, process_clone, clones, i)


# Figure 4
##-----------------------------------------------------------------------------
# 4a
Idents(gd) <- "disease"

UC_markers <- FindMarkers(gd, ident.1="UC")%>%
  rownames_to_column("gene")

geneLogSums <- log2(Matrix::rowMeans(GetAssayData(gd, "RNA", "counts")))
geneLogSums <- as.data.frame(geneLogSums)%>%
  rownames_to_column("gene")

UC_markers <- inner_join(UC_markers, geneLogSums, by="gene")%>%
  mutate(signif=case_when(
           p_val_adj<0.05 & avg_log2FC>0 ~ "up",
           p_val_adj<0.05 & avg_log2FC<0 ~ "down",
           .default="ns"
         ),
         signif=factor(signif, levels=c("ns", "down", "up"))
         )

genes_up <- c("TCF7", "IL7R","GZMK","CD5","CD27","CCR7","ZNF683","TNF",
              "IFNG","EOMES","S1PR1","SELL","NKG7","GZMB","TBX21","KLRG1", "ITGB2", 
              "PDCD1", "CTLA4", "LAG3")

genes_down <- c("ITGAE","ITGA1","CD160","ID3")

ggplot()+
  geom_point(UC_markers%>%filter(signif=="ns"),
             mapping=aes(x=geneLogSums, y=avg_log2FC, color=signif))+
  geom_point(UC_markers%>%filter(signif!="ns"),
             mapping=aes(x=geneLogSums, y=avg_log2FC, color=signif))+
  theme_classic()+
  geom_hline(yintercept=0, linetype="dashed", color="red", size=1)+
  labs(x=expression("log"[2]*"(Mean Expression)"),
       y=expression("log"[2]*"(Fold Change)"),
       color=element_blank())+
  scale_color_manual(values=c("gray", "#008B8B", "#8B0000"))+
  guides(color=guide_legend(override.aes=list(size=2), reverse=T))+
  theme(axis.text=element_text(color="black"), text=element_text(size=16),
        legend.position=c(0.15, 0.85), legend.text=element_text(size=16))+
  
  ggrepel::geom_label_repel(data=UC_markers%>%filter(gene %in% genes_up), 
                           aes(label=gene, x=geneLogSums, y=avg_log2FC),
                           min.segment.length=0, nudge_y=0, nudge_x=0,
                           size=3, max.overlaps=15)+
  ggrepel::geom_label_repel(data=UC_markers%>%filter(gene %in% genes_down), 
                           aes(label=gene, x=geneLogSums, y=avg_log2FC),
                           min.segment.length=0, nudge_y=-1, nudge_x=0,
                           size=3)

ggsave("../figures/FLASHseq_MAPlot.pdf", width=5, height=5)
ggsave("../figures/FLASHseq_MAPlot.png", width=5, height=5, dpi=600)

rm(genes_down, genes_up, UC_markers, geneLogSums)


##-----------------------------------------------------------------------------
# 4b
gsea <- readRDS("../../data_processed/FLASHseq/GSEA_res.Rds")

gsea@result%>%filter(ID %in% c("R-HSA-389948","R-HSA-202427","R-HSA-388841",
    "R-HSA-202430","R-HSA-6785807","R-HSA-449147","R-HSA-2132295",
    "R-HSA-877300", "R-HSA-500792"))%>%

ggplot(aes(x=NES, y=reorder(Description, NES), fill=-log10(p.adjust)))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_fill_gradientn(limits=c(1.3, 8.7), colors=c("#E8FFFF", "#00007A"))+
  scale_y_discrete(labels=label_wrap(25))+
  guides(fill=guide_colorbar(frame.colour="black", frame.linewidth=.5,
                             ticks.colour="black", ticks.linewidth=.5))+
  labs(x="Normalized Enrichment Score (NES)",y=NULL,
       fill=expression("-log"[10]*"(P)"), title="Gene Set Enrichment")+
  
  theme(panel.grid=element_blank(), panel.border=element_rect(color="black", size=1),
    axis.text=element_text(color="black"),plot.title=element_text(face="bold", hjust=.5)
        )

ggsave("../figures/FLASHseq_GSEA_NESplot.pdf", width=4.5, height=3.6)
ggsave("../figures/FLASHseq_GSEA_NESplot.png", width=4.5, height=3.6, dpi=600)



##-----------------------------------------------------------------------------
# 4c
gseaplot(gsea, geneSetID="R-HSA-389948", by="runningScore",
         title=gsea@result%>%filter(ID=="R-HSA-389948")%>%pull(Description))+
  annotate(geom="text", x=32000, y=.7, label=paste0("p=",
    gsea@result%>%filter(ID=="R-HSA-389948")%>%pull(p.adjust)%>%signif(2)),
    size=5)+
  theme(plot.title=element_text(face="bold"))

ggsave("../figures/FLASHseq_GSEAplot_PD1.pdf", width=7, height=4.5)
ggsave("../figures/FLASHseq_GSEAplot_PD1.png", width=7, height=4.5, dpi=600)

#
gseaplot(gsea, geneSetID="R-HSA-201722", by="runningScore",
         title=gsea@result%>%filter(ID=="R-HSA-201722")%>%pull(Description))+
  annotate(geom="text", x=32000, y=.45, label=paste0("p=",
    gsea@result%>%filter(ID=="R-HSA-201722")%>%pull(p.adjust)%>%signif(2)),
    size=5)+
  theme(plot.title=element_text(face="bold"))

ggsave("../figures/FLASHseq_GSEAplot_bCat_TCF.pdf", width=7, height=4.5)
ggsave("../figures/FLASHseq_GSEAplot_bCat_TCF.png", width=7, height=4.5, dpi=600)

# rm(gsea)


##-----------------------------------------------------------------------------
# 4d, 4e
heat <- DotPlot(gd, group.by="celltype",
          features=c("CD28", "ICOS", "TNFRSF4", "TNFRSF14", "TNFRSF9",
                     "PDCD1", "LAG3", "CD244", "CTLA4", "HAVCR2"))$data

p1 <- heat%>%filter(features.plot %in% c("CD28", "ICOS", "TNFRSF4", "TNFRSF14", "TNFRSF9"))%>%
  ggplot(aes(x=id, y=features.plot,fill=avg.exp.scaled))+
  geom_tile(color="white", size=.25)+
  theme_minimal()+
  scale_fill_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  guides(fill=guide_colorbar(frame.colour="black", frame.linewidth=.5,
                             ticks.colour="black", ticks.linewidth=.25
                             ))+
  labs(fill="Avg.\nExpr.\nScal.", title="Costimulatory\nMolecules")+
  theme(panel.grid.major=element_blank(), 
        panel.border=element_rect(color="black", fill=NA, size=1),
        axis.title=element_blank(), axis.text=element_text(color="black"),
        legend.key.width=unit(.5, "cm"), legend.key.height=unit(.45, "cm"),
        legend.justification="top")+
  RotatedAxis()


p2 <- heat%>%filter(features.plot %in% c("PDCD1", "LAG3", "CD244", "CTLA4", "HAVCR2"))%>%
  ggplot(aes(x=id, features.plot,fill=avg.exp.scaled))+
  geom_tile(color="white", size=.25)+
  theme_minimal()+
  scale_fill_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  labs(fill="Avg.\nExpr.\nScal.", title="Inhibitory\nReceptors")+
  theme(panel.grid.major=element_blank(), 
        panel.border=element_rect(color="black", fill=NA, size=1),
        axis.title=element_blank(), axis.text=element_text(color="black"))+
  RotatedAxis()

ggarrange(p1, p2, common.legend=T, legend="right", widths=c(1, .95))

ggsave("../figures/FLASHseq_HeatMap_1.pdf", width=4.3, height=2.6)
ggsave("../figures/FLASHseq_HeatMap_1.png", width=4.3, height=2.75, dpi=600)

rm(p1, p2, heat)


##-----------------------------------------------------------------------------
# 4f-4j
scores <- grep(colnames(gd@meta.data), pattern="stem1$", value=T)

p <- subset(gd, subset=celltype!="Cycling")%>%
  VlnPlot_scCustom(features=scores, group.by="celltype",
                 pt.size=0, colors_use=pal_cell, sort=T,
                 add.noise=F, num_columns=2)&
  scale_y_continuous(expand=expansion(mult=c(.05, .1)))&
  scale_x_discrete(labels=c("Memory-Like"="Memory-\nLike"))&
  theme(axis.title=element_blank(),
        axis.text.x=element_text(angle=0, hjust=0.5),
        plot.margin=margin(t=10),
        plot.subtitle=element_text(face="italic"),
        plot.title=element_text(hjust=0))&
  geom_boxplot(width=.25, fill="white", outlier.shape=NA, coef=0, color="black", size=.35)

p[[1]] <- p[[1]] + labs(title="Tcf7GFP+ Tumor Infil-\ntrating T Cell Signature",
                        subtitle="Siddiqui et al., 2019")

p[[2]] <- p[[2]] + labs(title="CXCR5+ CD8+ T Cell\nSignature",
                        subtitle="Im et al., 2016")

p[[3]] <- p[[3]] + labs(title="Tim3-Blimp-CD8+\nT Cell Signature",
                        subtitle="Wu et al., 2016")

p[[4]] <- p[[4]] + labs(title="Tcf1+ Memory-Like\nCD8+ T Cell Signature",
                        subtitle="Utzschneider et al., 2016")

p[[5]] <- p[[5]] + labs(title="Stem-Like alpha beta\nT Cell Signature in UC",
                        subtitle="Li et al., 2024")

for (i in seq_along(p)) {
  p[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE
}

p[[5]]+p[[4]]

ggsave("../figures/FLASHseq_VlnPlots_StemLike_Signatures_1.pdf", width=6, height=3.333)
ggsave("../figures/FLASHseq_VlnPlots_StemLike_Signatures_1.png", width=6, height=3.333)

p[[1]]+p[[2]]+p[[3]]

ggsave("../figures/FLASHseq_VlnPlots_StemLike_Signatures_2.pdf", width=9, height=3.333)
ggsave("../figures/FLASHseq_VlnPlots_StemLike_Signatures_2.png", width=9, height=3.333)

rm(p, i, scores)


##-----------------------------------------------------------------------------
# 4k, 4l
heat <- DotPlot(gd, group.by="celltype",
  features=c("TNF", "TNFSF10", "IFNG", "GZMB", "PRF1"))$data

p1 <- ggplot(heat, aes(x=id, y=features.plot,fill=avg.exp.scaled))+
  geom_tile(color="white", size=.25)+
  theme_minimal()+
  scale_fill_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  NoLegend()+
  labs(title="Cytokines and\nEffector\nMolecules")+
  theme(panel.grid.major=element_blank(), panel.border=element_rect(color="black", fill=NA, size=1),
        axis.title=element_blank(), axis.text=element_text(color="black"),
        legend.key.width=unit(.5, "cm"), legend.key.height=unit(.45, "cm"),
        legend.justification="top")+
  RotatedAxis()

heat <- DotPlot(gd, group.by="celltype",
    features=c("TOX", "ID3", "ID2", "FOXO1", "TCF7", "BCL6", "EOMES", "TBX21"))$data

p2 <- ggplot(heat, aes(x=id, y=features.plot,fill=avg.exp.scaled))+
  geom_tile(color="white", size=.25)+
  theme_minimal()+
  scale_fill_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  guides(fill=guide_colorbar(frame.colour="black", frame.linewidth=.5,
                             ticks.colour="black", ticks.linewidth=.25
                             ))+
  labs(fill="Avg.\nExpr.\nScal.", title="Transcription\nFactors")+
  theme(panel.grid.major=element_blank(), panel.border=element_rect(color="black", fill=NA, size=1),
        axis.title=element_blank(), axis.text=element_text(color="black"))+
  RotatedAxis()

ggarrange(ggarrange(p1, NULL, heights=c(1, 0.14), nrow=2), 
          p2, common.legend=T, legend="right", widths=c(1, .95))

ggsave("../figures/FLASHseq_HeatMap_2.pdf", width=4, height=3)
ggsave("../figures/FLASHseq_HeatMap_2.png", width=4, height=3, dpi=600)

rm(p1, p2, heat)


# Supplements & New Figures
##-----------------------------------------------------------------------------
do_DimPlot(gd, split.by="ID")&NoLegend()
ggsave("../figures/FLASHseq_SplitBy_Donor.pdf", width=10, height=7)
ggsave("../figures/FLASHseq_SplitBy_Donor.png", width=10, height=7, dpi=600)


##-----------------------------------------------------------------------------
p <- FeaturePlot_scCustom(gd, features=c("TRDV1", "TRDV2", "TRGV4", "TRGV9"), pt.size=.8,
                     na_color="#2A0134", split.by="disease",
  colors_use=c("#6900A8FF", "#B02991FF", "#E16462FF", "#FCA636FF", "#F0F921FF"))&
  NoAxes()&
  guides(color=guide_colorbar(ticks.colour="black", frame.colour="black",
                              legend.key.height=unit(2, "mm"))
         )

for (i in seq_along(p)) {
  p[[i]][[1]] <- p[[i]][[1]]+NoLegend()
}
p

ggsave("../figures/FLASHseq_FeaturePlot_TRgenes.pdf", width=7, height=9)
ggsave("../figures/FLASHseq_FeaturePlot_TRgenes.png", width=7, height=9, dpi=600)



##-----------------------------------------------------------------------------
FeaturePlot_scCustom(gd, features=c("TCF7", "PDCD1"), pt.size=.8,
                     na_color="#2A0134",
  colors_use=c("#6900A8FF", "#B02991FF", "#E16462FF", "#FCA636FF", "#F0F921FF"))&
  NoLegend()&NoAxes()

ggsave("../figures/FLASHseq_FeaturePlot_StemLikeGenes.pdf", width=4.5, height=2.75)
ggsave("../figures/FLASHseq_FeaturePlot_StemLikeGenes.png", width=4.5, height=2.75, dpi=600)


##-----------------------------------------------------------------------------
# 
chord.df <- getCirclize(gd, "CTstrict_filled", group.by="ID",proportion=F)%>%
              arrange(from)

pdf("../figures/FLASHseq_CloneOverlap_byDonor.pdf")

par(cex=2, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack=c("grid", "name"),
             self.link=1, grid.col=SCpubr:::generate_color_scale(sort(unique(gd$ID))))
title("Clonal Overlap (Donors)")

dev.off()
 #
png("../figures/FLASHseq_CloneOverlap_byDonor.png", res=600, width=4200, height=4200)

par(cex=2, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack=c("grid", "name"),
             self.link=1, grid.col=SCpubr:::generate_color_scale(sort(unique(gd$ID))))
title("Clonal Overlap (Donors)")

dev.off()

rm(chord.df)



##-----------------------------------------------------------------------------
gd@meta.data%>%
  ggplot(aes(x=cluster, fill=ID))+
  geom_bar(position="fill", color="white",size=.2)+
  scale_fill_manual(values=c(rep("#008B8B", 4), rep("#8B0000", 6)))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(), legend.position="none",
        axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text=element_text(color="black"), axis.ticks.length=unit(2, "mm"),
        axis.ticks=element_line(color="black"), axis.text.x=element_text(vjust=3),
        )

ggsave("../figures/FLASHseq_BarPlot_ClusterByDonor.pdf", width=2, height=3)
ggsave("../figures/FLASHseq_BarPlot_ClusterByDonor.png", width=2, height=3, dpi=600)


##-----------------------------------------------------------------------------
set.seed(1337)
subset(gd, subset=celltype!="Cycling")%>%
clonalDiversity(cloneCall="CTstrict_filled", group.by="celltype", exportTable=T,
                n.boots=1000)%>%
  
  mutate(celltype=factor(celltype, levels=levels(gd$celltype)))%>%
  ggplot(aes(x=celltype, y=shannon, color=celltype))+
  geom_point(size=3, shape=21, color="black", aes(fill=celltype))+
  scale_fill_manual(values=pal_cell)+
  scale_y_continuous(expand=expansion(mult=c(.2, .2)))+
  theme_classic()+
  RotatedAxis()+
  NoLegend()+
  labs(y=NULL, x=NULL, title="Shannon\nEntropy")+
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"))

ggsave("../figures/FLASHseq_ShannonEntropy_byCellType.pdf", width=1.25, height=2)
ggsave("../figures/FLASHseq_ShannonEntropy_byCellType.png", width=1.25, height=2, dpi=600)



##-----------------------------------------------------------------------------
set.seed(1337)
clonalDiversity(gd, cloneCall="CTstrict_filled", group.by="ID", exportTable=T,
                n.boots=1000)%>%
  mutate(disease=case_when(
    str_detect(ID, "HD") ~ "HD",
    .default="UC"
  ))%>%

  ggplot(aes(x=disease, y=shannon,fill=disease))+
  geom_boxplot(outlier.shape=NA, color="black")+
  geom_point(size=1.5, position=position_jitter(width=.1, height=0), shape=21, fill="white")+
  scale_fill_manual(values=pal_dis)+
  scale_x_discrete(labels=names_dis)+
  scale_y_continuous(expand=expansion(mult=c(.2, .2)))+
  theme_classic()+
  NoLegend()+
  RotatedAxis()+
  labs(y=NULL, x=NULL, title="Shannon\nEntropy")+
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"))+
  stat_compare_means(method="t.test", label.y=2.4, label.x=1.3,
                     aes(label=ifelse(..p.. < 0.0001, "p<0.0001", paste0("p=", ..p.format..)))
                     )

ggsave("../figures/FLASHseq_ShannonEntropy_byDisease.pdf", width=1.25, height=2.25)
ggsave("../figures/FLASHseq_ShannonEntropy_byDisease.png", width=1.25, height=2.25, dpi=600)




##-----------------------------------------------------------------------------
sessionInfo()









