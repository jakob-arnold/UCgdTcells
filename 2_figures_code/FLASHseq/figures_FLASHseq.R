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


##-----------------------------------------------------------------------------
gd <- readRDS("../../data_processed/FLASHseq/gd.Rds")
Idents(gd) <- "cluster"


##-----------------------------------------------------------------------------
pal_seur <- c("0"="#4682B4", "1"="#AF46B4", "2"="#B47846", "3"="#4BB446")
pal_clu <- c("C0"="#4682B4", "C1"="#AF46B4", "C2"="#B47846", "C3"="#4BB446")
pal_dis <- c("control"="#008B8B","UC"="#8B0000")
names_dis <- c("control"="HD", "UC"="UC (A)")

pal_trdv <- setNames(c("#8C564BFF", "#9467BDFF", "#D62728FF", "#2CA02CFF", 
                       "#FF7F0EFF", "#1F77B4FF"), levels(gd$TRDV))

pal_trgv <- setNames(c("#BCBD22FF", "#7F7F7FFF", "#E377C2FF", "#8C564BFF", "#9467BDFF",
    "#D62728FF", "#2CA02CFF", "#FF7F0EFF", "#1F77B4FF"), levels(gd$TRGV))


# Figure 1

##-----------------------------------------------------------------------------
# 1b
do_DimPlot(gd, label=T, group.by="seurat_clusters", legend.position="none",
           colors.use=pal_seur)
ggsave("../figures/FLASHseq_UMAP_by_Cluster.pdf", width=2.7, height=3)
ggsave("../figures/FLASHseq_UMAP_by_Cluster.png", width=2.7, height=3, dpi=600)


##-----------------------------------------------------------------------------
# 1c
do_DimPlot(gd, group.by="disease", legend.position="left")+
  scale_color_manual(values=pal_dis ,labels=names_dis)+
  theme(legend.position=c(.2, .15))
ggsave("../figures/FLASHseq_UMAP_by_disease.pdf", width=2.7, height=3)
ggsave("../figures/FLASHseq_UMAP_by_disease.png", width=2.7, height=3, dpi=600)


##-----------------------------------------------------------------------------
# 1d
Idents(gd) <- "cluster"

genes <- c("ITGAE","ITGA1","CD160","GZMA","ID3","KIT","IRF8","TCF7", "IL7R", 
           "PDCD1", "GZMK","CD5","CD27","CCR7","ZNF683","TNF","IFNG", "EOMES",
           "S1PR1","SELL","FCGR3A", "GNLY", "NKG7","GZMB", "PRF1","TBX21",
           "ZEB2","KLRG1", "ITGB2")

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
        legend.position=c(-.6, .5), legend.key.size=unit(4, "mm"), 
        legend.spacing.y=unit(0.5, "mm"), legend.margin=margin(t=-8)
        )+
  guides(fill=guide_legend(ncol=1, byrow=T))+
  geom_hline(yintercept=.5, linetype="dashed", linewidth=.75)

p/p_1+plot_layout(heights = c(1, 0.05))

ggsave("../figures/FLASHseq_DotPlot_1.pdf", width=2.8, height=6.75)
ggsave("../figures/FLASHseq_DotPlot_1.png", width=2.8, height=6.75, dpi=600)

rm(p, p_1, genes)


##-----------------------------------------------------------------------------
# 1e+1f
p <- VlnPlot_scCustom(gd, colors_use=pal_clu, pt.size=0, sort=T,
        features=c("cytokine_cytotoxic_signature1", "TRM_signature1"))&
  geom_boxplot(width=.3, fill="white", outlier.shape=NA, coef=0, color="black")&
  NoLegend()&
  scale_y_continuous(expand=expansion(c(0.05, 0.05)))&
  labs(y="Module Score", x=NULL)&
  theme(plot.title=element_text(hjust=0), axis.text.x=element_text(angle=0, hjust=.5))

for (i in seq_along(p)){
  p[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- F
}

p[[1]]+ggtitle(label="Cytotoxicity and Cytokine\nProduction Signature",
               subtitle=expression(italic("Szabo et al., 2019")))

ggsave("../figures/FLASHseq_VlnPlot_Cyto.pdf", width=3.5, height=3.5)
ggsave("../figures/FLASHseq_VlnPlot_Cyto.png", width=4, height=3.5, dpi=600)

p[[2]]+ggtitle(label="Tissue Residency\nSignature",
               subtitle=expression(italic("Kumar et al., 2017")))

ggsave("../figures/FLASHseq_VlnPlot_TRM.pdf", width=3.5, height=3.5)
ggsave("../figures/FLASHseq_VlnPlot_TRM.png", width=4, height=3.5, dpi=600)

rm(i, p)


##-----------------------------------------------------------------------------
# 1g+1h
do_DimPlot(gd, group.by = "cluster",pt.size = 2, legend.position = "right",
           reduction="scenic_umap", colors.use=pal_clu)+
  guides(color=guide_legend(ncol=2, override.aes=list(size=4)))+
  theme(legend.position=c(0.85, 0.9), text=element_text(size=24))

ggsave("../figures/FLASHseq_SCENIC_UMAP_byCluster.pdf", width=6, height=5)
ggsave("../figures/FLASHseq_SCENIC_UMAP_byCluster.png", dpi=600, width=6, height=5)

do_DimPlot(gd, group.by = "disease",pt.size = 2, legend.position = "right",
           reduction="scenic_umap")+
  scale_color_manual(values=pal_dis, labels=names_dis)+
  theme(legend.position=c(0.85, 0.9), text=element_text(size=24))

ggsave("../figures/FLASHseq_SCENIC_UMAP_byDisease.pdf", width=6, height=5)
ggsave("../figures/FLASHseq_SCENIC_UMAP_byDisease.png", dpi=600, width=6, height=5)


##-----------------------------------------------------------------------------
# 1i
# Import AUC scores
rss <- importAUCfromText(
  "../../data_processed/FLASHseq/pyscenic_output.csv")
rss@NAMES <- gsub(rss@NAMES, pattern="\\(\\+\\)", replacement="")

rss <- calcRSS(AUC=getAUC(rss),cellAnnotation=gd@meta.data$cluster)
rss <- plotRSS(rss)
rss <- rss$df

# DotPlot
regulons <- c("TBX21", "EOMES", "PRDM1", "STAT1", "STAT5B", "BCL3", "BCL6", "TCF7", 
              "LEF1", "KLF7", "JUN", "JUNB", "JUND", "FOS", "FOSB")

rss%>%
  filter(Topic %in% regulons)%>%
  mutate(Topic=factor(Topic, levels=rev(regulons)),
         cellType=factor(cellType, levels=paste0("C", 0:3)),
         Z=case_when(Z>4~4, .default=Z)
         )%>%
  
ggplot(aes(x=cellType, y=Topic, size=RSS, fill=Z))+
  geom_point(color="black", shape=21)+
  theme_minimal()+
  scale_fill_gradientn(colors = c("white", "navy"), 
    guide=guide_colorbar(frame.colour="black", frame.linewidth=.5,
      ticks.colour="black", ticks.linewidth=.5)
                      )+
  labs(fill="Z-Score")+
  theme(axis.text=element_text(color="black"),axis.title=element_blank(),
        axis.text.x=element_text(size=11))+
  RotatedAxis()

ggsave("../figures/FLASHseq_SCENIC_RSS_by_Cluster.pdf", width=2.5, height=4)
ggsave("../figures/FLASHseq_SCENIC_RSS_by_Cluster.png", dpi=600 ,width=2.5, height=4)

rm(rss, regulons)


# Figure 3

##-----------------------------------------------------------------------------
# 3a
do_DimPlot(gd, group.by="TRDV",
           idents.keep=unique(gd$TRDV[!is.na(gd$TRDV)]),na.value="gray80")+
  scale_color_manual(values=pal_trdv)+
  guides(color=guide_legend(nrow=3, override.aes=list(size=4)))+
  theme(legend.margin=margin(t=-35), legend.justification="left",
        legend.key.spacing.y=unit(.25, "mm"))
ggsave("../figures/FLASHseq_UMAP_by_TRDV.pdf", width=3.5, height=3)
ggsave("../figures/FLASHseq_UMAP_by_TRDV.png", width=3.5, height=3, dpi=600)

# 3d
do_DimPlot(gd, group.by="TRGV",
           idents.keep=unique(gd$TRGV[!is.na(gd$TRGV)]),na.value="gray80")+
  scale_color_manual(values=pal_trgv)+
  guides(color=guide_legend(nrow=3, override.aes=list(size=4)))+
  theme(legend.margin=margin(t=-25), legend.justification="left",
        legend.key.spacing.y=unit(.25, "mm"))
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
p <- do_DimPlot(gd, group.by="TRGV_TRDV", idents.keep=arranged$TRGV_TRDV[2:6],na.value="gray80")+
  scale_color_manual(values=pal_gvdv, labels=names_gvdv)+
  guides(color=guide_legend(nrow=3, override.aes=list(size=4)))+
  theme(legend.margin=margin(t=-25), legend.justification="left",
        legend.key.spacing.y=unit(.25, "mm"))
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

ggplot(gd@meta.data, aes(x=cluster, fill=fct_rev(TRDV)))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_trdv, na.value="gray80")+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_text(vjust=2), axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black"), legend.margin=margin(l=-7),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"))
ggsave("../figures/FLASHseq_BarPlot_TRDV_by_Cluster.pdf", width=3, height=3)
ggsave("../figures/FLASHseq_BarPlot_TRDV_by_Cluster.png", width=3, height=3, dpi=600)


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

ggplot(gd@meta.data, aes(x=cluster, fill=fct_rev(TRGV)))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_trgv, na.value="gray80")+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_text(vjust=2),
        axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"),
        legend.margin=margin(l=-7))
ggsave("../figures/FLASHseq_BarPlot_TRGV_by_Cluster.pdf", width=2.5, height=3.25)
ggsave("../figures/FLASHseq_BarPlot_TRGV_by_Cluster.png", width=2.5, height=3.25, dpi=600)


##-----------------------------------------------------------------------------
# 3h
gd@meta.data%>%
  mutate(TRGV_TRDV=case_when(
    TRGV_TRDV %in% levels(arranged$TRGV_TRDV)[1:5] ~ TRGV_TRDV,
    TRGV_TRDV %in% levels(arranged$TRGV_TRDV[-1:-5]) ~ "other"
  ),
          TRGV_TRDV=factor(TRGV_TRDV, levels=c(levels(arranged$TRGV_TRDV[1:5]), "other"))        
  )%>%
ggplot(aes(x=cluster, fill=TRGV_TRDV))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_gvdv, labels=names_gvdv, na.value="gray80")+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_text(vjust=2), axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=12), legend.margin=margin(l=-7),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"))
ggsave("../figures/FLASHseq_BarPlot_TRDV-TRGV_by_Cluster.pdf", width=3.25, height=4)
ggsave("../figures/FLASHseq_BarPlot_TRDV-TRGV_by_Cluster.png", width=3.25, height=4, dpi=600)


##-----------------------------------------------------------------------------
# # Test Function
# test <-   function(y.pos=NULL, x.pos=NULL){
#   stat_compare_means(bracket.size=0, method = "t.test", label.y=y.pos, label.x=x.pos,
#   aes(label = ifelse(..p.. < 0.0001, "p<0.0001", paste0("p=", ..p.format..)))
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
  stat_compare_means(bracket.size=0, method = "t.test", label.x=1.27,
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
  stat_compare_means(bracket.size=0, method = "t.test", label.x=1.27,
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
gd@meta.data%>%
  filter(!is.na(TRGV_TRDV))%>%
  group_by(disease_ID)%>%
  count(TRGV_TRDV)%>%
  mutate(percent=(n/sum(n)*100))%>%
  separate(disease_ID, into=c("disease", "treatment", "id"), sep="_")%>%
  filter(TRGV_TRDV %in% c("TRGV9_TRDV2", "TRGV4_TRDV1"))%>%
  mutate(TRGV_TRDV=factor(TRGV_TRDV, levels=c("TRGV9_TRDV2", "TRGV4_TRDV1")))%>%
ggplot(aes(x=disease, y=percent, fill=disease))+
  geom_boxplot(color="black", outlier.shape=NA)+
  geom_point()+
  scale_fill_manual(values=pal_dis)+
  scale_x_discrete(labels=names_dis)+
  theme_minimal()+
  labs(x=NULL, y="% of all subsets")+
  theme(axis.line.y=element_line(), axis.ticks.y=element_line(),
        panel.grid=element_blank(), axis.text=element_text(color="black"),
        legend.position="none", strip.text=element_text(face="bold", size=12))+
  facet_wrap(~TRGV_TRDV, scales="free_y")

ggsave("../figures/FLASHseq_BarPlot_TRGVs-TRDVs_by_Disease.pdf", width=3.5, height=3)
ggsave("../figures/FLASHseq_BarPlot_TRGVs-TRDVs_by_Disease.png", width=3.5, height=3, dpi=600)


##-----------------------------------------------------------------------------
# 3j
chord.df <- getCirclize(gd, "strict", group.by="cluster",proportion=F)%>%
              arrange(from)

pdf("../figures/FLASHseq_CloneOverlap_byCluster.pdf")

par(cex=2.5, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack = c("name", "grid"),
             self.link=1, grid.col=pal_clu)
title("Clonal Overlap")

dev.off()
 #
png("../figures/FLASHseq_CloneOverlap_byCluster.png", res=600, width=4200, height=4200)

par(cex=2.5, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack = c("name", "grid"),
             self.link=1, grid.col=pal_clu)
title("Clonal Overlap")

dev.off()

rm(chord.df)


##-----------------------------------------------------------------------------
# 3k
names <- c("1"="TRGV2.TRGJ2.TRGC2;CDR3 (TRDV NA)",
           "2"="TRDV2.TRDD3.TRDJ3.TRDC;CDR3 (TRGV NA)",
           "3"="TRGV9.TRGJP.TRGC1;T;CDR3-\nTRDV2.TRDD3.TRDJ1.TRDC;CDR3")

colors <- c("1"="red", "2"="cornflowerblue", "3"="green")

gd@meta.data <- gd@meta.data%>%
  mutate(clonotype_selected=case_when(
    str_detect(CTaa, "CACDTRDPSSWDTRQMFF") ~ "1",
    str_detect(CTaa, "CATWDEHYYKKLF") ~ "2",
    str_detect(CTaa, "CALWEVQGLGKKIKVF_CACDTMGDTDKLIF") ~ "3",
    .default = NA_character_
  ))

gd <- SetIdent(gd, value="clonotype_selected")

do_DimPlot(gd, group.by = "clonotype_selected", pt.size=1.5, na.value="gray80",
          idents.keep=c("1", "2", "3"), legend.position="right")+
  scale_color_manual(values=colors, labels=names)+
  guides(color=guide_legend(nrow=3, byrow=T, override.aes=list(size=3),
                            title="Clones"))+
  theme(plot.background=element_blank(),
        legend.margin=margin(l=-40), legend.text=element_text(face="italic", size=9))

ggsave("../figures/FLASHseq_UMAP_by_C2C3_Clones.pdf", width=5, height=3)
ggsave("../figures/FLASHseq_UMAP_by_C2C3_Clones.png", dpi=600, width=5, height=3)


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
         signif=factor(signif, levels = c("ns", "down", "up"))
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
  scale_color_manual(values = c("gray", "#008B8B", "#8B0000"))+
  guides(color=guide_legend(override.aes=list(size=2), reverse=T))+
  theme(axis.text=element_text(color="black"), text=element_text(size=16),
        legend.position=c(0.15, 0.85))+
  
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
  scale_fill_gradientn(limits = c(1.3, 8.5), colors=c("#E8FFFF", "#00007A"))+
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
  annotate(geom="text", x=32000, y=.7, label=paste0("p = ",
    gsea@result%>%filter(ID=="R-HSA-389948")%>%pull(p.adjust)%>%signif(2)),
    size=5)+
  theme(plot.title=element_text(face="bold"))

ggsave("../figures/FLASHseq_GSEAplot_PD1.pdf", width=7, height=4.5)
ggsave("../figures/FLASHseq_GSEAplot_PD1.png", width=7, height=4.5, dpi=600)

#
gseaplot(gsea, geneSetID="R-HSA-201722", by="runningScore",
         title=gsea@result%>%filter(ID=="R-HSA-201722")%>%pull(Description))+
  annotate(geom="text", x=32000, y=.45, label=paste0("p = ",
    gsea@result%>%filter(ID=="R-HSA-201722")%>%pull(p.adjust)%>%signif(2)),
    size=5)+
  theme(plot.title=element_text(face="bold"))

ggsave("../figures/FLASHseq_GSEAplot_bCat_TCF.pdf", width=7, height=4.5)
ggsave("../figures/FLASHseq_GSEAplot_bCat_TCF.png", width=7, height=4.5, dpi=600)

# rm(gsea)


##-----------------------------------------------------------------------------
# 4d, 4e
heat <- DotPlot(gd, group.by="cluster",
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
        legend.justification="top")


p2 <- heat%>%filter(features.plot %in% c("PDCD1", "LAG3", "CD244", "CTLA4", "HAVCR2"))%>%
  ggplot(aes(x=id, features.plot,fill=avg.exp.scaled))+
  geom_tile(color="white", size=.25)+
  theme_minimal()+
  scale_fill_gradient2(low="blue", mid="beige", high="red", midpoint=0)+
  labs(fill="Avg.\nExpr.\nScal.", title="Inhibitory\nReceptors")+
  theme(panel.grid.major=element_blank(), 
        panel.border=element_rect(color="black", fill=NA, size=1),
        axis.title=element_blank(), axis.text=element_text(color="black"))

ggarrange(p1, p2, common.legend=T, legend="right", widths=c(1, .95))

ggsave("../figures/FLASHseq_HeatMap_1.pdf", width=4.3, height=2)
ggsave("../figures/FLASHseq_HeatMap_1.png", width=4.3, height=2, dpi=600)

rm(p1, p2, heat)


##-----------------------------------------------------------------------------
# 4f-4j
scores <- grep(colnames(gd@meta.data), pattern="stem1$|effector1$", value=T)

p <- VlnPlot_scCustom(gd, features=scores, group.by="cluster",
                 pt.size=0, colors_use=pal_clu, sort=T,
                 add.noise=F, num_columns=2)&
  theme(axis.title=element_blank(),axis.title.y=element_text(angle=0, vjust=0.5))&
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.05)))&
  geom_boxplot(width=.4, fill="white", outlier.shape=NA, coef=0, color="black", size=.35)

p[[1]] <- p[[1]]+ggtitle("Stem-Like\nSignature")
p[[2]] <- p[[2]]+ggtitle("Effector\nSignature")

for (i in 3:length(p)) {
  p[[i]] <- p[[i]] + ggtitle(NULL)
}


p[[1]] <- p[[1]] + ylab(bquote(atop("Siddiqui"~italic("et al."), "2019")))
p[[3]] <- p[[3]] + ylab(bquote(atop(plain("Im")~italic("et al."), "2016")))
p[[5]] <- p[[5]] + ylab(bquote(atop(plain("Wu")~italic("et al."), "2016")))
p[[7]] <- p[[7]] + ylab(bquote(atop(plain("Utzschneider"), italic("et al.")~"2016")))
p[[9]] <- p[[9]] + ylab(bquote(atop(plain("Li")~italic("et al."), "2024")))


for (i in seq_along(p)) {
  p[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE
}

p

ggsave("../figures/FLASHseq_VlnPlots_Stem_Effector_Signatures.pdf", width=7, height=10)
ggsave("../figures/FLASHseq_VlnPlots_Stem_Effector_Signatures.png", width=7, height=10)

rm(p, i, scores)


##-----------------------------------------------------------------------------
# 4k, 4l
heat <- DotPlot(gd, group.by="cluster",
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
        legend.justification="top")

heat <- DotPlot(gd, group.by="cluster",
    features=c("TOX", "ID3", "ID2", "FOXO1", "TCF7", "BCL6", "EOMES", "TBX21"))$data

p2 <- ggplot(heat, aes(x=id, y=features.plot,fill=avg.exp.scaled))+
  geom_tile(color="white", size=.25)+
  theme_minimal()+
  scale_fill_gradient2(low="blue", mid="beige", high="red", midpoint=0)+
  guides(fill=guide_colorbar(frame.colour="black", frame.linewidth=.5,
                             ticks.colour="black", ticks.linewidth=.25
                             ))+
  labs(fill="Avg.\nExpr.\nScal.", title="Transcription\nFactors")+
  theme(panel.grid.major=element_blank(), panel.border=element_rect(color="black", fill=NA, size=1),
        axis.title=element_blank(), axis.text=element_text(color="black"))

ggarrange(ggarrange(p1, NULL, heights=c(1, 0.2), nrow=2), 
          p2, common.legend=T, legend="right", widths=c(1, .95))

ggsave("../figures/FLASHseq_HeatMap_2.pdf", width=4, height=2.5)
ggsave("../figures/FLASHseq_HeatMap_2.png", width=4, height=2.5, dpi=600)

rm(p1, p2, heat)


# Supplements & New Figures
##-----------------------------------------------------------------------------
do_DimPlot(gd, split.by="ID")&NoLegend()
ggsave("../figures/FLASHseq_SplitBy_Donor.pdf", width=10, height=7)
ggsave("../figures/FLASHseq_SplitBy_Donor.png", width=10, height=7, dpi=600)


##-----------------------------------------------------------------------------
FeaturePlot_scCustom(gd, features=c("TRDV1", "TRDV2", "TRGV4", "TRGV9"), na_color="#2A0134",
  colors_use=c("#6900A8FF", "#B02991FF", "#E16462FF", "#FCA636FF", "#F0F921FF"))&
  NoLegend()&NoAxes()

ggsave("../figures/FLASHseq_FeaturePlot_TRgenes.pdf", width=4.5, height=4.5)
ggsave("../figures/FLASHseq_FeaturePlot_TRgenes.png", width=4.5, height=4.5, dpi=600)


##-----------------------------------------------------------------------------
FeaturePlot_scCustom(gd, features=c("TCF7", "PDCD1"), na_color="#2A0134",
  colors_use=c("#6900A8FF", "#B02991FF", "#E16462FF", "#FCA636FF", "#F0F921FF"))&
  NoLegend()&NoAxes()

ggsave("../figures/FLASHseq_FeaturePlot_StemLikeGenes.pdf", width=4.5, height=2.75)
ggsave("../figures/FLASHseq_FeaturePlot_StemLikeGenes.png", width=4.5, height=2.75, dpi=600)


##-----------------------------------------------------------------------------
# 
chord.df <- getCirclize(gd, "strict", group.by="ID",proportion=F)%>%
              arrange(from)

pdf("../figures/FLASHseq_CloneOverlap_byDonor.pdf")

par(cex=2.5, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack = c("name", "grid"),
             self.link=1, grid.col=SCpubr:::generate_color_scale(sort(unique(gd$ID))))
title("Clonal Overlap")

dev.off()
 #
png("../figures/FLASHseq_CloneOverlap_byDonor.png", res=600, width=4200, height=4200)

par(cex=2.5, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack = c("name", "grid"),
             self.link=1, grid.col=SCpubr:::generate_color_scale(sort(unique(gd$ID))))
title("Clonal Overlap")

dev.off()

rm(chord.df)


##-----------------------------------------------------------------------------
sessionInfo()



















