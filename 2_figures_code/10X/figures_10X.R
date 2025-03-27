---
title: "PBMC scRNA-seq (10X): Code for Figures"
output: github_document
---




##-----------------------------------------------------------------------------
library(Seurat)
library(SCpubr)
library(scCustomize)
library(tidyverse)
library(scRepertoire)
library(circlize)
library(patchwork)


# Data
##-----------------------------------------------------------------------------
vd1 <- readRDS("../../data_processed/10X/Vd1.Rds")
vd2 <- readRDS("../../data_processed/10X/Vd2.Rds")


##-----------------------------------------------------------------------------
FindAllMarkers(vd1, only.pos=T)%>%
  filter(p_val_adj<0.05)%>%
  write.csv("Vd1_All_Markers_byCluster.csv")

##-----------------------------------------------------------------------------
FindAllMarkers(vd2, only.pos=T)%>%
  filter(p_val_adj<0.05)%>%
  write.csv("Vd2_All_Markers_byCluster.csv")



# Parameters
##-----------------------------------------------------------------------------
pal_clu <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gray70", "gold1",
  "deeppink1", "#FB9A99", "darkorange4"
            )
pal_SeurClu <- setNames(pal_clu, paste0(0:9))
pal_clu <- setNames(pal_clu, paste0("C",0:9))

pal_dis <- c("HD"="#008B8B","UC"="#8B0000")
names_dis <- c("UC"="UC (A)")

pal_trgv <- c("#BCBD22FF", "#7F7F7FFF", "#E377C2FF", "#8C564BFF", "yellow", "#9467BDFF",
    "#D62728FF", "#2CA02CFF", "#FF7F0EFF", "#1F77B4FF")
pal_trgv <- setNames(pal_trgv, levels(vd1$TRGV))


# Vd1
##-----------------------------------------------------------------------------
p <- do_DimPlot(vd1, label=T, group.by="seurat_clusters", legend.position="none",
           colors.use=pal_SeurClu, pt.size=0.2)

p[[1]][["layers"]][[3]][["data"]] <- p[[1]][["layers"]][[3]][["data"]]%>%
  mutate(umap_1=case_when(seurat_clusters=="8" ~ 1.5, .default=umap_1))
p

ggsave("../figures/10X_Vd1_UMAP_by_Cluster.pdf", width=4.5, height=4.5)
ggsave("../figures/10X_Vd1_UMAP_by_Cluster.png", dpi=600, bg="white", width=4.5, height=4.5)

do_DimPlot(vd1,group.by="disease",colors.use=pal_dis, pt.size=0.1)+
  scale_color_manual(values=pal_dis, labels=names_dis)+
  theme(legend.position=c(.2, .8))

ggsave("../figures/10X_Vd1_UMAP_by_Disease.pdf", width=3, height=3)
ggsave("../figures/10X_Vd1_UMAP_by_Disease.png", dpi=600, bg="white", width=3, height=3)


##-----------------------------------------------------------------------------
clonalOccupy(vd1, x.axis="seurat_clusters", label=F)&
  scale_fill_manual(values=c(
      "#fffe9e","#f79522", "#c43271", "#611062", "gray30"),
                    labels=c("Hyper (>100)", "Large (26-100)",
                    "Medium (5-25)", "Rare (2-5)","Single Clone"),
                    name="Clone Size")&
  labs(y="Number of Cells", x="Cluster")&
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
        legend.position=c(.8, .7), axis.title.x=element_text(),
        legend.key.size=unit(4, "mm"), legend.key.spacing.y=unit(.5, "mm"))

ggsave("../figures/10X_Vd1_CloneSize_byCluster.pdf", width=4, height=3)
ggsave("../figures/10X_Vd1_CloneSize_byCluster.png", dpi=600, bg="white", width=4, height=3)


##-----------------------------------------------------------------------------
chord.df <- getCirclize(vd1, "strict", group.by="seurat_clusters", proportion=F)%>%
              arrange(from)


##-----------------------------------------------------------------------------
pdf("../figures/10X_Vd1_CloneOverlap_byCluster.pdf")

par(cex=2.5, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack = c("name", "grid"),
             self.link=1, grid.col=pal_SeurClu)
title("Clonal Overlap")

dev.off()


png("../figures/10X_Vd1_CloneOverlap_byCluster.png", res=600, width=4200, height=4200)

par(cex=2.5, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack = c("name", "grid"),
             self.link=1, grid.col=pal_SeurClu)
title("Clonal Overlap")

dev.off()

##-----------------------------------------------------------------------------
chord.df <- getCirclize(vd1, "strict", group.by="diseaseID",proportion=F)%>%
              arrange(from)


##-----------------------------------------------------------------------------
pdf("../figures/10X_Vd1_CloneOverlap_byDonor.pdf")

par(cex=1.8, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack = c("name", "grid"),
             self.link=1, grid.col=setNames(scCustomize_Palette(length(unique(chord.df$from))), unique(chord.df$from)))
title("Clonal Overlap")

dev.off()


##-----------------------------------------------------------------------------
vd1$cluster <- factor(vd1$cluster, levels=c(
  "C8","C0","C1","C3","C6","C2","C4","C5", "C7"
))

genes <- c("GZMB", "PRF1", "TBX21", "CD74", "HLA-DRB1", "HLA-DRA", "ZNF683", "SOX4",
           "TCF7", "LEF1", "IL7R", "CD27", "CCR7", "SELL", "MKI67")

p <- DotPlot(vd1, group.by="cluster", features=genes)+
  scale_color_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  guides(color=guide_colorbar(frame.colour="black", frame.linewidth=.5, order=1,
          ticks.colour="black", ticks.linewidth=.5, title="Avg.\nScal.\nExpr."),
         size=guide_legend(override.aes=list(fill="white", shape=21),title="% Expr.",
                           order=2)
        )+
  geom_point(aes(size=pct.exp+12), shape=21, fill=NA, colour="black", stroke=.25)+
  theme(axis.title=element_blank(),axis.ticks.length.x=unit(1.5, "mm"),
        panel.grid.major=element_line(color="gray95"), legend.justification="top")+
  coord_flip()+RotatedAxis()


p_1 <- vd1@meta.data%>%mutate(cluster=factor(cluster, levels=levels(p$data$id)))

p_1 <- ggplot(p_1, aes(x=cluster, fill=disease))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_dis, labels=names_dis)+
  theme_void()+
  theme(text=element_text(size=18),legend.title=element_blank(),
        legend.position=c(1.3, .5), legend.key.size=unit(4, "mm"), 
        legend.spacing.y=unit(0.5, "mm"), legend.margin=margin(t=-8)
        )+
  guides(fill=guide_legend(ncol=1, byrow=T))+
  geom_hline(yintercept=.5, linetype="dashed", linewidth=.75)

p/p_1+plot_layout(heights = c(1, 0.1))

ggsave("../figures/10X_Vd1_DotPlot_byCluster.pdf", width=4.1, height=4.75)
ggsave("../figures/10X_Vd1_DotPlot_byCluster.png", dpi=600, bg="white", width=4.1, height=4.75)

rm(p, p_1, genes)


##-----------------------------------------------------------------------------
ggplot(vd1@meta.data, aes(x=cluster, fill=disease))+
  geom_bar(position="fill", color="white", size=.5)+
  scale_fill_manual(values=c(rep("#008B8B", 5),
                    rep("#8B0000", 7)))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(), legend.title=element_blank(),
        axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text=element_text(color="black"), axis.ticks.length=unit(2, "mm"),
        axis.ticks=element_line(color="black"), axis.text.x=element_text(vjust=3),
        legend.position="none")


ggsave("../figures/10X_Vd1_BarPlot_ClusterByDonor.pdf", width=3.25, height=3.75)
ggsave("../figures/10X_Vd1_BarPlot_ClusterByDonor.png",dpi=600, width=3.25, height=3.75)


##-----------------------------------------------------------------------------
ggplot(vd1@meta.data, aes(x=cluster, fill=fct_rev(TRGV)))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_trgv, na.value="gray80")+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(), legend.title=element_blank(),
        axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text=element_text(color="black"), axis.text.x=element_text(vjust=3),
        axis.ticks.length=unit(2, "mm"), axis.ticks=element_line(color="black"),
        legend.key.spacing.y = unit(1, "mm"), legend.margin=margin(l=-10))

ggsave("../figures/10X_Vd1_BarPlot_ClusterByTRGV.pdf", width=3.75, height=3.75)
ggsave("../figures/10X_Vd1_BarPlot_ClusterByTRGV.png", width=3.75, height=3.75, dpi=600)




# Vd2
##-----------------------------------------------------------------------------
do_DimPlot(vd2, label=T, group.by="seurat_clusters", legend.position="none",
           colors.use=pal_SeurClu, pt.size=0.2)

ggsave("../figures/10X_Vd2_UMAP_by_Cluster.pdf", width=4.5, height=4.5)
ggsave("../figures/10X_Vd2_UMAP_by_Cluster.png", dpi=600, bg="white", width=4.5, height=4.5)

do_DimPlot(vd2,group.by="disease",colors.use=pal_dis, pt.size=0.1)+
  scale_color_manual(values=pal_dis, labels=names_dis)+
  theme(legend.position=c(.15, 1))

ggsave("../figures/10X_Vd2_UMAP_by_Disease.pdf", width=3, height=3)
ggsave("../figures/10X_Vd2_UMAP_by_Disease.png", dpi=600, bg="white", width=3, height=3)


##-----------------------------------------------------------------------------
clonalOccupy(vd2, x.axis="cluster", label=F)&
  scale_fill_manual(values=c(
      "#fffe9e","#f79522", "#c43271", "#611062", "gray30"),
                    labels=c("Hyper (>100)", "Large (26-100)",
                    "Medium (5-25)", "Rare (2-5)","Single Clone"),
                    name="Clone Size")&
  labs(y="Number of Cells", x="Cluster")&
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"),
        legend.position=c(.8, .7), axis.title.x=element_text(),
        legend.key.size=unit(4, "mm"), legend.key.spacing.y=unit(.5, "mm"))

ggsave("../figures/10X_Vd2_CloneSize_byCluster.pdf", width=4, height=3)
ggsave("../figures/10X_Vd2_CloneSize_byCluster.png", dpi=600, bg="white", width=4, height=3)


##-----------------------------------------------------------------------------
chord.df <- getCirclize(vd2, "strict", group.by="seurat_clusters",proportion=F)%>%
              arrange(from)


##-----------------------------------------------------------------------------
pdf("../figures/10X_Vd2_CloneOverlap_byCluster.pdf")

par(cex=2.5, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack = c("name", "grid"),
             self.link=1, grid.col=pal_SeurClu)
title("Clonal Overlap")

dev.off()

png("../figures/10X_Vd2_CloneOverlap_byCluster.png", res=600, width=4200, height=4200)

par(cex=2.5, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack = c("name", "grid"),
             self.link=1, grid.col=pal_SeurClu)
title("Clonal Overlap")

dev.off()

##-----------------------------------------------------------------------------
chord.df <- getCirclize(vd2, "strict", group.by="diseaseID",proportion=F)%>%
              arrange(from)


##-----------------------------------------------------------------------------
pdf("../figures/10X_Vd2_CloneOverlap_byDonor.pdf")

par(cex=1.8, mar=c(0, 0 , 2, 0))
chordDiagram(chord.df, annotationTrack = c("name", "grid"),
             self.link=1, grid.col=setNames(
              scCustomize_Palette(length(unique(chord.df$from))), unique(chord.df$from)
               ))
title("Clonal Overlap")

dev.off()


##-----------------------------------------------------------------------------
vd2$cluster <- factor(vd2$cluster,
                      levels=c("C0","C1","C4","C5","C8","C9","C2","C3","C7","C6"))

genes <- c("GZMB", "GNLY", "PRF1", "NKG7", "TCF7", "IL7R", "GZMK", "LEF1", 
           "CD27", "CCR7", "S1PR1", "SELL", "RORC", "IL23R", "CCR6", "CXCR6", 
           "SCART1", "BLK"
           )

p <- DotPlot(vd2, group.by="cluster", features=genes)+
  scale_color_gradient2(low="navy", mid="beige", high="red3",, midpoint=0)+
  guides(color=guide_colorbar(frame.colour="black", frame.linewidth=.5, order=1,
          ticks.colour="black", ticks.linewidth=.5, title="Avg.\nScal.\nExpr."),
         size=guide_legend(override.aes=list(fill="white", shape=21),title="% Expr.",
                           order=2)
        )+
  geom_point(aes(size=pct.exp+12), shape=21, fill=NA, colour="black", stroke=.25)+
  theme(axis.title=element_blank(),axis.ticks.length.x=unit(1.5, "mm"),
        panel.grid.major=element_line(color="gray95"), legend.justification="top")+
  coord_flip()+RotatedAxis()


p_1 <- vd2@meta.data%>%mutate(cluster=factor(cluster, levels=levels(p$data$id)))

p_1 <- ggplot(p_1, aes(x=cluster, fill=disease))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_dis, labels=names_dis)+
  theme_void()+
  theme(text=element_text(size=18),legend.title=element_blank(),
        legend.position=c(1.25, .5), legend.key.size=unit(4, "mm"), 
        legend.spacing.y=unit(0.5, "mm"), legend.margin=margin(t=-8)
        )+
  guides(fill=guide_legend(ncol=1, byrow=T))+
  geom_hline(yintercept=.5, linetype="dashed", linewidth=.75)

p/p_1+plot_layout(heights = c(1, 0.1))

ggsave("../figures/10X_Vd2_DotPlot_byCluster.pdf", width=4.1, height=5.5)
ggsave("../figures/10X_Vd2_DotPlot_byCluster.png", dpi=1600, bg="white", width=4.1, height=5.5)

rm(p, p_1, genes)


##-----------------------------------------------------------------------------
ggplot(vd2@meta.data, aes(x=cluster, fill=diseaseID))+
  geom_bar(position="fill", color="white", size=.5)+
  scale_fill_manual(values=c(rep("#008B8B", 5),
                    rep("#8B0000", 7)))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(), legend.title=element_blank(),
        axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text=element_text(color="black"), axis.ticks.length=unit(2, "mm"),
        axis.ticks=element_line(color="black"), axis.text.x=element_text(vjust=3),
        legend.position="none")

ggsave("../figures/10X_Vd2_BarPlot_ClusterByDonor.pdf", width=3.75, height=3.75)
ggsave("../figures/10X_Vd2_BarPlot_ClusterByDonor.png", dpi=600, width=3.75, height=3.75)


##-----------------------------------------------------------------------------
ggplot(vd2@meta.data, aes(x=cluster, fill=fct_rev(TRGV)))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_trgv, na.value="gray80")+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(), legend.title=element_blank(),
        axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text=element_text(color="black"), axis.text.x=element_text(vjust=3),
        axis.ticks.length=unit(2, "mm"), axis.ticks=element_line(color="black"),
        legend.key.spacing.y = unit(1, "mm"), legend.margin=margin(l=-10))

ggsave("../figures/10X_Vd2_BarPlot_ClusterByTRGV.pdf", width=3.75, height=3.75)
ggsave("../figures/10X_Vd2_BarPlot_ClusterByTRGV.png",dpi=600, width=3.75, height=3.75)


##-----------------------------------------------------------------------------
sessionInfo()








































