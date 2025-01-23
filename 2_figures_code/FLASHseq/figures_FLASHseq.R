##-----------------------------------------------------------------------------
library(Seurat)
library(SCpubr)
library(scCustomize)
library(tidyverse)
library(scRepertoire)
library(circlize)
library(patchwork)
library(ggpubr)
library(rstatix)
library(SCENIC)
library(AUCell)
library(paletteer)


##-----------------------------------------------------------------------------
gd <- readRDS("../../data_processed/FLASHseq/gd.Rds")
Idents(gd) <- "cluster"

##-----------------------------------------------------------------------------
pal_seur <- c("0"="#4682B4", "1"="#AF46B4", "2"="#B47846", "3"="#4BB446")
pal_clust <- c("C0"="#4682B4", "C1"="#AF46B4", "C2"="#B47846", "C3"="#4BB446")
pal_dis <- c("control"="#008B8B","UC"="#8B0000")



# Med2 Retreat Poster
## Fig2

A
##-----------------------------------------------------------------------------
A <- do_DimPlot(gd, label = T, group.by = "seurat_clusters", label.size = 4,
                 pt.size = 1.5, legend.position = "none", colors.use=pal_seur)+
  annotate("segment", size = 1, lineend = "round",
    x = -3.75, xend = c(0.15, -3.75),
    y = -5, yend = c(-5, -1.9),
    arrow = arrow(type = "closed", length = unit(10, 'pt')))+
  annotate("text", size = 5,
    x = c(-2., -4.4), y = c(-5.5, -3.8),
    label = c("UMAP1", "UMAP2"), angle = c(0, 90))
A <- A[[1]]


B
##-----------------------------------------------------------------------------
B <-do_DimPlot(gd, group.by = "disease", pt.size = 1.5, colors.use=pal_dis)+
  theme(legend.position = c(0.33, .097), legend.background = element_blank(),
        legend.text = element_text(size=14))+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(labels=c("control"="ctrl"), values=pal_dis)
B <- B[[1]]


C
##-----------------------------------------------------------------------------
# test <- gd@meta.data%>%wilcox_test(cytokine_cytotoxic_signature1~cluster)%>%
#   mutate(xmin = c(1, 1, 1, 2.05, 2, 3.05),
#          xmax = c(1.95, 3, 4, 2.95, 4, 4 ),
#          y.position = c(0.75, 1.1, 1.45, 0.75, 1.85, 0.75))%>%
#     mutate(p.adj.signif = case_when(p.adj.signif == "ns" ~ "", 
#                                     p.adj.signif == "****" ~ "***",
#                                     .default = p.adj.signif))

gd <- SetIdent(gd, value="cluster")

C <- VlnPlot_scCustom(gd, features = "cytokine_cytotoxic_signature1", 
                       sort = "increasing", colors_use=pal_clust, pt.size = 0)+
  theme(text = element_text(size = 14), axis.text = element_text(size = 14),
        plot.title = element_text(size=16, face="plain"))+
  ggtitle("Cytotoxic and\nCytokine")+
  #ylim(-0.2, 1.8)+
  xlab(element_blank())+ylab(element_blank())+
  NoLegend()+
  scale_y_continuous(expand = expansion(c(0.05, 0.05)))+
  
 # stat_pvalue_manual(test, inherit.aes=F, size = 5, bracket.size = .75, vjust = 0.2)+
  #annotate("text", x = 2.5, y = 0.87, label = "ns", size = 4) 
  geom_boxplot(fill="white", color="black", width=.3, outlier.shape = NA, coef=0)

C[[1]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE

D
##-----------------------------------------------------------------------------
# test <- gd@meta.data%>%wilcox_test(TRM_signature1~seurat_clusters)%>%
#   mutate(xmin = c(1, 1, 1, 2.05, 2, 3.05),
#          xmax = c(1.95, 3, 4, 2.95, 4, 4 ),
#          y.position = c(0.25, 0.4, 0.55, 0.25, 0.7, 0.25))%>%
#     mutate(p.adj.signif = case_when(p.adj.signif == "ns" ~ "", 
#                                     p.adj.signif == "****" ~ "***",
#                                     .default = p.adj.signif))

D <- VlnPlot_scCustom(gd, features = "TRM_signature1", colors_use = c(
                        "#AF46B4", "#4682B4", "#B47846", "#4BB446"),
                       pt.size = 0, sort = "increasing")+
  theme(text = element_text(size = 14), axis.text = element_text(size = 14),
        plot.title = element_text(size = 16, face="plain"))+
  ggtitle("Tissue\nResidency")+
  #ylim(-0.2, .8)+
  xlab(element_blank())+ylab(NULL)+
  NoLegend()+
  scale_y_continuous(expand = expansion(c(0.05, 0.05)))+
  
  #stat_pvalue_manual(test, inherit.aes=F, size = 5, bracket.size = .75, vjust = 0.17)
   geom_boxplot(fill="white", color="black", width=.3, outlier.shape = NA, coef=0)
   
D[[1]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE


E
##-----------------------------------------------------------------------------
gd <- SetIdent(gd, value="cluster")

genes <- c( "ITGAE","ITGA1","CD160","GZMA","ID3","KIT","IRF8","TCF7",
                      "IL7R", "PDCD1", "GZMK","CD5","CD27","CCR7","ZNF683","TNF","IFNG",
                      "EOMES","S1PR1","SELL","FCGR3A", "GNLY", "NKG7","GZMB", "PRF1","TBX21","ZEB2","KLRG1",
                      "ITGB2")

E <- DotPlot(gd, group.by = "cluster", features = genes, scale = T, cluster.idents = T)+
  scale_color_gradient2(low = "navy", mid = "beige", high = "darkred", midpoint = 0)+
    
guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = .5,
                           ticks.colour = "black", ticks.linewidth = .5,
                           title="Avg.\nScal.\nExpr."),
        size = guide_legend(override.aes = list(fill = "white", shape = 21),
                             title="% Expr."
                             )
       )+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+
  theme(axis.title = element_blank(),
        axis.ticks.length.x = unit(1.5, "mm"),
        panel.grid.major = element_line(color="gray95"))+
  coord_flip()+RotatedAxis()


meta <- gd@meta.data%>%mutate(
  cluster=factor(cluster, levels = 
                              levels(p$data$id)))

E_1 <- ggplot(meta, aes(x=cluster, fill=disease))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_dis)+
  theme_void()+
  theme(text=element_text(size=18),,
        legend.title = element_blank(),
        legend.justification = c("left", "top"),
        legend.key.size = unit(4, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.margin = margin(t=-8),
        plot.background = element_blank())+
  guides(fill=guide_legend(ncol=1, byrow = T))+
  geom_hline(yintercept = .5, linetype="dashed", linewidth=.75)

E <- E/E_1+
  plot_layout(heights = c(1, 0.05))

F
##-----------------------------------------------------------------------------
clonotypes_percent <- gd@meta.data %>%
  filter(!is.na(TRGV_TRDV))%>%
  group_by(ID, disease) %>%
  summarize(TRGV4_TRDV1 = mean(TRGV_TRDV == "TRGV4_TRDV1") * 100,
            TRGV9_TRDV2 = mean(TRGV_TRDV == "TRGV9_TRDV2") * 100)%>%
  ungroup()

f <- ggplot(clonotypes_percent, aes(x=disease, y = TRGV4_TRDV1, fill=disease))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size=2)+
  scale_x_discrete(labels=c("ctrl", "UC"))+
  scale_fill_manual(values = pal_dis)+
  theme(text = element_text(size = 17), plot.title = 
          element_text(hjust=0.5, size =17),
        axis.text = element_text(color="black", size = 16),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(2.5, "mm"),
        axis.title.x = element_blank(),
        panel.grid = element_blank(), panel.background = element_blank()
          )+
  labs(y="% of All Subsets", fill=element_blank())+
  NoLegend()+ggtitle("TRGV4-TRDV1")



G
##-----------------------------------------------------------------------------
G <- ggplot(clonotypes_percent, aes(x=disease, y = TRGV9_TRDV2, fill=disease))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size=2)+
  scale_fill_manual(values = pal_dis)+
  scale_x_discrete(labels=c("ctrl", "UC"))+
  theme(text = element_text(size = 17), plot.title = 
          element_text(hjust=0.5, size = 17),
        axis.text = element_text(color="black", size = 16),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(2.5, "mm"),
        axis.title.x = element_blank(),
        panel.grid = element_blank(), panel.background = element_blank()
  )+
  labs(y="% of All Subsets", fill=element_blank())+
  NoLegend()+ggtitle("TRGV9-TRDV2")

##-----------------------------------------------------------------------------
scores <- c("utz16_stem1","utz16_effector1")

 H <- VlnPlot_scCustom(gd, features=scores, group.by = "cluster",
                 pt.size = 0, colors_use = pal_clust, sort=T,
                 num_columns = 2, adjust=1)&
  theme(axis.title = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(face="plain"))&
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))&
    geom_boxplot(fill="white", color="black", width=.3, outlier.shape = NA, coef=0)



H[[1]] <- H[[1]] + ggtitle("Stem-Like")
H[[2]] <- H[[2]] + ggtitle("Effector")

H <- H + plot_annotation(title = "Gene Signatures from\nUtzschneider et al. 2016",
                    theme = theme(plot.title = element_text(size = 16, hjust=.5)))

for (i in seq_along(H)){
  H[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE
}

##-----------------------------------------------------------------------------
gd@meta.data <- gd@meta.data%>%
  mutate(clonotype_selected=case_when(
    str_detect(CTaa, "CACDTRDPSSWDTRQMFF") ~ "1",
    str_detect(CTaa, "CATWDEHYYKKLF") ~ "2",
    str_detect(CTaa, "CALWEVQGLGKKIKVF_CACDTMGDTDKLIF") ~ "3",
    .default = NA_character_
  ))

names <- c("1"="TRGV2.TRGJ2.TRGC2;CDR3_13aa (TRDV NA)",
           "2"="TRDV2.TRDD3.TRDJ3.TRDC;CDR3_18aa (TRGV NA)",
           "3"="TRGV9.TRGJP.TRGC1;T;CDR3_16aa-\nTRDV2.TRDD3.TRDJ1.TRDC;CDR3_14aa")



gd <- SetIdent(gd, value="clonotype_selected")

I <- do_DimPlot(gd, group.by = "clonotype_selected", pt.size = 1, na.value="gray80",
          idents.keep=c("1", "2", "3"))&
  scale_color_manual(values=c(
    "#1F77B4FF", "green3", "darkorange3"), labels=names)&
  guides(color=guide_legend(nrow=3, byrow=T,
                            override.aes = list(size = 5)))&
  theme(text=element_text(size=12), plot.background = element_blank())


##-----------------------------------------------------------------------------
E
ggsave("../figures/FLASHseq_DotPlot_1.png", dpi=600, width=3.25, height=6.75)



##-----------------------------------------------------------------------------
blank <- ggplot()+theme_void()

ABCD <- ggarrange(A, B, C, D, labels = c("A", "B", "C", "D"), 
                font.label=list(size=24))



ABCDE <- ggarrange(ABCD, blank, E, widths=c(1, 0.065, 0.65), 
                   labels = c("", "E", ""), font.label=list(size=24),
                   nrow=1)

FG <- ggarrange(blank, f, blank, G, labels = c("F","", "G", ""), font.label=list(size=24),
                widths = c(0.21, 1))

i <- ggarrange(blank, I, blank, widths = c(0.4, 1, 0.2), nrow=1,
               font.label=list(size=24), labels=c("", "I",""))

HI <- ggarrange(blank, H,blank, i, labels=c("","H", "", ""), font.label=list(size=24),
                widths = c(0.1, 1))

ABCDEFG <- ggarrange(ABCDE, FG, HI, widths = c(1, 0.27, 0.54),
                     ncol=3)

ABCDEFG

ggsave("Poster_Med2_Fig2.png", dpi=600, bg="white", width = 14, height = 6.5)





Cluster Marker
##-----------------------------------------------------------------------------
gd <- SetIdent(gd, value="cluster")

c3_marker <- FindMarkers(gd, ident.1 ="C3", only.pos=T,
                           logfc.threshold = 1)%>%
  rownames_to_column("gene")%>%filter(p_val_adj<0.05)

c2_marker <- FindMarkers(gd, ident.1 ="C2", only.pos=T,
                           logfc.threshold = 1)%>%
  rownames_to_column("gene")%>%filter(p_val_adj<0.05)

write.csv(c3_marker, "../../genesets/FLASHseq_C3_marker.csv")
write.csv(c2_marker, "../../genesets/FLASHseq_C2_marker.csv")

c2_vs_c3_marker <- FindMarkers(gd_test, ident.1 = "C2", ident.2="C3")%>%
  filter(p_val_adj<0.05&avg_log2FC>0)%>%
  rownames_to_column("gene")

write.csv(c2_vs_c3_marker, "../../genesets/FLASHseq_C2vC3_marker.cs")



# SCENIC
##-----------------------------------------------------------------------------
library(SCENIC)
library(AUCell)
scenic <- importAUCfromText("../../data_processed/FLASHseq/pyscenic_output_hvg_jakob_no_ENSG_nes_2.8.csv")
scenic@NAMES <- gsub(scenic@NAMES, pattern="\\(\\+\\)", replacement="")

rss <- calcRSS(AUC=getAUC(scenic),cellAnnotation=gd@meta.data$cluster)
rss <- plotRSS(rss)
rss <- rss$df


RSS DotPlot
##-----------------------------------------------------------------------------
library(viridis)
regulons <- c("TBX21", "EOMES", "PRDM1", "STAT1", "STAT5B",
                      "BCL3", "BCL6", "TCF7", "LEF1", "KLF7",
                      "JUN", "JUNB", "JUND", "FOS", "FOSB")

rss%>%
  filter(Topic %in% regulons)%>%
  mutate(Topic=factor(Topic, levels=rev(regulons)),
         cellType=factor(cellType, levels=paste0("C", 0:3)),
         Z=case_when(Z>4 ~ 4, .default=Z)
         )%>%
  
ggplot(aes(x=cellType, y=Topic, size=RSS, fill=Z))+
  geom_point(color="black", shape=21)+
  theme_minimal()+
  
  scale_fill_gradientn(colors = c("white", "navy"), 
    guide = guide_colorbar(frame.colour = "black", frame.linewidth = .5,
                           ticks.colour = "black", ticks.linewidth = .5))+
  
  labs(fill = "Z-Score")+
  
  theme(text=element_text(size = 16), axis.text = element_text(color = "black"),
    axis.title = element_blank()
    )+
  RotatedAxis()

ggsave("../figures/FLASHseq_SCENIC_byCLuster.png", dpi=600, bg="white",
       width =3.25, height = 5)



##-----------------------------------------------------------------------------
scenic_mtx <- read.csv("../../data_processed/FLASHseq/pyscenic_output_hvg_jakob_no_ENSG_nes_2.8.csv")

rows <- scenic_mtx$Cell
scenic_mtx$Cell <- NULL
scenic_mtx <- as.matrix(scenic_mtx)
rownames(scenic_mtx) <- rows

scenic_mtx <- t(scenic_mtx)


##-----------------------------------------------------------------------------
scenic <- CreateAssayObject(data = scenic_mtx)

gd[["scenic"]] <- scenic


##-----------------------------------------------------------------------------
DefaultAssay(gd) <- "scenic"

gd <- ScaleData(gd)%>%RunPCA(features = rownames(gd), reduction.key = "scenic",
                             reduction.name="scenic_pca")

##-----------------------------------------------------------------------------
ElbowPlot(gd, 50, reduction = "scenic_pca")


##-----------------------------------------------------------------------------
gd <- RunUMAP(gd, dims=1:10, reduction.name = "umap_scenic", reduction.key = "SCENIC_UMAP", reduction = "scenic_pca")

##-----------------------------------------------------------------------------
do_DimPlot(gd, group.by = "cluster",pt.size = 2, legend.position = "right",
           reduction="umap_scenic", colors.use=pal_clust)+
  guides(color=guide_legend(ncol=2, override.aes=list(size=4)))+
  theme(legend.position=c(0.85, 0.9), text=element_text(size=24))

ggsave("../figures/FLASHseq_SCENIC_UMAP_byCluster.png", dpi=600, width=6, height=5)

do_DimPlot(gd, group.by = "disease",pt.size = 2, legend.position = "right",
           reduction="umap_scenic", colors.use=pal_dis)+
  theme(legend.position=c(0.85, 0.9), text=element_text(size=24))

ggsave("../figures/FLASHseq_SCENIC_UMAP_byDisease.png", dpi=600, width=6, height=5)




##-----------------------------------------------------------------------------
heat_scenic <- DotPlot_scCustom(gd, features = rev(regulons))
heat_scenic <- heat_scenic$data

ggplot(heat_scenic, aes(x=id, y=features.plot, fill=avg.exp.scaled))+
  geom_tile(color="black", size=.3)+
  scale_fill_gradient2(low="navy", mid="beige", high="darkred", midpoint = 0)+
  theme_minimal()+
  labs(fill="Avg. Scal.\nRegulon Act.")+
  theme(panel.grid = element_blank(), axis.title = element_blank())
ggsave("../figures/FLASHseq_SCENIC_Heatmap_Scaled.png", dpi=600, width=3,
       height=4.25, bg="white")



##-----------------------------------------------------------------------------
rss%>%
  filter(Topic %in% regulons)%>%
  mutate(Topic=factor(Topic, levels=rev(regulons)),
         cellType=factor(cellType, levels=paste0("C", 0:3)),
         Z=case_when(Z>3 ~ 3, .default=Z)
         )%>%
ggplot(aes(x=cellType, y=Topic, fill=Z))+
  geom_tile(color="black", size=.3)+
  theme_minimal()+
  scale_fill_gradientn(colors = c("white", "navy"), 
    guide = guide_colorbar(frame.colour = "black", frame.linewidth = .5,
                           ticks.colour = "black", ticks.linewidth = .5))+
  labs(fill = "Z-Score")+
  
  theme(axis.text = element_text(color = "black"),axis.title = element_blank(),
        legend.key.width = unit(4, "mm"))

ggsave("../figures/FLASHseq_SCENIC_Heatmap_Zscore.png", dpi=600, width=2.5,
       height=4, bg="white")


# VDJ
##-----------------------------------------------------------------------------
library(paletteer)
pal_TR <- paletteer_d("ggsci::category10_d3")

gd@meta.data <- gd@meta.data%>%
  mutate(
    TRDV=factor(TRDV, levels=rev(c("TRDV1", "TRDV2", "TRDV3", "TRAV29/DV5", 
                                   "TRAV36/DV7", "TRAV38-2/DV8"))),
    TRGV=factor(TRGV, levels=rev(c(paste0("TRGV", 2:5), paste0("TRGV", 7:11))))
    )


##-----------------------------------------------------------------------------
ggplot(gd@meta.data, aes(x=disease, fill=TRDV))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_TR, na.value="gray")+
  scale_x_discrete(labels=c("control"="Ctrl"))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(), legend.title=element_blank(),
        axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text=element_text(color="black"), axis.ticks.length=unit(2, "mm"),
        legend.key.spacing.y = unit(1, "mm"), legend.margin=margin(l=-10))
ggsave("../figures/FLASHseq_TRDV_by_Disease.png", dpi=600, width=2.75, height=4)


##-----------------------------------------------------------------------------
ggplot(gd@meta.data, aes(x=disease, fill=TRGV))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_TR, na.value="gray")+
  scale_x_discrete(labels=c("control"="Ctrl"))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(), legend.title=element_blank(),
        axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text=element_text(color="black"), axis.ticks.length=unit(2, "mm"),
        legend.key.spacing.y = unit(1, "mm"), legend.margin=margin(l=-10))
ggsave("../figures/FLASHseq_TRGV_by_Disease.png", dpi=600, width=2.75, height=4)

##-----------------------------------------------------------------------------
names <- c("1"="TRGV2.TRGJ2.TRGC2;CDR3_13aa (TRDV NA)",
           "2"="TRDV2.TRDD3.TRDJ3.TRDC;CDR3_18aa (TRGV NA)",
           "3"="TRGV9.TRGJP.TRGC1;T;CDR3_16aa-\nTRDV2.TRDD3.TRDJ1.TRDC;CDR3_14aa")



gd <- SetIdent(gd, value="clonotype_selected")

do_DimPlot(gd, group.by = "clonotype_selected", pt.size=2.5, na.value="gray80", 
           idents.keep=c("1", "2", "3"))&
  scale_color_manual(values=c("#1F77B4FF", "red3", "#4BB446"), labels=names)&
  guides(color=guide_legend(nrow=3, byrow=T,override.aes = list(size=5)))

ggsave("../figures/FLASHseq_UMAP_Clonotypes_Overlap_C2C3.png", dpi=600,
       width=5, height=5)


# MA Plot
##-----------------------------------------------------------------------------
Idents(gd) <- "disease"

UC_markers <- FindMarkers(gd, ident.1 = "UC")%>%
  rownames_to_column("gene")

geneLogSums <- log2(Matrix::rowMeans(GetAssayData(gd, "RNA", "counts")))
geneLogSums <- as.data.frame(geneLogSums)%>%
  rownames_to_column("gene")

UC_markers <- inner_join(UC_markers, geneLogSums, by = "gene")%>%
  mutate(signif = case_when(
           p_val_adj < 0.05 & avg_log2FC > 0 ~ "up",
           p_val_adj < 0.05 & avg_log2FC < 0 ~ "down",
           .default = "ns"
         ),
         signif=factor(signif, levels = c("ns", "down", "up"))
         )



##-----------------------------------------------------------------------------
genes_up <- c("TCF7", "IL7R","GZMK","CD5","CD27","CCR7","ZNF683","TNF",
              "IFNG","EOMES","S1PR1","SELL","NKG7","GZMB","TBX21","KLRG1", "ITGB2", 
              "PDCD1", "CTLA4", "LAG3")

genes_down <- c("ITGAE","ITGA1","CD160","ID3")

ggplot()+
  geom_point(UC_markers%>%filter(signif == "ns"),
             mapping=aes(x=geneLogSums, y=avg_log2FC, color=signif))+
    geom_point(UC_markers%>%filter(signif != "ns"),
             mapping=aes(x=geneLogSums, y=avg_log2FC, color=signif))+
  theme_classic()+
  geom_hline(yintercept=0, linetype="dashed", color="red", size=1)+
  labs(x=expression("log"[2]*"(Mean Expression)"),
       y=expression("log"[2]*"(Fold Change)"),
       color=element_blank())+
  scale_color_manual(values = c("gray", "#008B8B", "#8B0000"))+
  guides(color = guide_legend(override.aes = list(size = 2), reverse=T))+
  theme(axis.text = element_text(color="black"), text=element_text(size=16),
        legend.position = c(0.15, 0.85))+
  
  ggrepel::geom_label_repel(data=UC_markers%>%filter(gene %in% genes_up), 
                           aes(label = gene, x=geneLogSums, y=avg_log2FC),
                           min.segment.length = 0, nudge_y = 0, nudge_x = 0,
                           size = 3, max.overlaps=15)+
  ggrepel::geom_label_repel(data=UC_markers%>%filter(gene %in% genes_down), 
                           aes(label = gene, x=geneLogSums, y=avg_log2FC),
                           min.segment.length = 0, nudge_y = -1, nudge_x = 0,
                           size = 3)

ggsave("../figures/FLASHseq_MA_Plot.png", dpi=600, width=5, height=5)


# Stem-Like Heatmap
##-----------------------------------------------------------------------------
library(biomaRt)
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

genes_inh <- c("Lag3", "Pdcd1", "Ctla4", "Cd244a", "Havcr2")%>%mouse_to_human()
genes_cost <- c("Tnfsf14", "Cd28", "Tnfrsf4", "Icos", "Tnfrsf9")%>%mouse_to_human()
genes_chemo <- c("Cxcl11", "Xcl1", "Cxcl10", "Ccl5", "Ccl4", "Ccl3", "Csf1",
                 "Ccr5", "Cxcr4", "Cxcr5")%>%mouse_to_human()
genes_eff <- c("Fasl", "Ifng", "Tnfsf10", "Prf1", "Gzma", "Gzmb", "Il2", "Tnf")%>%
  mouse_to_human()
genes_tf <- c("Bcl6", "Id3", "Tcf7", "Foxo1", "Id2", "Eomes", "Tbx21", "Tox")%>%
  mouse_to_human()
genes_mem <- c("Sell", "Il7r", "Il2rb", "Klrg1")%>%mouse_to_human("Ccl3")


##-----------------------------------------------------------------------------
h_inh <- DotPlot(gd, features=genes_inh, group.by="cluster")$data
h_cost <- DotPlot(gd, features=genes_cost, group.by="cluster")$data
h_chemo <- DotPlot(gd, features=genes_chemo, group.by="cluster")$data%>%
  filter(avg.exp.scaled!="NaN")
h_chemo$features.plot <- droplevels(h_chemo$features.plot)

h_eff <- DotPlot(gd, features=genes_eff, group.by="cluster")$data
h_tf <- DotPlot(gd, features=genes_tf, group.by="cluster")$data
h_mem <- DotPlot(gd, features=genes_mem, group.by="cluster")$data

##-----------------------------------------------------------------------------
breaks_with_gap <- c(levels(h_cost$features.plot), "", levels(h_inh$features.plot))
ggplot(full_join(h_cost, h_inh), aes(x=id, y=features.plot, fill=avg.exp.scaled))+
  geom_tile()+
  theme_minimal()+
  scale_fill_gradient2(low="blue", mid="beige", high="red", midpoint=0)+
  guides(fill=guide_colorbar(title="Avg.\nExpr.\nScal.", frame.colour="black",
                             ticks.colour="black"))+
  theme(panel.grid.major=element_blank(), axis.title=element_blank())+
  scale_y_discrete(limits = breaks_with_gap, labels = breaks_with_gap)
ggsave("../figures/FLASHseq_HeatMap_1.png", dpi=600, width=3, height=4, bg="white")


##-----------------------------------------------------------------------------
breaks_with_gap <- c(levels(h_eff$features.plot), "", levels(h_chemo$features.plot))
ggplot(full_join(h_eff, h_chemo), aes(x=id, y=features.plot, fill=avg.exp.scaled))+
  geom_tile()+
  theme_minimal()+
  scale_fill_gradient2(low="blue", mid="beige", high="red", midpoint=0)+
  guides(fill=guide_colorbar(title="Avg.\nExpr.\nScal.", frame.colour="black",
                             ticks.colour="black"))+
  theme(panel.grid.major=element_blank(), axis.title=element_blank())+
  scale_y_discrete(limits = breaks_with_gap, labels = breaks_with_gap)
ggsave("../figures/FLASHseq_HeatMap_2.png", dpi=600, width=2.35, height=4, bg="white")

##-----------------------------------------------------------------------------
breaks_with_gap <- c(levels(h_mem$features.plot), "", levels(h_tf$features.plot))
ggplot(full_join(h_mem, h_tf), aes(x=id, y=features.plot, fill=avg.exp.scaled))+
  geom_tile()+
  theme_minimal()+
  scale_fill_gradient2(low="blue", mid="beige", high="red", midpoint=0)+
  guides(fill=guide_colorbar(title="Avg.\nExpr.\nScal.", frame.colour="black",
                             ticks.colour="black"))+
  theme(panel.grid.major=element_blank(), axis.title=element_blank())+
  scale_y_discrete(limits = breaks_with_gap, labels = breaks_with_gap)
ggsave("../figures/FLASHseq_HeatMap_3.png", dpi=600, width=2.35, height=3.25, bg="white")

# GSEA
##-----------------------------------------------------------------------------
library(ReactomePA)
library(enrichplot)
gsea <- readRDS("../../data_processed/FLASHseq/GSEA_res.Rds")


##-----------------------------------------------------------------------------
gseaplot(gsea, geneSetID="R-HSA-389948", by="runningScore",
         title=gsea@result%>%filter(ID=="R-HSA-389948")%>%pull(Description))
ggsave("../figures/FLASHseq_GSEAlot_PD1.png", dpi=600, width=7, height=4.5)

gseaplot(gsea, geneSetID="R-HSA-201722", by="runningScore",
         title=gsea@result%>%filter(ID=="R-HSA-201722")%>%pull(Description))
ggsave("../figures/FLASHseq_GSEAlot_bCat_TCF.png", dpi=600, width=7, height=4.5)

##-----------------------------------------------------------------------------
gd <- readRDS("../../data_processed/FLASHseq/gd.Rds")
Idents(gd) <- "cluster"

##-----------------------------------------------------------------------------
pal_seur <- c("0"="#4682B4", "1"="#AF46B4", "2"="#B47846", "3"="#4BB446")
pal_clu <- c("C0"="#4682B4", "C1"="#AF46B4", "C2"="#B47846", "C3"="#4BB446")

pal_dis <- c("control"="#008B8B","UC"="#8B0000")
names_dis <- c("control"="HD", "UC"="UC (A)")



# Fig 1
##-----------------------------------------------------------------------------
# 1a
do_DimPlot(gd, label=T, group.by="seurat_clusters", legend.position="none",
           colors.use=pal_seur)
ggsave("../figures/FLASHseq_UMAP_by_Cluster.pdf", width=2.7, height=3)

##-----------------------------------------------------------------------------
# 1b
do_DimPlot(gd, group.by="disease", legend.position="left")+
  scale_color_manual(values=pal_dis ,labels=names_dis)+
  theme(legend.position=c(.2, .15))
ggsave("../figures/FLASHseq_UMAP_by_disease.pdf", width=2.7, height=3)


##-----------------------------------------------------------------------------
# 1d
Idents(gd) <- "cluster"

genes <- c("ITGAE","ITGA1","CD160","GZMA","ID3","KIT","IRF8","TCF7", "IL7R", 
           "PDCD1", "GZMK","CD5","CD27","CCR7","ZNF683","TNF","IFNG", "EOMES",
           "S1PR1","SELL","FCGR3A", "GNLY", "NKG7","GZMB", "PRF1","TBX21",
           "ZEB2","KLRG1", "ITGB2")

p <- DotPlot(gd, group.by="cluster", features=genes)+
  scale_color_gradient2(low="navy", mid="beige", high="darkred", midpoint=0)+
  guides(color=guide_colorbar(frame.colour="black", frame.linewidth=.5,
          ticks.colour="black", ticks.linewidth=.5, title="Avg.\nScal.\nExpr."),
         size=guide_legend(override.aes=list(fill="white", shape=21),title="% Expr.")
        )+
  geom_point(aes(size=pct.exp), shape=21, colour="black", stroke=.5)+
  theme(axis.title=element_blank(),axis.ticks.length.x=unit(1.5, "mm"),
        panel.grid.major=element_line(color="gray95"), legend.justification="top")+
  coord_flip()+RotatedAxis()

help("RotatedAxis")

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

p[[2]]+ggtitle(label="Tissue Residency\nSignature",
               subtitle=expression(italic("Kumar et al., 2017")))

ggsave("../figures/FLASHseq_VlnPlot_TRM.pdf", width=3.5, height=3.5)

rm(i, p)


##-----------------------------------------------------------------------------
# 1g+1h
do_DimPlot(gd, reduction="umap_scenic", colors.use=pal_clu, pt.size=1.5,
           plot.title="SCENIC Regulons")+
  theme(legend.position=c(.85, .17), legend.background=element_blank(),
        plot.title=element_text(hjust=.5))
ggsave("../figures/FLASHseq_SCENIC_UMAP_by_Cluster.pdf", width=3.5, height=3.5)

do_DimPlot(gd, reduction="umap_scenic", pt.size=1.5, group.by="disease")+
  scale_color_manual(values=pal_dis, labels=names_dis)+
  theme(legend.position=c(.87, .17), legend.background=element_blank())
ggsave("../figures/FLASHseq_SCENIC_UMAP_by_Disease.pdf", width=3.5, height=3.5)


##-----------------------------------------------------------------------------
# 1i

# Import AUC scores
rss <- importAUCfromText(
  "../../data_processed/FLASHseq/pyscenic_output_hvg_jakob_no_ENSG_nes_2.8.csv")
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

rm(rss, regulons)

# Fig 3
##-----------------------------------------------------------------------------
# Parameters
pal_trdv <- setNames(c("#8C564BFF", "#9467BDFF", "#D62728FF", "#2CA02CFF", 
                       "#FF7F0EFF", "#1F77B4FF"), levels(gd$TRDV))

pal_trgv <- setNames(c("#BCBD22FF", "#7F7F7FFF", "#E377C2FF", "#8C564BFF", "#9467BDFF",
    "#D62728FF", "#2CA02CFF", "#FF7F0EFF", "#1F77B4FF"), levels(gd$TRGV)) 

arranged <- gd@meta.data%>%count(TRGV_TRDV)%>%arrange(desc(n))
pal_gvdv <- setNames(c("#1F77B4FF", "#8C564BFF", "darkgreen", "darkred", "darkorange3"),
                     arranged$TRGV_TRDV[2:6])


##-----------------------------------------------------------------------------
# 3b
do_DimPlot(gd, group.by="TRDV",
           idents.keep=unique(gd$TRDV[!is.na(gd$TRDV)]),na.value="gray80")+
  scale_color_manual(values=pal_trdv)+
  guides(color=guide_legend(byrow=T, nrow=3, override.aes=list(size=4)))+
  theme(legend.margin=margin(t=-35), legend.justification="left")
ggsave("../figures/FLASHseq_UMAP_by_TRDV.pdf", width=3.5, height=3)

# 3d
do_DimPlot(gd, group.by="TRGV",
           idents.keep=unique(gd$TRGV[!is.na(gd$TRGV)]),na.value="gray80")+
  scale_color_manual(values=pal_trgv)+
  guides(color=guide_legend(byrow=T, nrow=3, override.aes=list(size=4)))+
  theme(legend.margin=margin(t=-25), legend.justification="left")
ggsave("../figures/FLASHseq_UMAP_by_TRGV.pdf", width=3.5, height=3)

# 3g
p <- do_DimPlot(gd, group.by="TRGV_TRDV", idents.keep=arranged$TRGV_TRDV[2:6],na.value="gray80")+
  scale_color_manual(values=pal_gvdv)+
  guides(color=guide_legend(byrow=T, nrow=3, override.aes=list(size=4)))+
  theme(legend.margin=margin(t=-25), legend.justification="left")
p[[1]]$data$TRGV_TRDV <- factor(p[[1]]$data$TRGV_TRDV, levels=arranged$TRGV_TRDV[2:6])
p

ggsave("../figures/FLASHseq_UMAP_by_TRGV-TRDV.pdf", width=3.5, height=3)
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
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.text=element_text(color="black", size=12),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"))
ggsave("../figures/FLASHseq_BarPlot_TRDV_by_Disease.pdf", width=3, height=4)

ggplot(gd@meta.data, aes(x=cluster, fill=fct_rev(TRDV)))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_trdv, na.value="gray80")+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.text=element_text(color="black", size=12),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"))
ggsave("../figures/FLASHseq_BarPlot_TRDV_by_Cluster.pdf", width=3.5, height=4)

##-----------------------------------------------------------------------------
# 3e
ggplot(gd@meta.data, aes(x=disease, fill=fct_rev(TRGV)))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_trgv, na.value="gray80")+
  scale_x_discrete(labels=names_dis)+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.text=element_text(color="black", size=12),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"))
ggsave("../figures/FLASHseq_BarPlot_TRGV_by_Disease.pdf", width=2.5, height=4)

ggplot(gd@meta.data, aes(x=cluster, fill=fct_rev(TRGV)))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_trgv, na.value="gray80")+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.text=element_text(color="black", size=12),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"))
ggsave("../figures/FLASHseq_BarPlot_TRGV_by_Cluster.pdf", width=3, height=4)


##-----------------------------------------------------------------------------
# 3h
gd@meta.data%>%
  mutate(TRGV_TRDV=factor(TRGV_TRDV, levels=arranged$TRGV_TRDV[2:6])
  )%>%
ggplot(aes(x=cluster, fill=TRGV_TRDV))+
  geom_bar(position="fill")+
  theme_classic()+
  scale_fill_manual(values=pal_gvdv, na.value="gray80")+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.text=element_text(color="black", size=12),
        legend.justification="top", legend.key.spacing.y=unit(1, "mm"))
ggsave("../figures/FLASHseq_BarPlot_TRDV-TRGV_by_Cluster.pdf", width=3.5, height=4)

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

ggsave("../figures/FLASHseq_BarPlot_TRGVs-TRDVs_by_Disease.pdf", width=4, height=3)


##-----------------------------------------------------------------------------
# 3j


##-----------------------------------------------------------------------------
# 3k










































