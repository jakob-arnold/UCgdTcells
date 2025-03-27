---
title: "Integrated Datasets (Li 2021, Smillie 2019): Code for Figures"
output: github_document
---



##-----------------------------------------------------------------------------
library(Seurat)
library(SCpubr)
library(scCustomize)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(enrichplot)

# Data
##-----------------------------------------------------------------------------
entero <- readRDS("../../data_processed/Li21_Smillie19_Integrated/enterocytes.Rds")
gd <- readRDS("../../data_processed/Li21_Smillie19_Integrated/gdTcells.Rds")


# Cluster Marker
##-----------------------------------------------------------------------------
gd%>%
  FindAllMarkers(only.pos=T)%>%
  filter(p_val_adj<0.05)%>%
  write.csv("Li21Sm19_Marker_AllClusters.csv")



# Palettes/Parameters
##-----------------------------------------------------------------------------
pal_dis <- c("HD"="#008B8B","UC"="#8B0000")
pal_SeuClu <- c("0"="#4682B4","1"= "#B44946","2"= "#238212","3"= "#CDAD00", 
            "4"=  "#6846B4","5"= "gray60")
pal_clu <- c("C0"="#4682B4","C1"= "#B44946","C2"= "#238212","C3"= "#CDAD00", 
            "C4"=  "#6846B4","C5"= "gray60")


# Enterocytes Pseudobulk
##-----------------------------------------------------------------------------
entero_pb <- AggregateExpression(entero, group.by = c("disease", "ID"),
                                 return.seurat=T)

# gd T Cells

##-----------------------------------------------------------------------------
do_DimPlot(gd, group.by="seurat_clusters", colors.use=pal_SeuClu, label=T,
           legend.position="none")

ggsave("../figures/LiSm_gd_UMAP_byCluster.pdf", width=3, height=3.5)
ggsave("../figures/LiSm_gd_UMAP_byCluster.png", dpi=600, width=3, height=3.5)

##-----------------------------------------------------------------------------
do_DimPlot(gd, group.by="disease")+
  scale_color_manual(values=pal_dis, labels=c("UC"="UC (A)"))+
  theme(legend.position=c(.8, .9))

ggsave("../figures/LiSm_gd_UMAP_byDisease.pdf", width=3, height=3.5)
ggsave("../figures/LiSm_gd_UMAP_byDisease.png", dpi=600, width=3, height=3.5)




##-----------------------------------------------------------------------------
gd <- SetIdent(gd, value="cluster")
gd$cluster <- factor(gd$cluster, levels=c("C0", "C3", "C2", "C1", "C4", "C5"))

genes <- c( "ITGAE","ITGA1","CD160","GZMA","ID3","KIT","IRF8","GZMK","CD5",
            "CD27","ZNF683","FCGR3A", "GNLY", "NKG7","GZMB","PRF1","TBX21",
            "ZEB2",  "KLRG1", "ITGB2")

p <- DotPlot(gd, group.by="cluster", features=genes, scale=T)+
  scale_color_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  guides(color=guide_colorbar(frame.colour="black", frame.linewidth=.5, order=1,
                           ticks.colour="black", ticks.linewidth=.5,
                           title="Avg.\nScal.\nExpr."),
        size=guide_legend(override.aes=list(fill="white", shape=21),
                             title="% Expr.", order=2)
       )+
  geom_point(aes(size=pct.exp+6), shape=21,fill=NA, colour="black", stroke=0.25)+
  theme(axis.title=element_blank(),
        axis.ticks.length.x=unit(1.5, "mm"), legend.justification="top",
        panel.grid.major=element_line(color="gray95"))+
  coord_flip()+RotatedAxis()


p_1 <- gd@meta.data%>%mutate(cluster=factor(cluster, levels=levels(p$data$id)))%>%
ggplot(aes(x=cluster, fill=disease))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_dis, labels=c("UC"="UC (A)"))+
  theme_void()+
  theme(text=element_text(size=18),,
        legend.title=element_blank(),
        legend.justification=c("left", "top"),
        legend.key.size=unit(4, "mm"),
        legend.spacing.y=unit(0.5, "mm"),
        plot.background=element_blank())+
  guides(fill=guide_legend(ncol=1, byrow=T))+
  geom_hline(yintercept=.5, linetype="dashed", linewidth=.75)

p/p_1+plot_layout(heights = c(1, 0.1))
ggsave("../figures/LiSm_gd_DotPlot_1.pdf", width=3.5, height=5)
ggsave("../figures/LiSm_gd_DotPlot_1.png", dpi=600, width=3.5, height=5)

##-----------------------------------------------------------------------------
genes <- c("SELL", "S1PR1", "EOMES", "IFNG", "TNF", "CCR7", "PDCD1", "IL7R", "TCF7")

p <- DotPlot(gd, group.by="cluster", features=genes, scale=T)+
  scale_color_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  guides(color=guide_colorbar(frame.colour="black", frame.linewidth=.5, order=1,
                           ticks.colour="black", ticks.linewidth=.5,
                           title="Avg.\nScal.\nExpr."),
        size=guide_legend(override.aes=list(fill="white", shape=21),
                             title="% Expr.", order=2)
       )+
  geom_point(aes(size=pct.exp+4.25), shape=21,fill=NA, colour="black", stroke=0.25)+
  theme(axis.title=element_blank(),
        axis.ticks.length.x=unit(1.5, "mm"), legend.justification="top",
        panel.grid.major=element_line(color="gray95"))+
  coord_flip()+RotatedAxis()


p_1 <- gd@meta.data%>%mutate(cluster=factor(cluster, levels=levels(p$data$id)))%>%
ggplot(aes(x=cluster, fill=disease))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_dis, labels=c("UC"="UC (A)"))+
  theme_void()+
  theme(text=element_text(size=14),,
        legend.title=element_blank(),
        legend.position=c(-.35, .6),
        legend.key.size=unit(4, "mm"),
        legend.spacing.y=unit(0.5, "mm"),
        plot.background=element_blank())+
  guides(fill=guide_legend(ncol=1, byrow=T))+
  geom_hline(yintercept=.5, linetype="dashed", linewidth=.75)

p/p_1+plot_layout(heights = c(1, 0.1))
ggsave("../figures/LiSm_gd_DotPlot_2.pdf", width=3.25, height=3.5)
ggsave("../figures/LiSm_gd_DotPlot_2.png", dpi=600, width=3.25, height=3.5)

gd$cluster <- as.character(gd$cluster)
rm(p, p_1, genes)



##-----------------------------------------------------------------------------
p <- VlnPlot_scCustom(gd, sort=T, colors_use=pal_clu, group.by="cluster",
                      pt.size=0, num_columns=1, add.noise=F,
                      features=c("Vg4Vd1_1", "Vg9Vd2_1", "cyto1", "li_stem1"))&
  scale_y_continuous(expand=expansion(c(0.05, 0.05)))&
  geom_boxplot(width=.3, fill="white", outlier.shape=NA, coef=0, color="black")&
  theme(plot.subtitle=element_text(face="italic"),
        plot.title=element_text(hjust=0, size=15))&
  labs(x=NULL)
  
p[[1]] <- p[[1]]+ggtitle("TRGV4-TRDV1 Signature")
p[[2]] <- p[[2]]+ggtitle("TRGV9-TRDV2 Signature")
p[[3]] <- p[[3]]+ggtitle("Cytotoxicity and Cytokine\nProduction Signature")
p[[4]] <- p[[4]]+labs(title="Stem-Like Signature",
                      subtitle="Li et al., 2024")


for (i in seq_along(p)){
  p[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE
}
p

ggsave("../figures/LiSm_gdT_VlnPlots.pdf", width=3.25, height=9)
ggsave("../figures/LiSm_gdT_VlnPlots.png", dpi=600, width=3.25, height=10)
rm(p, i)

##-----------------------------------------------------------------------------
scores <- grep(colnames(gd@meta.data), patter="stem1$|effector1$", value=T)

p <- VlnPlot_scCustom(gd, features=scores, group.by = "cluster",
                 pt.size = 0, colors_use = pal_clu, sort=T,
                 add.noise=F, num_columns = 2)&
  theme(axis.title = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5))&
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))&
  geom_boxplot(width=.4, fill="white", outlier.shape=NA, coef=0)

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

ggsave("../figures/LiSm_gdT_Stem_Effector_Signatures.pdf", width=6, height=9)
ggsave("../figures/LiSm_gdT_Stem_Effector_Signatures.png", dpi=600, width=6, height=9)

rm(i, scores)

##-----------------------------------------------------------------------------
do_DimPlot(gd, group.by="cohort")+
  theme(legend.position=c(.75, 1), legend.background=element_blank())

ggsave("../figures/LiSm_gd_UMAP_byCohort.pdf", width=3, height=3.5)
ggsave("../figures/LiSm_gd_UMAP_byCohort.png", dpi=600, width=3, height=3.5)

##-----------------------------------------------------------------------------
gd@meta.data%>%
  ggplot(aes(x=cluster, fill=diseaseID))+
  geom_bar(position="fill", color="white",size=.2)+
  scale_fill_manual(values=c(rep("#008B8B", 16), rep("#8B0000", 21)))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(), legend.position="none",
        axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text=element_text(color="black"), axis.ticks.length=unit(2, "mm"),
        axis.ticks=element_line(color="black"), axis.text.x=element_text(vjust=3),
        )

ggsave("../figures/LiSm_gd_BarPlot_ClusterByDonor.pdf", width=2.5, height=3)
ggsave("../figures/LiSm_gd_BarPlot_ClusterByDonor.png", width=2.5, height=3, dpi=600)



# Enterocytes
##-----------------------------------------------------------------------------
do_DimPlot(entero, group.by="disease", pt.size=.25, border.size=3)+
  scale_color_manual(values=pal_dis, labels=c("UC"="UC (A)"))+
  theme(legend.position=c(.8, .95))

ggsave("../figures/LiSm_entero_UMAP_byDisease.pdf", width=3, height=3)
ggsave("../figures/LiSm_entero_UMAP_byDisease.png", dpi=600, width=3, height=3)




## Diff. Genes
##-----------------------------------------------------------------------------
Idents(entero) <- "disease"

UC_markers <- FindMarkers(entero, ident.1="UC")%>%rownames_to_column("gene")

geneLogSums <- log2(Matrix::rowMeans(GetAssayData(entero, "RNA", "counts")))
geneLogSums <- as.data.frame(geneLogSums)%>%rownames_to_column("gene")

UC_markers <- inner_join(UC_markers, geneLogSums, by="gene")%>%
  mutate(signif=case_when(
           p_val_adj < 0.05 & avg_log2FC > 0 ~ "up",
           p_val_adj < 0.05 & avg_log2FC < 0 ~ "down",
           .default = "ns"
         ),
         signif=factor(signif, levels = c("ns", "down", "up"))
         )


##-----------------------------------------------------------------------------
genes_up <- c("BTN2A1", "BTN3A1", "BTN3A3", "SLC37A3", "ICAM1", "IRF1",
              "IRF8", "JAK1", "JAK2", "STAT1", "STAT2", "PRKAA1", "PRKAA2")
genes_down <- c("BTNL3", "BTNL8", "HNF4A", "CDX2", "ABCA1")

ggplot()+
  geom_point(UC_markers%>%filter(signif == "ns"),
             mapping=aes(x=geneLogSums, y=avg_log2FC, color=signif), size=.5)+
    geom_point(UC_markers%>%filter(signif != "ns"),
             mapping=aes(x=geneLogSums, y=avg_log2FC, color=signif), size=.1)+
  theme_classic()+
  geom_hline(yintercept=0, linetype="dashed", color="red", size=1)+
  labs(x=expression("log"[2]*"(Mean Expression)"),
       y=expression("log"[2]*"(Fold Change)"),
       color=element_blank())+
  scale_color_manual(values = c("gray", "#008B8B", "#8B0000"))+
  guides(color = guide_legend(override.aes=list(size = 2), reverse=T))+
  theme(axis.text = element_text(color="black"), text=element_text(size=16),
        legend.position = c(0.85, 0.85))+
  
  ggrepel::geom_label_repel(data=UC_markers%>%filter(gene %in% genes_up), 
                           aes(label = gene, x=geneLogSums, y=avg_log2FC),
                           min.segment.length=0, nudge_y=1.5, nudge_x=0,
                           size = 3, max.overlaps=15)+
  ggrepel::geom_label_repel(data=UC_markers%>%filter(gene %in% genes_down), 
                           aes(label = gene, x=geneLogSums, y=avg_log2FC),
                           min.segment.length=0, nudge_y=-1, nudge_x=1.5,
                           size=3)

ggsave("../figures/LiSm_entero_MA_Plot.pdf", width=5, height=4)
ggsave("../figures/LiSm_entero_MA_Plot.png", dpi=600, width=5, height=4)

rm(UC_markers, geneLogSums, genes_down, genes_up)


##-----------------------------------------------------------------------------
p <- VlnPlot_scCustom(entero_pb, group.by = "disease", pt.size=0,
            colors_use=pal_dis, num_columns=3, adjust=1, add.noise=F,
            features=c("BTNL3", "HNF4A", "BTN3A1", "BTNL8", "CDX2", "BTN3A3"))&
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))&
  scale_x_discrete(labels=c("UC"="UC (A)"))&
  geom_boxplot(fill="white",width=.5, outlier.shape=NA, coef=0, color="black")&
  geom_jitter(width=.05, height=0, size=1.5, fill="white", shape=21)&
  theme(axis.title=element_blank(),plot.background=element_blank(),
        plot.title=element_text(face="plain"))




for (i in seq_along(p)){
  p[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE
}

test <-   function(pos=NULL){
  stat_compare_means(bracket.size=0, method = "t.test", label.y=pos,
  aes(label = ifelse(..p.. < 0.0001, "p<0.0001", paste0("p=", ..p.format..)))
                    )
                        }
p[[1]]<-p[[1]]+test(1.6)
p[[2]]<-p[[2]]+test(1.7)
p[[3]]<-p[[3]]+test(0.4)
p[[4]]<-p[[4]]+test(1.72)
p[[5]]<-p[[5]]+test(2)
p[[6]]<-p[[6]]+test(0.6)

p
ggsave("../figures/LiSm_Entero_VlnPlots.pdf", width=6, height=5)
ggsave("../figures/LiSm_Entero_VlnPlots.png", dpi=600, width=6, height=5, bg="white")

rm(i, p, test)



##-----------------------------------------------------------------------------
genes <- c("SLC37A3", "ICAM1")

p <- VlnPlot_scCustom(entero_pb, group.by = "disease", pt.size=0,
            colors_use=pal_dis, num_columns=2, adjust=1, add.noise=F,
            features=genes)&
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))&
  scale_x_discrete(labels=c("UC"="UC (A)"))&
  geom_boxplot(fill="white",width=.5, outlier.shape=NA, coef=0, color="black")&
  geom_jitter(width=.05, height=0, size=1.5, fill="white", shape=21)&
  theme(axis.title=element_blank(),plot.background=element_blank(),
        plot.title=element_text(face="plain"))




for (i in seq_along(p)){
  p[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE
}

test <-   function(y.pos=NULL, x.pos=NULL){
  stat_compare_means(bracket.size=0, method = "t.test", label.y=y.pos, label.x=x.pos,
  aes(label = ifelse(..p.. < 0.0001, "p<0.0001", paste0("p=", ..p.format..)))
                    )
                        }
p[[1]]<-p[[1]]+test(x.pos=0.8)
p[[2]]<-p[[2]]+test(x.pos=0.8)

p
ggsave("../figures/LiSm_Entero_VlnPlots_2.pdf", width=4, height=2.5)
ggsave("../figures/LiSm_Entero_VlnPlots_2.png", dpi=600, width=4, height=2.5, bg="white")

rm(p, i, test, genes)


## GSEA
##-----------------------------------------------------------------------------
gsea <- readRDS("../../data_processed/Li21_Smillie19_Integrated/GSEA.Rds")

# TCA
gseaplot(gsea, geneSetID="R-HSA-71403", by="runningScore",
         title=gsea@result%>%filter(ID=="R-HSA-71403")%>%pull(Description))+
  annotate(geom="text", x=7500, y=.15, label=paste0("p = ",
    gsea@result%>%filter(ID=="R-HSA-71403")%>%pull(p.adjust)%>%signif(2)),
    size=5)+
  scale_x_continuous(expand=expansion(mult=c(0.05, .1)))+
  theme(plot.title=element_text(face="bold"))
  
ggsave("../figures/LiSm_entero_GSEA_TCA.pdf", width=6, height=4)
ggsave("../figures/LiSm_entero_GSEA_TCA.png", dpi=600,  width=6, height=4)

# Mito Fatty
gseaplot(gsea, geneSetID="R-HSA-77289", by="runningScore",
         title=gsea@result%>%filter(ID=="R-HSA-77289")%>%pull(Description))+
  annotate(geom="text", x=7500, y=.15, label=paste0("p = ",
    gsea@result%>%filter(ID=="R-HSA-77289")%>%pull(p.adjust)%>%signif(2)),
    size=5)+
  scale_x_continuous(expand=expansion(mult=c(0.05, .1)))+
  theme(plot.title=element_text(face="bold"))

ggsave("../figures/LiSm_entero_GSEA_BetaOx.pdf", width=6, height=4)
ggsave("../figures/LiSm_entero_GSEA_BetaOx.png", dpi=600,  width=6, height=4)

rm(gsea)


##-----------------------------------------------------------------------------
FeaturePlot_scCustom(entero, features=c("BTNL3", "BTNL8", "BTN2A1", "BTN3A1", "BTN3A3"), 
                     na_color="#2A0134", pt.size=.1,
                     colors_use=c("#6900A8FF", "#B02991FF", "#E16462FF","#FCA636FF", "#F0F921FF"))&
  NoLegend()&NoAxes()

ggsave("../figures/LiSm_entero_FeaturePlot_BTNs.pdf", width=5, height=6.5)
ggsave("../figures/LiSm_entero_FeaturePlot_BTNs.png", width=5, height=6.5, dpi=600)


##-----------------------------------------------------------------------------
do_DimPlot(entero, group.by="cohort", pt.size=.25)+
  theme(legend.position=c(.75, .98), legend.background=element_blank())

ggsave("../figures/LiSm_entero_UMAP_byCohort.png", dpi=600, width=3, height=3.5)
ggsave("../figures/LiSm_entero_UMAP_byCohort.pdf", width=3, height=3.5)


##-----------------------------------------------------------------------------
entero@meta.data%>%
  unite("disease_cohort", disease, cohort, sep="_")%>%
  ggplot(aes(x=disease_cohort, fill=diseaseID))+
  geom_bar(position="fill", color="white",size=.2)+
  scale_fill_manual(values=c(rep("#008B8B", 16), rep("#8B0000", 22)))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(), legend.position="none",
        axis.line.x=element_blank(), 
        axis.text=element_text(color="black"), axis.ticks.length=unit(2, "mm"),
        axis.ticks=element_line(color="black"),
        )+
  RotatedAxis()

ggsave("../figures/LiSm_entero_BarPlot_DiseaseCohortByDonor.pdf", width=2, height=3)
ggsave("../figures/LiSm_entero_BarPlot_DiseaseCohortByDonor.png", width=2, height=3, dpi=600)


##-----------------------------------------------------------------------------
sessionInfo()












