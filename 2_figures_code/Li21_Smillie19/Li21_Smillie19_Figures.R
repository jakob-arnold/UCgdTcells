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


# Palettes/Parameters
##-----------------------------------------------------------------------------
pal_dis <- c("HD"="#008B8B","UC"="#8B0000")
pal_SeuClu <- c("0"="#F0E442","1"="#56B4E9","2"="#009E73","3"="#D55E00", 
            "4"=  "#999999","5"= "#CC79A7")
pal_clu <- c("C0"="#F0E442","C1"="#56B4E9","C2"="#009E73","C3"="#D55E00", 
            "C4"=  "#999999","C5"= "#CC79A7")
names_dis <- c("UC"="UC-A")


# Enterocytes Pseudobulk
##-----------------------------------------------------------------------------
entero_pb <- AggregateExpression(entero, group.by = c("disease", "ID"),
                                 return.seurat=T)


# Fig 6
##ABH
##-----------------------------------------------------------------------------
# a
 a <- do_DimPlot(gd,
                 group.by="seurat_clusters",
                 colors.use=pal_SeuClu,
                 label=T,
                 legend.position="none",
                 label.size=2.12,
                 pt.size=.25
                 )

b <- do_DimPlot(gd,
                group.by="disease",
                pt.size=.25,
                legend.icon.size=2.12
                )+
  scale_color_manual(values=pal_dis, labels=c("UC"="UC-A"))+
  theme(legend.position=c(.8, .9),
        legend.background=element_blank(),
        legend.key.spacing.y=unit(-2, "mm"),
        legend.text=element_text(size=6),
        text=element_text(size=6)
        )

ab <- ggarrange(a, b,
                labels=c("A", "B"),
                font.label=list(size=9, face="bold", family="sans")
                )

# h
h <- do_DimPlot(entero,
                group.by="disease",
                pt.size=.1,
                border.size=4.5,
                legend.icon.size=2.12
                )+
  scale_color_manual(values=pal_dis, labels=c("UC"="UC-A"))+
  labs(title="Enterocytes")+
  theme(legend.position=c(.8, .95),
        legend.background=element_blank(),
        legend.key.spacing.y=unit(-2, "mm"),
        text=element_text(size=6),
        legend.text=element_text(size=6),
        plot.title=element_text(hjust=.5, size=6)
        )

abh <- ggarrange(ab, h,
                 ncol=1, nrow=2, font.label=list(size=9, face="bold", family="sans"),
                 labels=c("", "H"),
                 heights=c(7, 9)
                 )


## C
##-----------------------------------------------------------------------------
gd <- SetIdent(gd, value="cluster")
gd$cluster <- factor(gd$cluster, levels=c("C0", "C3", "C2", "C1", "C4", "C5"))

genes <- c( "ITGAE","ITGA1","CD160","GZMA","ID3","KIT","IRF8","GZMK","CD5",
            "CD27","ZNF683","FCGR3A", "GNLY", "NKG7","GZMB","PRF1","TBX21",
            "ZEB2",  "KLRG1", "ITGB2")

c1 <- DotPlot(gd,
              group.by="cluster",
              features=genes,
              dot.scale=3.25
              )+
  scale_color_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  guides(color=guide_colorbar(frame.colour="black",
                              frame.linewidth=.5,
                              ticks.colour="black",
                              ticks.linewidth=.25,
                              barwidth=.6,
                              barheight=unit(2, "cm"),
                              order=1,
                              title="Avg. Scal. Expr."
                              ),
        size=guide_legend(override.aes=list(fill="white", shape=21),
                          title="% Expr.",
                          order=2
                          )
       )+
  geom_point(aes(size=pct.exp+20),
             shape=21,
             fill=NA,
             colour="black",
             stroke=0.25
             )+
  theme(axis.title=element_blank(),
        axis.ticks.length.x=unit(1.5, "mm"),
        legend.justification="top",
        legend.key.spacing.y=unit(-2, "mm"),
        panel.grid.major=element_line(color="gray95"),
        text=element_text(size=6),
        axis.text=element_text(size=6),
        axis.text.y=element_text(face="italic"),
        legend.ticks.length=unit(c(-1, 0), "mm"),
        legend.text=element_text(size=6)
        )+
  coord_flip()+
  RotatedAxis()

# C2
genes <- c("SELL", "S1PR1", "EOMES", "IFNG", "TNF", "CCR7", "PDCD1", "IL7R", "TCF7")

c2_1 <- DotPlot(gd,
              group.by="cluster",
              features=genes,
              dot.scale=3.75
              )+
  scale_color_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  guides(color=guide_colorbar(frame.colour="black",
                              frame.linewidth=.5,
                              ticks.colour="black",
                              ticks.linewidth=.25,
                              barheight=unit(2, "cm"),
                              barwidth=.6,
                              title="Avg. Scal. Expr.",
                              order=1
                              ),
        size=guide_legend(override.aes=list(fill="white", shape=21),
                          title="% Expr.",
                          order=2
                          )
       )+
  geom_point(aes(size=pct.exp+4),
             shape=21,
             fill=NA,
             colour="black",
             stroke=0.25)+
  theme(axis.title=element_blank(),
        axis.ticks.length.x=unit(1.5, "mm"),
        legend.justification=c("top"),
        legend.key.spacing.y=unit(-2, "mm"),
        legend.key.width=unit(4, "mm"),
        panel.grid.major=element_line(color="gray95"),
        text=element_text(size=6),
        axis.text=element_text(size=6),,
        axis.text.y=element_text(face="italic"),
        legend.ticks.length=unit(c(-1, 0), "mm"),
        legend.text=element_text(size=6)
        )+
  coord_flip()+
  RotatedAxis()


c2_2 <- gd@meta.data%>%mutate(cluster=factor(cluster, levels=levels(c2_1$data$id)))%>%
  ggplot(aes(x=cluster, fill=disease))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_dis, labels=c("UC"="UC-A"))+
  theme_void()+
  theme(text=element_text(size=6),,
        legend.title=element_blank(),
        legend.position=c(-.35, .6),
        legend.key.size=unit(2, "mm"),
        legend.spacing.y=unit(0.5, "mm"),
        legend.text=element_text(size=6),
        plot.background=element_blank()
        )+
  guides(fill=guide_legend(ncol=1, byrow=T))+
  geom_hline(yintercept=.5, linetype="dashed", linewidth=.35)

c <- ggarrange(c1, c2_1, c2_2,
               ncol=1, nrow=3, heights=c(3.5, 2, .75),
               font.label=list(size=9, family="sans", face="bold"),
               labels=c("C", "", ""),
               align="hv"
               )

## D-G
##-----------------------------------------------------------------------------
defg <- VlnPlot_scCustom(gd,
                         sort=T,
                         colors_use=pal_clu,
                         group.by="cluster",
                         pt.size=0,
                         num_columns=1,
                         add.noise=F,
                         features=c("Vg4Vd1_1", "Vg9Vd2_1", "cyto1", "li_stem1")
                         )&
  scale_y_continuous(expand=expansion(c(0.05, 0.05)))&
  geom_boxplot(width=.3, fill="white", outlier.shape=NA, coef=0, color="black", size=.25)&
  theme(plot.subtitle=element_text(face="italic", size=6),
        plot.title=element_text(hjust=0, size=6, face="plain"),
        text=element_text(size=6),
        axis.text=element_text(size=6)
        )&
  labs(x=NULL, y="Module Score")
  
defg[[1]] <- defg[[1]]+ggtitle("TRGV4-TRDV1 Signature")
defg[[2]] <- defg[[2]]+ggtitle("TRGV9-TRDV2 Signature")
defg[[3]] <- defg[[3]]+ggtitle("Cytotoxicity and Cytokine\nProduction Signature")
defg[[4]] <- defg[[4]]+labs(title="Stem-Like Signature")


for (i in seq_along(defg)){
  defg[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE
}

defg <- ggarrange(defg[[1]], defg[[2]], defg[[3]], defg[[4]],
                  ncol=1, nrow=4, font.label=list(size=9, face="bold", family="sans"),
                  labels=c("D", "E", "F", "G"),
                  align="v"
                  )

## A-H
##-----------------------------------------------------------------------------
abcdefgh <- ggarrange(abh, c, defg,
                      ncol=3, nrow=1, widths=c(4, 3.2, 2.8)
                      )


## I
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

genes_up <- c("BTN2A1", "BTN3A1", "BTN3A3", "SLC37A3", "ICAM1", "IRF1",
              "IRF8", "JAK1", "JAK2", "STAT1", "STAT2", "PRKAA1", "PRKAA2")
genes_down <- c("BTNL3", "BTNL8", "HNF4A", "CDX2", "ABCA1")

I <- ggplot()+
  geom_point(UC_markers%>%filter(signif == "ns"),
             mapping=aes(x=geneLogSums, y=avg_log2FC, color=signif),
             size=.5
             )+
  geom_point(UC_markers%>%filter(signif != "ns"),
             mapping=aes(x=geneLogSums, y=avg_log2FC, color=signif),
             size=.1
             )+
  theme_classic()+
  geom_hline(yintercept=0, linetype="dashed", color="red", size=.5)+
  labs(x=expression("log"[2]*"(Mean Expression)"),
       y=expression("log"[2]*"(Fold Change)"),
       color=element_blank())+
  scale_color_manual(values = c("gray", "#008B8B", "#8B0000"))+
  guides(color=guide_legend(override.aes=list(size = 2), reverse=T))+
  theme(axis.text=element_text(color="black", size=6),
        text=element_text(size=6),
        legend.position = c(0.85, 0.85),
        legend.key.spacing.y=unit(-2, "mm"),
        legend.text=element_text(size=6)
        )+
  
  ggrepel::geom_label_repel(data=UC_markers%>%filter(gene %in% genes_up), 
                           aes(label = gene, x=geneLogSums, y=avg_log2FC, fontface="italic"),
                           min.segment.length=0, nudge_y=1.5, nudge_x=1,
                           size=1.9, max.overlaps=20)+
  ggrepel::geom_label_repel(data=UC_markers%>%filter(gene %in% genes_down), 
                           aes(label = gene, x=geneLogSums, y=avg_log2FC, fontface="italic"),
                           min.segment.length=0, nudge_y=-2, nudge_x=-1.5,
                           size=1.9)


## J-O
##-----------------------------------------------------------------------------
jklmo <- VlnPlot_scCustom(entero_pb,
                          group.by="disease",
                          pt.size=0,
                          colors_use=pal_dis,
                          num_columns=3,
                          adjust=1,
                          add.noise=F,
                          features=c("BTNL3", "HNF4A", "BTN3A1", "BTNL8", "CDX2", "BTN3A3")
                          )&
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.15)))&
  scale_x_discrete(labels=c("UC"="UC-A"))&
  geom_boxplot(fill="white", width=.5, outlier.shape=NA, coef=0, color="black", size=.25)&
  geom_jitter(width=.05, height=0, size=1, fill="white", shape=21)&
  labs(y="Expression")&
  theme(axis.title.x=element_blank(),
        plot.background=element_blank(),
        plot.title=element_text(face="italic", size=6),
        text=element_text(size=6),
        axis.text=element_text(size=6)
        )&
  stat_compare_means(method="t.test", label.x=1.5, bracket.size=0, size=2.82, label="p.signif")

for (i in seq_along(jklmo)){
  jklmo[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE
}

jklmo <- ggarrange(jklmo[[1]], jklmo[[2]], jklmo[[3]], jklmo[[4]], jklmo[[5]], jklmo[[6]],
                   ncol=3, nrow=2, font.label=list(size=9, family="sans", face="bold"),
                   labels=c("J", "L", "N", "K", "M", "O"),
                   align="v"
                   )


## I-O
##-----------------------------------------------------------------------------
ijklmo <- ggarrange(I, jklmo,
                    ncol=2, nrow=1, font.label=list(size=9, face="bold", family="sans"),
                    labels=c("I", ""),
                    widths=c(2, 3)
                    )


## PQ
##-----------------------------------------------------------------------------
genes <- c("SLC37A3", "ICAM1")

pq <- VlnPlot_scCustom(entero_pb,
                       group.by="disease",
                       pt.size=0,
                       colors_use=pal_dis,
                       num_columns=2,
                       adjust=1,
                       add.noise=F,
                       features=genes
                       )&
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))&
  scale_x_discrete(labels=c("UC"="UC-A"))&
  geom_boxplot(fill="white",width=.5, outlier.shape=NA, coef=0, color="black", size=.25)&
  geom_jitter(width=.05, height=0, size=1, fill="white", shape=21)&
  labs(y="Expression")&
  theme(axis.title.x=element_blank(),
        plot.background=element_blank(),
        plot.title=element_text(face="italic", size=6),
        text=element_text(size=6),
        axis.text=element_text(size=6)
        )

for (i in seq_along(pq)){
  pq[[i]][["layers"]][[1]][["stat_params"]][["trim"]] <- FALSE
}

test <-   function(y.pos=NULL, x.pos=NULL){
  stat_compare_means(bracket.size=0, method = "t.test", label.y=y.pos, label.x=x.pos, size=2.12,
  aes(label = ifelse(..p.. < 0.0001, "p<0.0001", paste0("p=", ..p.format..)))
                    )
                        }
pq[[1]]<-pq[[1]]+test(x.pos=0.8)
pq[[2]]<-pq[[2]]+test(x.pos=0.8)

pq <- ggarrange(pq[[1]], pq[[2]],
                font.label=list(size=9, family="sans", face="bold"),
                labels=c("P", "Q")
                )

## RS
##-----------------------------------------------------------------------------
gsea <- readRDS("../../data_processed/Li21_Smillie19_Integrated/GSEA.Rds")

# TCA
r <- gseaplot(gsea,
              geneSetID="R-HSA-71403",
              by="runningScore",
              title=gsea@result%>%filter(ID=="R-HSA-71403")%>%pull(Description)
         )+
  annotate(geom="text",
           x=7000,
           y=.15,
           label=paste0("p = ", gsea@result%>%filter(ID=="R-HSA-71403")%>%pull(p.adjust)%>%signif(2)),
           size=2.12
           )+
  theme_bw()+
  scale_x_continuous(expand=expansion(mult=c(0.05, .1)))+
  scale_y_continuous(expand=expansion(mult=c(0.05, .1)))+
  labs(y="Running\nEnrichment Score")+
  theme(plot.title=element_text(size=6),
        text=element_text(size=6),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black")
        )

r[["layers"]][[1]][["geom"]][["default_aes"]][["linewidth"]] <- .25
r[["layers"]][[2]][["aes_params"]][["size"]] <- .5

s <- gseaplot(gsea,
              geneSetID="R-HSA-77289",
              by="runningScore",
              title=gsea@result%>%filter(ID=="R-HSA-77289")%>%pull(Description)
         )+
  annotate(geom="text",
           x=7000,
           y=.15,
           label=paste0("p = ", gsea@result%>%filter(ID=="R-HSA-77289")%>%pull(p.adjust)%>%signif(2)),
           size=2.12
           )+
  theme_bw()+
  scale_x_continuous(expand=expansion(mult=c(0.05, .1)))+
  scale_y_continuous(expand=expansion(mult=c(0.05, .1)))+
  labs(y="Running\nEnrichment Score")+
  theme(plot.title=element_text(size=6),
        text=element_text(size=6),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black")
        )

s[["layers"]][[1]][["geom"]][["default_aes"]][["linewidth"]] <- .25
s[["layers"]][[2]][["aes_params"]][["size"]] <- .5

rs <- ggarrange(r, s,
                font.label=list(size=9, face="bold", family="sans"),
                labels=c("R", "S")
                )


## P-S
##-----------------------------------------------------------------------------
pqrs <- ggarrange(pq, rs,
                  widths=c(2, 3)
                  )


## Combined
##-----------------------------------------------------------------------------
ggarrange(abcdefgh, ijklmo, pqrs,
          ncol=1, nrow=3, heights=c(9, 5, 2.5)
          )

ggsave("../figures/Fig6.pdf",
       width=7.3,
       height=9
       )


# Fig S6
##-----------------------------------------------------------------------------
pb <- readRDS("../../data_processed/Smillie19/PB.Rds")

a <- do_DimPlot(gd,
                group.by="cohort",
                pt.size=.25,
                legend.icon.size=2.12
                )+
  theme(legend.position=c(.8, .9),
        legend.background=element_blank(),
        legend.key.spacing.y=unit(-2, "mm"),
        legend.text=element_text(size=6),
        text=element_text(size=6),
        plot.margin=margin(5, 5, 10, 10)
        )

b <-gd@meta.data%>%
  ggplot(aes(x=cluster, fill=diseaseID))+
  geom_bar(position="fill", color="white",size=.2)+
  scale_fill_manual(values=c(rep("#008B8B", 16), rep("#8B0000", 21)))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(),
        text=element_text(size=6),
        legend.position="none",
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text=element_text(color="black", size=6),
        axis.ticks.length=unit(2, "mm"),
        axis.ticks=element_line(color="black"),
        axis.text.x=element_text(vjust=3),
        plot.margin=margin(5, 50, 5, 5)
        )

c <- do_DimPlot(entero,
                group.by="cohort",
                pt.size=.1,
                border.size=4.5,
                legend.icon.size=2.12
                )+
  theme(legend.position=c(.8, 1),
        legend.background=element_blank(),
        legend.key.spacing.y=unit(-2, "mm"),
        legend.text=element_text(size=6),
        text=element_text(size=6),
        plot.margin=margin(5, 5, 10, 10)
        )

d <- entero@meta.data%>%
  unite("disease_cohort", disease, cohort, sep="_")%>%
  ggplot(aes(x=disease_cohort, fill=diseaseID))+
  geom_bar(position="fill", color="white",size=.2)+
  scale_fill_manual(values=c(rep("#008B8B", 16), rep("#8B0000", 22)))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(),
        text=element_text(size=6),
        legend.position="none",
        axis.line.x=element_blank(), 
        axis.text=element_text(color="black", size=6),
        axis.ticks.length=unit(2, "mm"),
        axis.ticks=element_line(color="black"),
        plot.margin=margin(5, 50, 5, 5)
        )+
  RotatedAxis()

p <- VlnPlot_scCustom(
  pb,
  group.by="compartment",
  add.noise=F,
  split.by="Health",
  features=c("BTNL3", "BTNL8", "BTN3A1", "BTN3A3"),
  pt.size=0
  )&
  geom_boxplot(
    width=.7,
    position=position_dodge(width=0.9),
    outlier.shape=NA)&
  geom_point(
    position=position_dodge(width=0.9),
    shape=21,
    size=1,
    show.legend=F
    )&
  scale_fill_manual(values= c("Healthy"="#008B8B","Inflamed"="#8B0000"),
                    c("Healthy"="HD", "Inflamed"="UC-A"))&
  scale_y_continuous(expand=expansion(mult=c(.05, .15)))&
  labs(y="Expression")&
  theme(
    axis.title.x=element_blank(),
    text=element_text(size=6),
    axis.text=element_text(size=6),
    plot.title=element_text(size=6, face="italic")
  )&
  stat_compare_means(
    bracket.size=1,
    size=2.12,
    method="t.test",
    # label.x=1.27,
    aes(label=ifelse(
      ..p.. < 0.001,
      "p<0.001",
      paste0("p=", ..p.format..))),
    vjust=-1
    )

for(i in seq_along(p)){
  p[[i]][["layers"]][[1]] <- NULL
}

q <- VlnPlot_scCustom(entero,
                 group.by="disease",
                 features=c("meva1", "chol1"),
                 add.noise=F,
                 pt.size=0,
                 colors_use=pal_dis
                 )&
  geom_boxplot(fill="white", color="black", coef=0, outliers=F, size=.25, width=.2)&
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.1)))&
  labs(y="Module Score", x=NULL)&
  theme(
    text=element_text(size=6),
    plot.title=element_text(size=6, face="plain"),
    axis.text=element_text(size=6))&
    stat_compare_means(
      bracket.size=0,
      method="t.test",
      size=2.12,
      vjust=-1,
      aes(label=ifelse(..p.. < 0.001,"p<0.001", paste0("p=", ..p.format..)))
      )

q[[1]] <- q[[1]]+ggtitle("Mevalonate")
q[[2]] <- q[[2]]+ggtitle("Cholesterol")

abc <- ggarrange(a, b, c,
                 ncol=3, nrow=1,
                 labels=c("A", "B", "C"),
                 font.label=list(size=9, family="sans")
                 )

def <- ggarrange(d, p[[1]], p[[2]],
                 ncol=3, nrow=1,
                 labels=c("D", "E", "F"),
                 font.label=list(size=9, family="sans"),
                 align="h"
                 )

gj <- ggarrange(p[[3]], p[[4]], q[[1]], q[[2]],
                ncol=4, nrow=1,
                widths=c(1, 1, .5, .5),
                labels=c("G", "H", "I", "J"),
                font.label=list(size=9, family="sans"),
                align="hv"
                )


ggarrange(abc, def, gj,
          ncol=1, nrow=3
          )

ggsave("../figures/FigS6.pdf",
       width=7.3,
       height=7,
       bg="white"
       )


# Session Info
##-----------------------------------------------------------------------------
sessionInfo()












