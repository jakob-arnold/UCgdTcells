---
title: "PBMC scRNA-seq + FACS: Code for Figures"
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
library(CATALYST)
library(ggpubr)
library(ggplotify)


# Data
##-----------------------------------------------------------------------------
vd1 <- readRDS("../../data_processed/10X/Vd1.Rds")
vd2 <- readRDS("../../data_processed/10X/Vd2.Rds")
vd1_fcs <- readRDS("../../data_processed/FACS/vd1.Rds")


# Parameters
##-----------------------------------------------------------------------------
pal_clu <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
             "#CC79A7", "#999999", "#00446B", "#E3DD84"
            )
pal_SeurClu <- setNames(pal_clu, paste0(0:9))
pal_clu <- setNames(pal_clu, paste0("C",0:9))

pal_dis <- c("HD"="#008B8B","UC"="#8B0000")
names_dis <- c("UC"="UC-A")

pal_trgv <-c("#56B4E9", "#999999", "#CC79A7", "#D55E00", "#009E73", "#F0E442",
             "#0072B2", "#E69F00", "#5BEBC5", "#444444"
             )
pal_trgv <- setNames(pal_trgv, levels(vd1$TRGV))


# Figure 8
## A-C
##-----------------------------------------------------------------------------
b <- do_DimPlot(vd2,
                label=T,
                label.size=2.12,
                group.by="seurat_clusters",
                legend.position="none",
                colors.use=pal_SeurClu,
                pt.size=0.1,
                border.size=4
                )

b[[1]][["layers"]][[3]][["geom"]][["default_aes"]][["alpha"]] <- .7


c <- do_DimPlot(vd2,
           group.by="disease",
           colors.use=pal_dis,
           pt.size=0.1,
           border.size=4,
           legend.icon.size=2.12
           )+
  scale_color_manual(values=pal_dis, labels=names_dis)+
  theme(legend.position=c(.25, 1),
        text=element_text(size=6),
        legend.background=element_blank(),
        legend.key.spacing.y=unit(-2, "mm"),
        legend.text=element_text(size=6)
        )

abc <- ggarrange(NULL, b, c,
                 ncol=3, nrow=1, widths=c(2, 1, 1),
                 font.label=list(size=9, face="bold", family="sans"),
                 labels=c("A", "B", "C")
                 )

## A-K
##-----------------------------------------------------------------------------
def <- ggarrange(NULL, NULL, NULL,
                 ncol=3, nrow=1, widths=c(1, 2, 1),
                 font.label=list(size=9, face="bold", family="sans"),
                 labels=c("D", "E", "F")
                 )

ghijk <- ggarrange(NULL, NULL, NULL, NULL, NULL,
                   ncol=5, nrow=1, widths=c(3, 1.75, 1.75, 1.75, 1.75),
                   font.label=list(size=9, face="bold", family="sans"),
                   labels=c("G", "H", "I", "J", "K")
                   )

a_k <- ggarrange(abc, def, ghijk,
                 ncol=1, nrow=3
                 )


## LM
##-----------------------------------------------------------------------------
l <- do_DimPlot(vd1,
                label=T,
                group.by="seurat_clusters",
                legend.position="none",
                colors.use=pal_SeurClu,
                pt.size=0.1,
                label.size=2.12,
                border.size=4
                )

l[[1]][["layers"]][[3]][["data"]] <- l[[1]][["layers"]][[3]][["data"]]%>%
  mutate(umap_1=case_when(seurat_clusters=="8" ~ 1.5, .default=umap_1))
l[[1]][["layers"]][[3]][["geom"]][["default_aes"]][["alpha"]] <- .7


m <- do_DimPlot(vd1,
                group.by="disease",
                colors.use=pal_dis,
                pt.size=0.1,
                legend.icon.size=2.12,
                border.size=4
                )+
  scale_color_manual(values=pal_dis, labels=names_dis)+
  theme(legend.position=c(.2, .8),
        legend.key.spacing.y=unit(-2, "mm"),
        text=element_text(size=6),
        legend.text=element_text(size=6)
        )

lm <- ggarrange(l, m,
                ncol=1, nrow=2,
                labels=c("L", "M"),
                font.label=list(size=9, face="bold", family="sans")
                )

## N
##-----------------------------------------------------------------------------
vd1$cluster <- factor(vd1$cluster, levels=c(
  "C8","C0","C1","C3","C6","C2","C4","C5", "C7"
  ))

genes <- c("GZMB", "PRF1", "TBX21", "CD74", "HLA-DRB1", "HLA-DRA", "ZNF683", "SOX4",
           "TCF7", "LEF1", "IL7R", "CD27", "CCR7", "SELL", "MKI67")

n <- DotPlot(vd1,
             group.by="cluster",
             features=genes,
             dot.scale=3
             )+
  scale_color_gradient2(low="navy", mid="beige", high="red3", midpoint=0)+
  guides(color=guide_colorbar(frame.colour="black",
                              frame.linewidth=.5,
                              ticks.colour="black",
                              ticks.linewidth=.25,
                              barwidth=.6,
                              barheight=unit(2, "cm"),
                              order=1,
                              title="Avg. Scal. Expr."),
         size=guide_legend(override.aes=list(fill="white", shape=21),
                           title="% Expr.",
                           order=2)
        )+
  geom_point(aes(size=pct.exp+12), shape=21, fill=NA, colour="black", stroke=.25)+
  theme(axis.title=element_blank(),
        axis.ticks.length.x=unit(1.5, "mm"),
        panel.grid.major=element_line(color="gray95"),
        legend.justification="top",
        text=element_text(size=6),
        axis.text=element_text(size=6),
        axis.text.y=element_text(face="italic"),
        legend.key.spacing.y=unit(-2, "mm"),
        legend.ticks.length=unit(c(-1, 0), "mm"),
        legend.text=element_text(size=6)
        )+
  coord_flip()+
  RotatedAxis()


n_1 <- vd1@meta.data%>%mutate(cluster=factor(cluster, levels=levels(n$data$id)))

n_1 <- ggplot(n_1, aes(x=cluster, fill=disease))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_dis, labels=names_dis)+
  theme_void()+
  theme(text=element_text(size=6),
        legend.title=element_blank(),
        legend.position=c(1.3, .5),
        legend.key.size=unit(2, "mm"), 
        legend.spacing.y=unit(0.5, "mm"),
        legend.margin=margin(t=-8),
        legend.text=element_text(size=6)
        )+
  guides(fill=guide_legend(ncol=1, byrow=T))+
  geom_hline(yintercept=.5, linetype="dashed", linewidth=.25)

n <- n/n_1+plot_layout(heights = c(1, 0.1))


## O
##-----------------------------------------------------------------------------
o <- plotDR(vd1_fcs,
            color_by="meta6")+
  scale_color_manual(values=pal_SeurClu)+
  theme_void()+
  geom_point(size=.25)+
  guides(color=guide_legend(nrow=3, title=NULL, override.aes=list(size=2.5)))+
  NoLegend()
  
labels <- o$data%>%
  group_by(meta6) %>%
  summarize(
    x=median(x, na.rm=TRUE),
    y=median(y, na.rm=TRUE)
  )

o <- o+
  geom_label(data=labels, aes(x=x, y=y, label=meta6),
             color="black",
             fill = alpha("white", 0.7),
             size=2.12
             )


## PQ
##-----------------------------------------------------------------------------
clusters <- vd1_fcs@metadata[["cluster_codes"]]%>%
  mutate(cluster_id=som100)

vd1_meta <- as.data.frame(vd1_fcs@colData@listData)%>%
  left_join(clusters, by="cluster_id")

p <- ggplot(vd1_meta, aes(x=disease, fill=meta6))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_SeurClu, labels=paste0("C", 1:6))+
  scale_x_discrete(labels=c("UC"="UC-A"))+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  theme_classic()+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(),
        axis.text.x=element_text(hjust=1, vjust=1, angle=45),
        axis.ticks=element_line(color="black"),
        axis.ticks.length=unit(1.5, "mm"),
        axis.text=element_text(color="black", size=6),
        text=element_text(size=6),
        legend.margin=margin(l=0),
        legend.justification="top",
        legend.key.size=unit(2.12, "mm"),
        legend.key.spacing.y=unit(1, "mm"),
        legend.box.spacing=unit(1, "mm"),
        legend.text=element_text(size=6)
        )

df <- plotExprHeatmap(vd1_fcs, scale="last", by="cluster_id", k="meta6")@matrix%>%
  as.data.frame()%>%
  rownames_to_column("cluster")%>%
  pivot_longer(cols=rownames(vd1_fcs))%>%
  mutate(cluster=factor(cluster, levels=c("2","5","4","1","3","6")),
         name=factor(name, levels=c("GZMB", "GZMK", "GZMA", "TCF-1", "CD127", "CD28"))
         )

q <- df%>%
  ggplot(aes(x=cluster, y=name, fill=value))+
  geom_tile(color="white", size=.5)+
  theme_minimal()+
  scale_fill_gradient2(low="navy", mid="beige", high="red3", midpoint=.5)+
  scale_x_discrete(labels=paste0(paste0("C", levels(df$cluster))))+
  guides(fill=guide_colorbar(frame.colour="black",
                              frame.linewidth=.5,
                              ticks.colour="black",
                              ticks.linewidth=.25,
                              barwidth=.6,
                              barheight=unit(2, "cm"),
                              order=1,
                             ))+
  labs(fill="Med.\nScal.\nExpr.")+
  theme(panel.grid.major=element_blank(), 
        panel.border=element_rect(color="black", fill=NA, size=1),
        axis.title=element_blank(),
        axis.text=element_text(color="black", size=6),
        text=element_text(size=6),
        legend.justification="top",
        legend.ticks.length=unit(c(-1, 0), "mm"),
        legend.box.spacing=unit(1, "mm"),
        legend.text=element_text(size=6)
        )

pq <- ggarrange(p, q,
          widths=c(1.2, 2),
          font.label=list(size=9, family="sans", face="bold"),
          labels=c("P", "Q"),
          align="h"
          )


## OPQ
##-----------------------------------------------------------------------------
opq <- ggarrange(o, pq,
                 ncol=1, nrow=2,
                 font.label=list(size=9, family="sans", face="bold"),
                 labels=c("O", "")
                 )


## L-Q
##-----------------------------------------------------------------------------
l_q <- ggarrange(lm , n, opq,
                ncol=3, nrow=1,
                font.label=list(size=9, family="sans", face="bold"),
                labels=c("", "N", ""),
                widths=c(4, 6, 7)
                )


## Combined
##-----------------------------------------------------------------------------
ggarrange(a_k, l_q,
          ncol=1, nrow=2,
          heights=c(1.5, 1)
          )

ggsave("../figures/Fig8.pdf",
       width=7.3,
       height=9
       )


# Fig S8
##-----------------------------------------------------------------------------
vd2@meta.data <- vd2@meta.data%>%
  mutate(diseaseID=case_when(
  diseaseID=="HD_R1_1" ~ "HD1",
  diseaseID=="HD_R1_2" ~ "HD2",
  diseaseID=="HD_R1_3" ~ "HD3",
  diseaseID=="HD_R2_1" ~ "HD4",
  diseaseID=="HD_R2_2" ~ "HD5",
  diseaseID=="UC_R1_4" ~ "UC1",
  diseaseID=="UC_R1_5" ~ "UC2",
  diseaseID=="UC_R1_6" ~ "UC3",
  diseaseID=="UC_R2_3" ~ "UC4",
  diseaseID=="UC_R2_4" ~ "UC5",
  diseaseID=="UC_R2_5" ~ "UC6",
  diseaseID=="UC_R2_6" ~ "UC7"
  ))

vd1@meta.data <- vd1@meta.data%>%
  mutate(diseaseID=case_when(
  diseaseID=="HD_R1_1" ~ "HD1",
  diseaseID=="HD_R1_2" ~ "HD2",
  diseaseID=="HD_R1_3" ~ "HD3",
  diseaseID=="HD_R2_1" ~ "HD4",
  diseaseID=="HD_R2_2" ~ "HD5",
  diseaseID=="UC_R1_4" ~ "UC1",
  diseaseID=="UC_R1_5" ~ "UC2",
  diseaseID=="UC_R1_6" ~ "UC3",
  diseaseID=="UC_R2_3" ~ "UC4",
  diseaseID=="UC_R2_4" ~ "UC5",
  diseaseID=="UC_R2_5" ~ "UC6",
  diseaseID=="UC_R2_6" ~ "UC7"
  ))



##-----------------------------------------------------------------------------
b <- readRDS("../../data_processed/10X/Vd2_QC_plots.Rds")&
  theme(text=element_text(size=6),
        axis.text=element_text(size=6),
        plot.title=element_text(size=6, face="plain")
        )
b <- ggarrange(b[[1]], b[[2]], b[[3]],
               ncol=3, nrow=1
               )

c <- readRDS("../../data_processed/10X/Vd1_QC_plots.Rds")&
  theme(text=element_text(size=6),
        axis.text=element_text(size=6),
        plot.title=element_text(size=6, face="plain")
        )

bc <- ggarrange(b, c,
                ncol=1, nrow=2,
                labels=c("B", "C"),
                font.label=list(size=9, family="sans")
                )

abc <- ggarrange(NULL, bc,
                 ncol=2, nrow=1,
                 widths=c(2.5, 1),
                 labels=c("A", ""),
                 font.label=list(size=9, family="sans")
                 )

vd2$cluster <- factor(vd2$cluster,
                      levels=c("C0","C1","C4","C5","C8","C9","C2","C3","C7","C6"))

genes <- c("GZMB", "GNLY", "PRF1", "NKG7", "TCF7", "IL7R", "GZMK", "LEF1", 
           "CD27", "CCR7", "S1PR1", "SELL", "RORC", "IL23R", "CCR6", "CXCR6", 
           "SCART1", "BLK"
           )

d <- DotPlot(vd2, group.by="cluster", features=genes)+
  scale_color_gradient2(low="blue3", mid="beige", high="red3", midpoint=0)+
  scale_radius(range=c(0.5, 4))+
  guides(color=guide_colorbar(frame.colour="black",
                              frame.linewidth=.5,
                              order=1,
                              ticks.colour="black",
                              ticks.linewidth=.5,
                              title="Avg. Scal.\nExpr.",
                              barwidth=.6,
                              barheight=unit(1.5, "cm")
                              ),
         size=guide_legend(override.aes=list(fill="white", shape=21),
                           title="% Expr.",
                           order=2)
        )+
  geom_point(aes(size=(pct.exp+6)), shape=21, fill=NA, colour="black", stroke=.25)+
  theme(axis.title=element_blank(),
        panel.grid.major=element_line(color="gray95"),
        axis.text=element_text(size=6),
        axis.text.x=element_text(face="italic"),
        text=element_text(size=6),
        legend.ticks.length = unit(c(-1, 0), 'mm'),
        legend.key.spacing.y=unit(-1, "mm"),
        legend.box.spacing=unit(1, "mm"),
        legend.title=element_text(margin=margin(b=2.5)),
        legend.spacing.y=unit(2, "mm"),
        legend.justification="top"
        )+
  RotatedAxis()

vd2$cluster <- as.character(vd2$cluster)

e <- ggplot(vd2@meta.data, aes(x=cluster, fill=diseaseID))+
  geom_bar(position="fill", color="white", size=.5)+
  scale_fill_manual(values=c(rep("#008B8B", 5),
                    rep("#8B0000", 7)))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        axis.line.x=element_blank(),
        axis.text=element_text(color="black", size=6),
        axis.ticks=element_line(color="black"),
        axis.text.x=element_text(angle=90, vjust=.5),
        legend.position="none",
        text=element_text(size=6),
        plot.margin=margin(5, 30, 5, 5)
        )

de <- ggarrange(d, e,
                ncol=2, nrow=1,
                labels=c("D", "E"),
                font.label=list(size=9, family="sans"),
                widths=c(2, 1)
                )

chord.df <- getCirclize(vd2, "strict", group.by="seurat_clusters",proportion=F)%>%
              arrange(from)

f <- as.ggplot(function() {
  par(cex=.5)
  chordDiagram(chord.df,
               annotationTrack = c("name", "grid"),
               self.link=1,
               grid.col=pal_SeurClu
               )})+
  theme(aspect.ratio=1)

chord.df <- getCirclize(vd2, "strict", group.by="diseaseID",proportion=F)%>%
              arrange(from)

g <- as.ggplot(function() {
  par(cex=.5)
  chordDiagram(chord.df,
               annotationTrack = c("name", "grid"),
               self.link=1,
               grid.col=setNames(scCustomize_Palette(length(unique(chord.df$from))), unique(chord.df$from))
               )})+
  theme(aspect.ratio=1)

h <- clonalOccupy(vd2, x.axis="cluster", label=F)&
  scale_fill_manual(values=c(
      "#fffe9e","#f79522", "#c43271", "#611062", "gray30"),
                    labels=c("Hyper (>100)", "Large (26-100)",
                    "Medium (5-25)", "Rare (2-5)","Single Clone"),
                    name="Clone Size")&
  labs(y="Number of Cells", x=NULL)&
  theme(axis.text=element_text(color="black", size=6),
        axis.ticks=element_line(color="black"),
        legend.position=c(.8, .7),
        legend.key.size=unit(2.12, "mm"),
        text=element_text(size=6),
        legend.key.spacing.y=unit(.5, "mm")
        )

fgh <- ggarrange(f, g, h,
                 ncol=3, nrow=1,
                 widths=c(1, 1, 1.5),
                 labels=c("F", "G", "H"),
                 font.label=list(size=9, family="sans")
                 )

d_h <- ggarrange(de, fgh,
                 ncol=1, nrow=2)

d_i <- ggarrange(d_h, NULL,
                 ncol=2, nrow=1,
                 widths=c(2.5, 1),
                 labels=c("", "I"),
                 font.label=list(size=9, family="sans")
                 )

chord.df <- getCirclize(vd1, "strict", group.by="seurat_clusters",proportion=F)%>%
              arrange(from)

j <- as.ggplot(function() {
  par(cex=.5)
  chordDiagram(chord.df,
               annotationTrack = c("name", "grid"),
               self.link=1,
               grid.col=pal_SeurClu
               )})+
  theme(aspect.ratio=1)

chord.df <- getCirclize(vd1, "strict", group.by="diseaseID",proportion=F)%>%
              arrange(from)

k <- as.ggplot(function() {
  par(cex=.5)
  chordDiagram(chord.df,
               annotationTrack = c("name", "grid"),
               self.link=1,
               grid.col=setNames(scCustomize_Palette(length(unique(chord.df$from))), unique(chord.df$from))
               )})+
  theme(aspect.ratio=1)

l <- clonalOccupy(vd1, x.axis="cluster", label=F)&
  scale_fill_manual(values=c(
      "#fffe9e","#f79522", "#c43271", "#611062", "gray30"),
                    labels=c("Hyper (>100)", "Large (26-100)",
                    "Medium (5-25)", "Rare (2-5)","Single Clone"),
                    name="Clone Size")&
  labs(y="Number of Cells", x=NULL)&
  theme(axis.text=element_text(color="black", size=6),
        axis.ticks=element_line(color="black"),
        legend.position=c(.8, .7),
        legend.key.size=unit(2.12, "mm"),
        text=element_text(size=6),
        legend.key.spacing.y=unit(.5, "mm")
        )

m <- ggplot(vd1@meta.data, aes(x=cluster, fill=diseaseID))+
  geom_bar(position="fill", color="white", size=.5)+
  scale_fill_manual(values=c(rep("#008B8B", 5),
                    rep("#8B0000", 7)))+
  theme_classic()+
  scale_y_continuous(expand=expansion(c(0, 0)), name="Relative Proportion")+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        axis.line.x=element_blank(),
        axis.text=element_text(color="black", size=6),
        axis.ticks=element_line(color="black"),
        axis.text.x=element_text(angle=90, vjust=.5),
        legend.position="none",
        text=element_text(size=6),
        plot.margin=margin(5, 30, 5, 5)
        )

j_m <- ggarrange(j, k, l, m,
                 ncol=4, nrow=1,
                 widths=c(2.5, 2.5, 3.75, 3.5),
                 labels=c("J", "K", "L", "M"),
                 font.label=list(size=9, family="sans")
                 )

ggarrange(abc, d_i, j_m,
          ncol=1, nrow=3,
          heights=c(1, 1.75, 1, 1)
          )

ggsave("../figures/FigS8.pdf",
       width=7.3, height=9,
       bg="white"
       )


# Session Info
##-----------------------------------------------------------------------------
sessionInfo()






































