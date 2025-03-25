---
title: "PBMC FACS of Vd1 Cells: Code for Figures"
output: github_document
---



##-----------------------------------------------------------------------------
library(tidyverse)
library(CATALYST)


##-----------------------------------------------------------------------------
vd1 <- readRDS("../../data_processed/FACS/vd1.Rds")


##-----------------------------------------------------------------------------
clusters <- vd1@metadata[["cluster_codes"]]%>%
  mutate(cluster_id=som100)

vd1_meta <- as.data.frame(vd1@colData@listData)%>%
  left_join(clusters, by="cluster_id")


##-----------------------------------------------------------------------------
pal_clu <- c("#4682B4", "#6846B4", "#B44946", "#238212", "#CDAD00", "#46B488")

pal_dis <- c(HD="#008B8B", UC="#8B0000")



##-----------------------------------------------------------------------------
plotDR(vd1, color_by="meta6")+
  scale_color_manual(values=pal_clu, labels=paste0("C", 1:6))+
  theme_void()+
  geom_point(size=.5)+
  guides(color=guide_legend(nrow=3, title=NULL,override.aes=list(size=2.5)))+
  theme(legend.position=c(.55, .85))

ggsave("../figures/FACS_UMAP_byCluster.pdf", width=3, height=3)
ggsave("../figures/FACS_UMAP_byCluster.png", width=3, height=3, dpi=600)



##-----------------------------------------------------------------------------
plotDR(vd1, color_by="disease")+
  scale_color_manual(values=pal_dis, labels=c("UC"="UC (A)"))+
  theme_void()+
  geom_point(size=.5)+
  labs(color=NULL)+
  theme(legend.position=c(.5, .8))

ggsave("../figures/FACS_UMAP_Disease.pdf", width=3, height=3)
ggsave("../figures/FACS_UMAP_Disease.png", width=3, height=3, dpi=600)



##-----------------------------------------------------------------------------
ggplot(vd1_meta, aes(x=disease, fill=meta6))+
  geom_bar(position="fill")+
  scale_fill_manual(values=pal_clu, labels=paste0("C", 1:6))+
  scale_x_discrete(labels=c("UC"="UC (A)"))+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  theme_classic()+
  labs(x=NULL, y="Proportion of Cells", fill=NULL)+
  theme(axis.line.x=element_blank(), axis.text.x=element_text(hjust=1, vjust=1, angle=45),
        axis.ticks=element_line(color="black"), axis.ticks.length=unit(1.5, "mm"),
        axis.text=element_text(color="black", size=12), legend.margin=margin(l=-10),
        legend.justification="top", legend.key.spacing.y=unit(.5, "mm"))

ggsave("../figures/FACS_BarPlot_DiseaseByCluster.pdf", width=2, height=4)
ggsave("../figures/FACS_BarPlot_DiseaseByCluster.png", width=2, height=4, dpi=600)



##-----------------------------------------------------------------------------
df <- plotExprHeatmap(vd1, scale="last", by="cluster_id", k="meta6")@matrix%>%
  as.data.frame()%>%
  rownames_to_column("cluster")%>%
  pivot_longer(cols=rownames(vd1))%>%
  mutate(cluster=factor(cluster, levels=c("2","5","4","1","3","6")),
         name=factor(name, levels=c("GZMB", "GZMK", "GZMA", "TCF-1", "CD127", "CD28"))
         )


##-----------------------------------------------------------------------------
df%>%
  ggplot(aes(x=cluster, y=name, fill=value))+
  geom_tile(color="white", size=.5)+
  theme_minimal()+
  scale_fill_gradient2(low="navy", mid="beige", high="red3", midpoint=.5)+
  scale_x_discrete(labels=paste0(paste0("C", levels(df$cluster))))+
  guides(fill=guide_colorbar(frame.colour="black", frame.linewidth=.5,
                             ticks.colour="black", ticks.linewidth=.25
                             ))+
  labs(fill="Med.\nScal.\nExpr.")+
  theme(panel.grid.major=element_blank(), 
        panel.border=element_rect(color="black", fill=NA, size=1),
        axis.title=element_blank(), axis.text=element_text(color="black"),
        legend.key.width=unit(.5, "cm"), legend.key.height=unit(.45, "cm"),
        legend.justification="top")

ggsave("../figures/FACS_HeatMap_byCluster.pdf", width=4, height=3)
ggsave("../figures/FACS_HeatMap_byCluster.png", width=4, height=3, dpi=600)


##-----------------------------------------------------------------------------
ggplot(vd1_meta, aes(x=meta6, fill=sample_id))+
  geom_bar(position="fill", color="white", size=.25)+
  scale_x_discrete(labels=paste0("C", 1:6))+
  scale_y_continuous(expand=expansion(c(0, 0)))+
  scale_fill_manual(values=c(rep("#008B8B", 7), rep("#8B0000", 9)))+
  theme_classic()+
  labs(x=NULL, y="Proportion of Cells")+
  theme(axis.line.x=element_blank(), axis.text.x=element_text(hjust=1, vjust=1, angle=45),
        axis.ticks=element_line(color="black"), axis.ticks.length=unit(1.5, "mm"),
        axis.text=element_text(color="black", size=12), legend.position="none")

ggsave("../figures/FACS_BarPlot_DiseaseByDonor.pdf", width=3, height=4)
ggsave("../figures/FACS_BarPlot_DiseaseByDonor.png", width=3, height=4, dpi=600)



##-----------------------------------------------------------------------------
sessionInfo()






