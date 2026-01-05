---
title: "Thomas et al. 2024: Code for Figures"
output: github_document
---



##-----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scCustomize)
library(ggpubr)
library(patchwork)


# Data
##-----------------------------------------------------------------------------
entero <- readRDS("../../data_processed/Thomas_24/thomas24_entero_pb.Rds")
gd <- readRDS("../../data_processed/Thomas_24/thomas24_gd_pb.Rds")


# Figure 7
##-----------------------------------------------------------------------------
# Blank Canvas
abc <- ggarrange(NULL, NULL, NULL,
                 ncol=3, nrow=1,
                 labels=c("A", "B", "C"),
                 font.label=list(size=9, face="bold", family="sans")
                 )

defgh <- ggarrange(NULL, NULL, NULL, NULL, NULL,
                  ncol=5, nrow=1,
                  labels=c("D", "E", "F", "G", "H"),
                  font.label=list(size=9, face="bold", family="sans")
                 )

abcdefgh <- ggarrange(abc, defgh,
                     ncol=1, nrow=2
                     )


##-----------------------------------------------------------------------------
scores <- c(
  "IEL1",
  "Memory-Like1",
  "Cytotoxic1"
  )

gd_plots <- VlnPlot_scCustom(
  gd,
  features=scores,
  group.by="Treatment",
  pt.size=0,
  add.noise=F
  )&
  geom_boxplot(width=.5, outlier.shape=NA, color="black")&
  geom_point(shape=21, fill="white", size=1)&
  labs(x=NULL, y="Module Score")&
  theme(text=element_text(size=6),
        axis.text=element_text(size=6),
        plot.title=element_text(face="plain", size=6)
        )
  
for (i in seq_along(gd_plots)){
  gd_plots[[i]][["data"]] <- gd_plots[[i]][["data"]]%>%
    rownames_to_column("id")%>%
    separate(id, into=c("rem", "tp", "patient"), sep="_")%>%
    mutate(rem=case_when(rem=="Remission" ~ "R", .default="NR"),
           rem=factor(rem, levels=c("R", "NR")),
           ident=factor(ident, levels=c("Pre", "Post"))
           )
  
  gd_plots[[i]][["layers"]][[1]] <- NULL
}

gd_plots[[1]] <- gd_plots[[1]]+ggtitle("IEL Score")
gd_plots[[2]] <- gd_plots[[2]]+ggtitle("Memory-Like Score")
gd_plots[[3]] <- gd_plots[[3]]+ggtitle("Cytotoxic Score")

gd_plots <- gd_plots&
  geom_line(aes(group=patient), size=.25)&
  facet_wrap(~rem, scales="free_x")&
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.2)))&
  theme(strip.background=element_rect(fill="white", color="black", size=1),
        strip.text=element_text(size=6)
        )&
  stat_compare_means(
    method="t.test",
    paired=T,
    label.x=1.5,
    size=2.82,
    label="p.signif"
  )


##-----------------------------------------------------------------------------
genes <- c(
  "BTNL3",
  "BTNL8",
  "CDX2",
  "HNF4A",
  "BTN2A1",
  "BTN3A1",
  "BTN3A3"
  )

entero_plots <- VlnPlot_scCustom(
  entero,
  features=genes,
  group.by="Treatment",
  pt.size=0,
  add.noise=F
  )&
  geom_boxplot(width=.5, outlier.shape=NA, color="black")&
  geom_point(shape=21, fill="white", size=1)&
  labs(x=NULL, y="Expression")&
  theme(text=element_text(size=6),
        axis.text=element_text(size=6),
        plot.title=element_text(face="italic", size=6)
        )
  
for (i in seq_along(entero_plots)){
  entero_plots[[i]][["data"]] <- entero_plots[[i]][["data"]]%>%
    rownames_to_column("id")%>%
    separate(id, into=c("rem", "tp", "patient"), sep="_")%>%
    mutate(rem=case_when(rem=="Remission" ~ "R", .default="NR"),
           rem=factor(rem, levels=c("R", "NR")),
           ident=factor(ident, levels=c("Pre", "Post"))
           )
  
  entero_plots[[i]][["layers"]][[1]] <- NULL
  }

entero_plots <- entero_plots&
  geom_line(aes(group=patient), size=.25)&
  facet_wrap(~rem, scales="free_x")&
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.2)))&
  theme(strip.background=element_rect(fill="white", color="black", size=1),
        strip.text=element_text(size=6)
        )&
  stat_compare_means(
    method="t.test",
    paired=T,
    label.x=1.5,
    size=2.82,
    label="p.signif"
  )


##-----------------------------------------------------------------------------
ijklm <- ggarrange(gd_plots[[1]], gd_plots[[2]], gd_plots[[3]], entero_plots[[6]], entero_plots[[7]],
                   labels=c("I", "J", "K", "L", "M"),
                   font.label=list(size=9, face="bold", family="sans"),
                   align="hv",
                   ncol=5, nrow=1
                    )


##-----------------------------------------------------------------------------
ggarrange(abcdefgh, ijklm,
          ncol=1, nrow=2,
          heights=c(5.8, 2.1)
          )

ggsave("../figures/Fig7.pdf",
       width=7.25,
       height=6.162
       )


# Session Info
##-----------------------------------------------------------------------------
sessionInfo()












