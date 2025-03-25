##-----------------------------------------------------------------------------
library(tidyverse)
library(HDCytoData)
library(CATALYST)


##-----------------------------------------------------------------------------
fcs <- read.flowSet(path="../../data_raw/FACS")
markernames(fcs)
markernames(fcs)[] <- c("GZMB", "CD127", "CD28", "TCF-1", "GZMA", "GZMK")
markernames(fcs)


##-----------------------------------------------------------------------------
uc <- c("7516", "6919", "7355", "6344", "6588", "7658", "6593", "6324", "5189")

meta <- data.frame(file_name=list.files("../../data_raw/FACS"))%>%
  separate(file_name, into=c("panel", "TRDV", "id"), sep="_|\\.fcs", remove=F)%>%
   mutate(disease=case_when(
     id %in% uc ~ "UC",
     .default="HD"
   ))%>%
  mutate(ID=paste(disease, id, sep="_"))

vd1 <- prepData(
  fcs,
  md=meta,
  md_cols=list(file = "file_name", id="ID", factors=c("disease", "TRDV")),
  transform=T, cofactor=150, by_time=F, FACS=T
)

vd1@colData@listData[["ID"]] <- vd1@colData@listData[["sample_ID"]]

summary(vd1$disease)
summary(vd1$sample_id)


# Clustering
##-----------------------------------------------------------------------------
vd1 <- runDR(vd1, "UMAP", features=rownames(vd1), cells=1000)


##-----------------------------------------------------------------------------
set.seed(1337)
vd1 <- cluster(vd1, features=rownames(vd1), seed=1337)


##-----------------------------------------------------------------------------
saveRDS(vd1, "../../data_processed/FACS/vd1.Rds")


##-----------------------------------------------------------------------------
sessionInfo()

















