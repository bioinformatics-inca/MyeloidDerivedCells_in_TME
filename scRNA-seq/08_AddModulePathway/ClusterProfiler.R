#!/usr/bin/env R
### Integration
###-----------------------------------------------------------------------------###
DOC="Patwhay Analysis"
AUTHOR="Nayara Toledo"
###-----------------------------------------------------------------------------###

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(EnsDb.Hsapiens.v86)
  library(ReactomePA)
  library(clusterProfiler)
})

set.seed(1024)

data <- readRDS("data/seurat_obj.RDS")

deg <- FindAllMarkers(data, test.use = 'MAST', only.pos = TRUE)

top100 <- deg %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
df <- top100[, 8:6]
dfsample <- split(df$Gene, df$cluster)
dfsample <- lapply(dfsample, bitr, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
genelist <- list()
for (i in 1:length(dfsample)) {
  genelist[[i]] <-  dfsample[[i]][["ENTREZID"]]
}
names(genelist) = names(dfsample)

Rclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway")

pathcDCs$data$Description = factor(pathcDCs$data$Description)