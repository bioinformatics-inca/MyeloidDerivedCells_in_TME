#!/usr/bin/env R
### Pseudotime
###-----------------------------------------------------------------------------###
DOC="Monocle3"
AUTHOR="Gabriela Rapozo"
###-----------------------------------------------------------------------------###

suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(SeuratWrappers)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  library(EnsDb.Hsapiens.v86)
})

data <- readRDS("data/mono_macro_moDC.RDS")
Idents(data) <- factor(data$celltype, levels = levels)

cds <- as.cell_data_set(data)

reacreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition
list_cluster <- data@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings ----
cds@int_colData@listData$reducedDims$UMAP <-
  data@reductions$UMAP@cell.embeddings
cds <- learn_graph(cds, learn_graph_control=list(ncenter=300))
cds <-
  order_cells(cds,
              reduction_method = 'UMAP',
              root_cells = colnames(cds[, clusters(cds) == c("Mono_CD14_FOS-")]))

saveRDS("data/cds.RDS")