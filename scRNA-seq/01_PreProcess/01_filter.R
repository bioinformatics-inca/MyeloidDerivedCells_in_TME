###!/usr/bin/env Rscript
### Integration_01
###-----------------------------------------------------------------------------###
DOC="script_tme_find_cluster_markers"
AUTHOR="MarcoPretti"
###-----------------------------------------------------------------------------###
### declare local variables
.libPaths("~/tme_ovario/lib/Rpackages/4.0/")

# load packages

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(future)
})

# Set environment
#plan(strategy = "multicore", workers = 24)
options(future.globals.maxSize= 6144*1024^2)
setwd("~/macrophages_sc/data/pre_integration_out_convert/")
set.seed(1)

####
# Read pre-built Seurat objects into a list

print(paste("Read in individual objects", Sys.time()))
files = list.files()
seurat.list <- lapply(files, readRDS)
names(seurat.list) <- unlist(lapply(strsplit(files, "\\."), "[", 1))


# Filter out cells, different thresholds depending on technology
print("Filter out low quality cells")
seurat.list.filt=list()
for (x in names(seurat.list)) {
  if (head(seurat.list[[x]]@meta.data$harm_tech,1) == "10x") {
    seurat.list.filt[[x]] <- subset(seurat.list[[x]], subset = percent.mt < 20)
  }
  if (head(seurat.list[[x]]@meta.data$harm_tech,1) == "smartseq2") {
      seurat.list.filt[[x]] <- subset(seurat.list[[x]], subset = percent.mt < 20)
      }
  if (head(seurat.list[[x]]@meta.data$harm_tech,1) == "inDrop") {
    seurat.list.filt[[x]] <- subset(seurat.list[[x]], subset = percent.mt < 20)
  }
}; rm(x)

# Count cells pre and pos filtering
sum(unlist(lapply(seurat.list, ncol)))
sum(unlist(lapply(seurat.list.filt, ncol)))

rm(seurat.list)
#saveRDS(seurat.list.filt, file = "~/macrophages_sc/results/pre_integration_scanpy/seurat.list.filt.rds")

seurat.list = seurat.list.filt
 
# metadata
print(paste("Write metadata", Sys.time()))
for (x in 1:length(seurat.list)) {
  write.csv(
    as.data.frame(seurat.list[[x]]@meta.data),
    file = paste(
      '~/macrophages_sc/data/pre_integration_csv/metadata/',
      names (seurat.list[x]),
      '_metadata.csv',
      sep = ""
    ),
    row.names = T
  )
}

# counts
print(paste("Write csv", Sys.time()))
for (x in 1:length(seurat.list)) {
  print(x)
  write.csv(
    (as.data.frame(as.matrix(seurat.list[[x]]@assays$RNA@counts))),
    file = paste(
      '~/macrophages_sc/data/pre_integration_csv/counts/',
      names (seurat.list[x]),
      '_counts.csv',
      sep = ""
    ),
    row.names = T
  )
} 

rm(x)

print(paste("Done", Sys.time()))

