## Script to demonstrate a Compass run
## Tutorial:
# https://yoseflab.github.io/Compass/Compass-Postprocessing-Tutorial.html

# load required libraries
library(dplyr)
library(data.table)
library(Seurat)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)

# Set working folder
setwd("~/macrophages_sc/results/scanpy/version7/")
objeto <- readRDS("~/macrophages_sc/results/scanpy/version7/downstream/macrophages/downstream_macrophages_ensembl_090822.RDS")

# Choose on matrix to read
mat = GetAssayData(objeto, slot = "counts")
colnames(mat) <- objeto$celltype

# Count cells on cluster
table(colnames(mat))

# Change ENSG ID to HNSG
rownames(mat) <- mapIds(EnsDb.Hsapiens.v86, rownames(mat), 'SYMBOL','GENEID')
mat = mat[!is.na(rownames(mat)),]


# Divide per cluster to allow faster processing by Compass
col = unique(colnames(mat))
mat.list <- lapply(col, function(x){
  df = mat[,colnames(mat) %in% x]
})
names(mat.list) <- col


# Load Compass homemade functions + those from the tool
source("~/macrophages_sc/bin/compass/Compass.R")

# Run Compass for the entire dataset
# dir="macrophages_sc/results/compass/version7/all_MC100"
# agg=100
# dir.create(dir)
# date()
# prep.compass(mat = mat, dir = dir, agg = agg, name = paste0("MC_",agg))
# date()

# Run Compass on each cluster separately and for different microcluster sizes
for (agg in c(50,20,10)) {
  dir=paste0("macrophages_sc/results/compass/version7/MC",agg) # version7
  for(i in names(mat.list)){
    tmp = mat.list[[i]]
    colnames(tmp) <- paste(colnames(tmp), seq(1, ncol(tmp)), sep = '_')
    dir.create(dir)
    prep.compass(mat = tmp, dir = dir, agg = agg, name = paste0("MC_",agg,'_',i))
  }
}
