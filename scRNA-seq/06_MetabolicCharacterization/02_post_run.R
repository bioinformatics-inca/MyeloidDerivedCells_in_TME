## Post compass run
## To gather the microclustered files into one

# load required libraries
library(dplyr)
library(data.table)
library(Seurat)
library(tibble)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
# If mouse data change gene_id_col_name = "HGNC.symbol" in the script to MGI.symbol


# Set working folder
setwd("~/macrophages_sc/results/scanpy/version7/")
objeto <- readRDS("~/macrophages_sc/results/scanpy/version7/downstream/macrophages/downstream_macrophages_ensembl_090822.RDS")

# # Choose on matrix to read
mat = GetAssayData(objeto, slot = "counts")
colnames(mat) <- objeto$celltype

# Change ENSG ID to HNSG
rownames(mat) <- mapIds(EnsDb.Hsapiens.v86, rownames(mat), 'SYMBOL','GENEID')
mat = mat[!is.na(rownames(mat)),]


# Run for both Microclustering parameters: 10 and 20
for (agg in c(20)) {
  
  sample.info = data.frame(cluster=rep(unique(colnames(mat)), each=agg)) %>% 
    mutate(rows=paste0(cluster, '_', 1:agg)) %>% column_to_rownames('rows')
  
  # Aggregate all clusters
  dir=paste0("~/macrophages_sc/results/compass/version7/MC",agg,"/")
  setwd(dir)
  
  micro <- list.files(pattern = 'micropools.tsv', recursive = TRUE)
  
  files <- list.files(pattern = 'reactions.tsv', recursive = TRUE)
  reac <- lapply(files, function(read){
    dt <- fread(read); dt = dt %>% column_to_rownames('V1')
    print(paste(read, 'ok!'))
    return(dt)
  })
  
  files <- list.files(pattern = 'secretions.tsv', recursive = TRUE)
  secr <- lapply(files, function(read){
    dt <- fread(read); dt = dt %>% column_to_rownames('V1') 
    print(paste(read, 'ok!'))
    return(dt)
  })
  
  files <- list.files(pattern = 'uptake.tsv', recursive = TRUE)
  uptake <- lapply(files, function(read){
    dt <- fread(read); dt = dt %>% column_to_rownames('V1') 
    print(paste(read, 'ok!'))
    return(dt)
  })
  
  files <- list.files(pattern = 'linear_gene_expression_matrix.tsv', recursive = TRUE)
  linear <- lapply(files, function(read){
    dt <- fread(read); 
    dt = dt[!duplicated(dt[['SYMBOL']]),] # remove duplicated rownames
    dt = dt %>% column_to_rownames('SYMBOL') 
    print(paste(read, 'ok!'))
    return(dt)
  })
  
  for (i in 1:length(files)) {
    # Read the microclustering cluster names
    micro_1 <- fread(micro[i]) %>% mutate(microcluster=paste0('cluster_', microcluster))
    # Rename the lists replacing 'cluster_' by their actual cluster name matching the original name
    colnames(reac[[i]]) <- micro_1$V1[match(colnames(reac[[i]]), micro_1$microcluster)]
    colnames(secr[[i]]) <- micro_1$V1[match(colnames(secr[[i]]), micro_1$microcluster)]
    colnames(uptake[[i]]) <- micro_1$V1[match(colnames(uptake[[i]]), micro_1$microcluster)]
    colnames(linear[[i]]) <- micro_1$V1[match(colnames(linear[[i]]), micro_1$microcluster)]
  }

  reac = do.call('cbind', reac)
  secr = do.call('cbind', secr)
  uptake = do.call('cbind', uptake)
  linear = do.call('cbind', linear)
  
  # Save files to a folder
  dir.create(paste0(dir, '/combined'))
  
  fwrite(reac, paste0(dir, '/combined','/reactions.tsv'), sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)
  fwrite(secr, paste0(dir, '/combined','/secretions.tsv'), sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)
  fwrite(uptake, paste0(dir, '/combined','/uptake.tsv'), sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)
  fwrite(linear, paste0(dir, '/combined','/linear_gene_expression_matrix.tsv'), sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)
}

# save the cell_metadata
lapply(micro, function(i){
  file <- fread(i)
}) %>%
  do.call('rbind',.) %>%
  mutate(cluster=gsub('(\\S+)_\\d+','\\1',V1)) %>%
  dplyr::select(-microcluster) %>% column_to_rownames('V1') %>%
  write.csv(paste0(dir,"combined/cell_metadata.csv"), row.names = T, quote = F)

