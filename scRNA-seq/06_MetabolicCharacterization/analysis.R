## Script to demonstrate a Compass run
## Tutorial:
# https://yoseflab.github.io/Compass/Compass-Postprocessing-Tutorial.html

# load required libraries
library(dplyr)
library(data.table)
library(Seurat)

# Set working folder
setwd("~/macrophages_sc/results/scanpy/version7/")
objeto <- readRDS('~/macrophages_sc/results/scanpy/version7/second_level_anno/third_level_anno/downstream_macrophages.rds')

# Choose on matrix to read
colnames(mat) <- objeto$celltype

mat[1:5,1:5]
dim(mat)

# Count cells on cluster
table(colnames(mat))

# Divide per cluster to allow raw processing by Compass
col = unique(colnames(mat))
mat.list <- lapply(col, function(x){
  df = mat[,colnames(mat) %in% x]
})
names(mat.list) <- col

# Load Compass homemade functions
source("~/macrophages_sc/bin/compass/Compass.R")

# Run Compass
dir="/data04/projects04/MarianaBoroni/macrophages_sc/results/compass/version7/MC10"
#setwd(dir)
agg=10
for(i in names(mat.list)){
  tmp = mat.list[[i]]
  colnames(tmp) <- paste(colnames(tmp), seq(1, ncol(tmp)), sep = '_')
  dir.create(dir)
  prep.compass(mat = tmp, dir = dir, agg = agg, name = paste0("MC_",agg,'_',i))
}

dir="/data04/projects04/MarianaBoroni/macrophages_sc/results/compass/version7/MC20"
#setwd(dir)
agg=20
for(i in names(mat.list)){
  tmp = mat.list[[i]]
  colnames(tmp) <- paste(colnames(tmp), seq(1, ncol(tmp)), sep = '_')
  dir.create(dir)
  prep.compass(mat = tmp, dir = dir, agg = agg, name = paste0("MC_",agg,'_',i))
}


sample.info = data.frame(cluster=rep(unique(colnames(mat)), each=agg)) %>% 
  mutate(rows=paste0(cluster, '_', 1:agg)) #%>% column_to_rownames('rows')

# Aggregate all clusters
setwd("~/macrophages_sc/results/compass/MC20/")

files <- list.files(pattern = 'reactions.tsv', recursive = TRUE)
reac <- lapply(files, function(read){
  dt <- fread(read); dt = dt %>% column_to_rownames('V1')
  print(paste(read, 'ok!'))
})
reac = do.call('cbind', reac)

files <- list.files(pattern = 'secretions.tsv', recursive = TRUE)
secr <- lapply(files, function(read){
  dt <- fread(read); dt = dt %>% column_to_rownames('V1') 
  print(paste(read, 'ok!'))
  return(dt)
})
secr = do.call('cbind', secr)

files <- list.files(pattern = 'uptake.tsv', recursive = TRUE)
uptake <- lapply(files, function(read){
  dt <- fread(read); dt = dt %>% column_to_rownames('V1') 
  print(paste(read, 'ok!'))
  return(dt)
})
uptake = do.call('cbind', uptake)

files <- list.files(pattern = 'linear_gene_expression_matrix.tsv', recursive = TRUE)
linear <- lapply(files, function(read){
  dt <- fread(read); dt = dt %>% column_to_rownames('SYMBOL') 
  print(paste(read, 'ok!'))
  return(dt)
})
linear = do.call('cbind', linear)

# Save files to a folder
dir="~/macrophages_sc/results/compass/MC20/"
dir.create(paste0(dir, '/combined'))

fwrite(reac, paste0(dir, '/combined','/reactions.tsv'), sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)
fwrite(secr, paste0(dir, '/combined','/secretions.tsv'), sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)
fwrite(uptake, paste0(dir, '/combined','/uptake.tsv'), sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)
fwrite(linear, paste0(dir, '/combined','/linear_gene_expression_matrix.tsv'), sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)


### >>>> TO OTHER SCRIPT
## Analyze Compass result
sample.info %>%
  write.csv(paste0(dir,"/cell_metadata.csv"), row.names = F, quote = F)

compass_settings <- CompassSettings$new(
  user_data_directory = paste0(dir,'/combined'),
  cell_id_col_name = 'cluster' ,
  gene_id_col_name = "HGNC.symbol"
)

compass_settings <- pos.compass(dir = paste0(dir, "/compass_MC_10_Macro_FOLR2_22_Sep_17_12/"), # Macrophages
                                sample.info = sample.info %>% dplyr::filter(grepl('FOLR', cluster)), cell_id_col_name = "cluster")

compass_data <- CompassData$new(compass_settings)
compass_analyzer <- CompassAnalyzer$new(compass_settings)

## LAUNCH 'analysis_shark.R'
rm(mat, objeto, reac, secr, uptake)