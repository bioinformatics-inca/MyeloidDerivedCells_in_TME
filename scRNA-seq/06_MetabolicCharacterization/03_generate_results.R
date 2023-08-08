library(data.table)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

agg='MC50' # MC50, MC20
dir_=paste0('~/macrophages_sc/results/compass/version7/',agg)
setwd(dir_)

# do not change
reactions <- fread('combined/reactions.tsv') %>% column_to_rownames('V1')
cell_metadata <- read.csv('combined/cell_metadata.csv', header = TRUE, row.names = 1)
reaction_metadata = fread('/scr/marcopretti/R/x86_64-pc-linux-gnu-library/4.0/compassR/extdata/RECON2/reaction_metadata.csv')
labeled_reactions = fread('~/macrophages_sc/results/compass/labeled_reactions.tsv')

# run once
# dplyr::select(compass_data$reaction_metadata,
#               c('reaction_no_direction','reaction_name','subsystem','EC_number','confidence')) %>%
#   write.table('~/macrophages_sc/results/compass/labeled_reactions.tsv', sep='\t', row.names = FALSE, quote = FALSE)

# Converts the raw penalties outputs of compass into scores per reactions
# where higher numbers indicate more activity
get_reaction_consistencies <- function(compass_reaction_penalties, min_range=1e-3){
  df = -log(compass_reaction_penalties + 1)
  df = df[,apply(df, 2, max) - apply(df, 2, min) >= min_range]
  df = df - min(apply(df, 1, min))
  return(df)
}

reaction_consistencies = get_reaction_consistencies(reactions)

# Exclude reactions with no variance
reaction_consistencies = reaction_consistencies[apply(reaction_consistencies, 1, var) > 0,]

# Function to perform wilcoxon test between groups
wilcox_test_compass <- function(reaction_consistencies, group_1, group_2){
  require(effsize)
  t1 = reaction_consistencies %>% dplyr::select(contains(group_1))
  t2 = reaction_consistencies %>% dplyr::select(contains(group_2))
  
  wilcox_tests <- lapply(rownames(t1), function(x){
    wilcox <- wilcox.test(as.numeric(t1[x,]), as.numeric(t2[x,]))
    cohen_d <- cohen.d(as.numeric(t1[x,]), as.numeric(t2[x,]))
    adjusted_pval <- p.adjust(wilcox$p.value, method = 'fdr', n=nrow(t1))
    wilcox = c(row=x, wilcox_stat=as.numeric(wilcox$statistic), wilcox_pval=wilcox$p.value, 
               cohens_d=cohen_d$estimate, adjusted_pval=adjusted_pval, metadata_r_id=x)
  })
  return(do.call('rbind', wilcox_tests) %>% as.data.frame() %>% column_to_rownames('row'))
}


## Comparar par a par cada cluster e construir uma matriz
# paired combinations
combinations <- as.data.frame(combn(unique(cell_metadata$clus), 2))

lista <- lapply(combinations, function(p){
  print(p)
  group_1 = p[1]; group_2 = p[2]
  wilcox_test_compass(reaction_consistencies, group_1 = group_1, group_2 = group_2)
})
names(lista) <- lapply(combinations, paste, collapse ='__') %>% unlist() %>% as.character()

saveRDS(lista, paste0("~/macrophages_sc/results/compass/version7/",agg,"/combined/wilcoxon_results.rds"))
lista <- readRDS(paste0("~/macrophages_sc/results/compass/version7/",agg,"/combined/wilcoxon_results.rds"))


## Select reactions with confidence
# More specifically, we remove reactions with confidence other than 0 or 4 (4 = most confident; 0 = unassigned confidence)
# and filter out reactions in the citric acid cycle subsystem which are outside of the mitochondria
labeled_reactions_filtered = dplyr::filter(labeled_reactions, confidence %in% c(0,4))
labeled_reactions_filtered = labeled_reactions_filtered[!is.na(EC_number),]

## Save paired comparisons
lista[[1]] %>% head
names(lista)

source('~/macrophages_sc/bin/compass/Compass.R')
fig_dir=paste0(dir_, '/Figures/Pairwise_comparisons')
dir.create(fig_dir, recursive = TRUE)
setwd(fig_dir)
for (comparison in names(lista)) {
  png(filename = paste0(comparison, '.png'), width = 1200, height = 1100, res = 150)
  print(plot.Compass.path(wilcoxon_results = lista[[comparison]], comparison))
  dev.off() 
}

plots <- lapply(names(lista), function(comparison){
  plot.Compass.path(wilcoxon_results = lista[[comparison]], comparison)
})

pdf('Pairwise_plots.pdf', onefile = TRUE, width = 12, height = 10)
plots
dev.off()

# Filter significant p values
lista.sub = lapply(lista, function(sub){
  sub = sub %>% 
    dplyr::filter(adjusted_pval<.1) %>%
    mutate(reaction_no_direction=gsub('(\\S+)_\\S+','\\1',metadata_r_id)) %>%
    inner_join(labeled_reactions, by='reaction_no_direction')
})


# Create summary with median cohens_d per comparison
lista.summ <- lapply(lista.sub, function(df){
  df %>% group_by(subsystem) %>% 
    summarise(med_cohens=median(as.numeric(cohens_d), na.rm=TRUE), size=n())
})

lista.all = rbindlist(lista.summ)
lista.all$sample = rep(names(lista.summ), unlist(lapply(lista.summ, nrow)))
lista.all[,c('REF','SUBJ')] = tstrsplit(lista.all$sample, '__')

pathways_summary <- sapply(unique(c(lista.all$REF, lista.all$SUBJ)), function(cluster){
  sapply(unique(lista.all$subsystem), function(sub){
    dplyr::select(lista.all, -c(sample)) %>% # reconsiderar o size depois
      dplyr::filter(REF %in% cluster | SUBJ %in% cluster) %>%
      dplyr::filter(subsystem==sub) %>% pull(med_cohens) %>% mean(na.rm=TRUE)
  })
})

pathways_summary_size <- sapply(unique(c(lista.all$REF, lista.all$SUBJ)), function(cluster){
  sapply(unique(lista.all$subsystem), function(sub){
    dplyr::select(lista.all, -c(sample)) %>% # reconsiderar o size depois
      dplyr::filter(REF %in% cluster | SUBJ %in% cluster) %>%
      dplyr::filter(subsystem==sub) %>% pull(size) %>% sd(na.rm=TRUE) # pull(size) %>% mean() 
  })
})

save(lista.all, pathways_summary, pathways_summary_size, 
     file = paste0(dir_,'/combined/lista_pathways.RData'))
