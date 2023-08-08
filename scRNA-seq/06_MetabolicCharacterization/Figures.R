pacman::p_load(data.table, pheatmap, cowplot)

source('~/macrophages_sc/bin/compass/Compass.R')
MC='MC20'

###
#
setwd(paste0('~/macrophages_sc/results/compass/version7/', MC))
load('combined/lista_pathways.RData')
lista <- readRDS('combined/wilcoxon_results.rds')

# do not change
reactions <- fread('combined/reactions.tsv') %>% column_to_rownames('V1')
cell_metadata <- read.csv('combined/cell_metadata.csv', header = TRUE, row.names = 1)
reaction_metadata = fread('/scr/marcopretti/R/x86_64-pc-linux-gnu-library/4.0/compassR/extdata/RECON2/reaction_metadata.csv')
labeled_reactions = fread('~/macrophages_sc/results/compass/labeled_reactions.tsv')

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

## Select reactions with confidence
# More specifically, we remove reactions with confidence other than 0 or 4 (4 = most confident; 0 = unassigned confidence)
# and filter out reactions in the citric acid cycle subsystem which are outside of the mitochondria
labeled_reactions_filtered = dplyr::filter(labeled_reactions, confidence %in% c(0,4))
labeled_reactions_filtered = labeled_reactions_filtered[!is.na(EC_number),]

### Figure
col = colorRampPalette(c('blue','white','red'))
#p1 <- plots[[1]]

p1 <- plot.Compass.path(wilcoxon_results = lista[[1]], names(lista)[1])+
  geom_rect(xmin=15.5, xmax=16.5, ymin=-.3, ymax=1.6, alpha=0, size=.05, color='black')+
  theme(axis.text.y = element_blank())+
  ggtitle(label = gsub('__', ' versus ', names(lista)[1]))

# seta = data.frame(x=1, y=2) %>%
#   ggplot()+geom_curve(aes(x=x, xend=x, y=y, yend=y-2), color='black',
#                       arrow = arrow(length = unit(.05, 'npc')),
#                       show.legend = FALSE, size=1, curvature = -.5)+theme_void()
# seta

# Middle
del = dcast(dplyr::select(lista.all, -c(size, sample)), subsystem+REF ~SUBJ, value.var='med_cohens')
del = del %>% dplyr::filter(subsystem=='Glutamate metabolism') %>% dplyr::select(-subsystem) %>%
  column_to_rownames('REF')
del[is.na(del)] <- 0
p2 <- pheatmap::pheatmap(del, cluster_rows = F, cluster_cols = F, color = col(100),
                         main = 'Glutamate metabolism', angle_col = 315)

# Left
pathways_summary_heatmap = pathways_summary
pathways_summary_heatmap[is.na(pathways_summary_heatmap)] = 0
p3 <- pheatmap(pathways_summary_heatmap, scale = 'column', color = col(100), treeheight_col = 5, treeheight_row = 10,
               angle_col = 315, fontsize_row = 8, fontsize_col = 12)

pathways_summary_plot = pathways_summary %>% as.data.frame() %>% rownames_to_column('subsystem') %>% melt(value.name = 'med_cohens')
pathways_summary_size_plot = pathways_summary_size %>% as.data.frame() %>% rownames_to_column('subsystem') %>% melt(value.name = 'med_size')

png('Figures/Supp_Figure_compass.png', 2000, 1800, res=150)
plot_grid(plot_grid(plot_grid(p1, p2$gtable, nrow=2, rel_heights = c(1,.75))), 
          p3$gtable, nrow=1, rel_widths = c(.7, 1))
dev.off()


# Long heatmap not clustered
png('Figures/Supp_Figure_compass_heat_not_clustered.png', 2000, 1800, res=150)
plot_grid(plot_grid(plot_grid(p1, p2$gtable, nrow=2, rel_heights = c(1,.75))), 
          pheatmap(pathways_summary_heatmap, scale = 'column', color = col(100), treeheight_col = 5, treeheight_row = 10,
                   angle_col = 315, fontsize_row = 8, fontsize_col = 12, cluster_rows = FALSE, cluster_cols = FALSE)$gtable, nrow=1, rel_widths = c(.7, 1))

dev.off()

## END long version

########## Short version
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(Seurat)
col = colorRampPalette(c('blue','white','red'))

ind = c('Glycolysis/gluconeogenesis', # glycolysis
        'Citric acid cycle', # TCA cycle
        'Oxidative phosphorylation', # OXPHOS
        'NAD metabolism', # NADPH oxidase activity
        'Arginine and Proline Metabolism', # use of arginine
        'Transport, mitochondrial', # mitochondrial respiration
        'Fatty acid oxidation', 'Fatty acid synthesis', # fatty acid utilization
        'Glycerophospholipid metabolism','Glycosphingolipid metabolism','Sphingolipid metabolism') # production of inflammatory lipids

pathways_summary_short = pathways_summary[ind,]
pathways_summary_short = pathways_summary_short[,sort(colnames(pathways_summary_short))]
p1=pheatmap(t(pathways_summary_short), scale = 'row', color = col(100), treeheight_col = 5, treeheight_row = 10,
            angle_col = 315, fontsize_row = 8, fontsize_col = 8, cluster_rows = F, cluster_cols = F)

png('~/macrophages_sc/results/compass/Figures/pathways_summary_short.png', 800, 800, res=150)
p1
dev.off()

p1=pheatmap(t(pathways_summary_short), scale = 'row', color = col(100), treeheight_col = 5, treeheight_row = 10,
            angle_col = 315, fontsize_row = 8, fontsize_col = 8)

png('~/macrophages_sc/results/compass/Figures/pathways_summary_short_clustered.png', 800, 800, res=150)
p1
dev.off()

#
macro <- readRDS("~/macrophages_sc/results/scanpy/version7/downstream/macrophages/downstream_macrophages_ensembl_090822.RDS")
macro_rows = mapIds(EnsDb.Hsapiens.v86, rownames(macro), 'SYMBOL','GENEID')

macro_mat = GetAssayData(macro)
macro_mat = CreateSeuratObject(macro_mat)
Idents(macro_mat) = Idents(macro)
macro_mat = NormalizeData(macro_mat)

macro_mat = AverageExpression(macro_mat, features = names(ind))[[1]] %>%
  `rownames<-`(ind) %>% t()
macro_mat = macro_mat[,colSums(macro_mat) > 0]
macro_mat = macro_mat[sort(colnames(pathways_summary_short)),]

## Aggregate gene 'classes'
# ECM
macro_mat_summ <- data.frame(row.names = rownames(macro_mat),
                             FGF=Matrix::rowMeans(macro_mat[,grepl('^FGF\\d+$', colnames(macro_mat))]),
                             VEGF=Matrix::rowMeans(macro_mat[,grepl('^VEGF', colnames(macro_mat))]),
                             MMP=Matrix::rowMeans(macro_mat[,grepl('^MMP\\d+$', colnames(macro_mat))]),
                             # class 1 ppt
                             HLA_1=Matrix::rowMeans(macro_mat[,grepl('^HLA-[ABC]', colnames(macro_mat))]),
                             B2M=macro_mat[,grepl('^B2M$', colnames(macro_mat))],
                             TAP=Matrix::rowMeans(macro_mat[,grepl('^TAP[1B]', colnames(macro_mat))]),
                             # class 2 Ag ppt
                             HLA_2=Matrix::rowMeans(macro_mat[,grepl('^HLA-D[PQR]', colnames(macro_mat))]))


# Normalize by columns (z-score) and then z-score by row in the heatmap
macro_mat_summ_z <- apply(macro_mat_summ, 2, scale)
rownames(macro_mat_summ_z) = rownames(macro_mat_summ)

# Inflammatory intermediates ?
# NA: production of itaconate

# Rows not clustered
p0 = macro_mat_summ_z %>% pheatmap(cluster_rows = F, cluster_cols = T, scale = 'row', color = col(100), 
                                   angle_col = 315, treeheight_col = 20)
p1=pheatmap(t(pathways_summary_short), scale = 'row', color = col(100), treeheight_col = 5, treeheight_row = 10,
            angle_col = 315, fontsize_row = 8, fontsize_col = 8, cluster_rows = F, cluster_cols = T)

png('../../Figures/Main_compass_AgPpt_EMT_clusteredCols.png', 1500, 800, res=150)
plot_grid(p0$gtable,
          p1$gtable, align = 'hv', scale = .8, rel_widths = c(1,1.2))
dev.off()

p0 = macro_mat_summ_z %>% pheatmap(cluster_rows = T, cluster_cols = T, angle_col = 315, scale = 'row', color = col(100))
p1=pheatmap(t(pathways_summary_short), scale = 'row', color = col(100), treeheight_col = 5, treeheight_row = 10,
            angle_col = 315, fontsize_row = 8, fontsize_col = 8, cluster_rows = T, cluster_cols = T)

png('../../Figures/Main_compass_AgPpt_EMT_clusteredColsRows.png', 1500, 800, res=150)
plot_grid(p0$gtable, p1$gtable, scale = .9)
dev.off()

###########
as.data.frame(pca_res$x) %>% rownames_to_column('Cluster') %>%
  mutate(`Cell type`=gsub('(\\S+)_\\S+','\\1',Cluster)) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Cell type`))+geom_point()+
  ggrepel::geom_text_repel(aes(label=Cluster))


