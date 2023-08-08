#!/usr/bin/env R
### Integration
###-----------------------------------------------------------------------------###
DOC="Monocle3"
AUTHOR="Gabriela Rapozo"
###-----------------------------------------------------------------------------###


cds <- readRDS("data/cds.RDS")
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 20)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))

pt.matrix1 <- as.matrix(exprs(cds)[match(genes, rownames(rowData(cds))), order(pseudotime(cds))])
pt.matrix <-
  t(apply(pt.matrix1, 1, function(x) {
    smooth.spline(x, df = 3)$y
  }))
pt.matrix <- t(apply(pt.matrix, 1, function(x) {
  (x - mean(x)) / sd(x)
}))
rownames(pt.matrix) <- genes
ensembl <-
  mapIds(
    keys = row.names(pt.matrix),
    EnsDb.Hsapiens.v86,
    column = "SYMBOL",
    keytype = "GENEID",
    multiVals = "first"
  )
gene_id <- as.data.frame(ensembl)
rownames(pt.matrix) <- gene_id$ensembl
colnames(pt.matrix) <- colnames(pt.matrix1) #nao perder colnames

pseudotime <- pseudotime(cds) %>% as.data.frame() #pseudotime
pseudotime$cell <- rownames(pseudotime)

celltype <- cds@clusters$UMAP %>% as.data.frame() #celltype
celltype$cell <- rownames(celltype)

merge <-
  merge(pseudotime, celltype, by = 'cell') #pseudotime column '.', celltype columns clusters
merge <- merge[order(merge$.), ]

colnames(pt.matrix) <- merge$clusters #set colnames to clusters

ann_colors <- list(
  Celltype = c(
    "Mac_Alv-like" = "#8B87FE",
    "Mac_Angio" = "#DC4A82",
    "Mac_Hypo" = "#FFA8DE",
    "Mac_IFN" = "#6B6AF4",
    "Mac_LA" = "#BD26DE",
    "Mac_Rec" = "#BF8FFF",
    "Mac_Reg" = "#00AEEE",
    "RTM-like_MT" = "#579C99",
    'Mono_CD14_FOS+' = "#ffaf65",
    'Mono_CD14_FOS-' = "#ef3142",
    'Mono_CD16' = "#755c46",
    'Mono_IL1B' = "#d35b00",
    'MonoInter_FOS+' = "#b88b97",
    'MonoInter_FOS-' = "#b0003b",
    'cDC2_CD14' = "#2e525e",
    'cDC2_FCGR3A' = "#cfa8d5"
  )
)

ha = HeatmapAnnotation(
  Celltype = colnames(pt.matrix),
  Pseudotime = anno_barplot(merge$., gp = gpar(fill = 1:8616, col = '#912568')),
  col = list(Celltype = pal)
)

ha2 = rowAnnotation(foo = anno_mark(at = c(1:9), labels = gene_id[1:9, 1]))

hthc <- Heatmap(
  pt.matrix,
  name = "z-score",
  col = colorRamp2(seq(
    from = -2, to = 2, length = 11
  ), viridis::plasma(n = 11)),
  show_row_names = TRUE,
  top_annotation = ha,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  row_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot = 0,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  raster_by_magick = T,
  raster_quality = 2
)

hthc = plot(hthc)
