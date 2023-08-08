#!/usr/bin/env R
### Integration
###-----------------------------------------------------------------------------###
DOC="AddModuleScore"
AUTHOR="GiovannaMaklouf e GabrielaRapozo"
###-----------------------------------------------------------------------------###

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(ggpubr)
  library(UCell)
})

data <- readRDS("data/seurat_obj.RDS")

create_violin_plot <- function(data, features, score_name, title) {
  data <- UCell::AddModuleScore_UCell(data, features = list(features), name = score_name)
  order <- reorder(data$celltype, -data@meta.data[[score_name]])
  data$celltype <- factor(data$celltype, levels = levels(order))
  p <- ggviolin(data@meta.data, x = "celltype", y = score, fill = "celltype",
                palette = macro_pal,
                add = "boxplot", add.params = list(fill = "white")) + 
    stat_compare_means(label = "p.signif", method = "wilcox", ref.group = ".all.", hide.ns = FALSE, show.legend = FALSE) +
    stat_compare_means(method = "kruskal.test", label.y = max(data@meta.data[[score]]) + 0.1, label.x = 1.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = mean(data@meta.data[[score]]), linetype = 2) +
    NoLegend() +
    ggtitle(title) +
    ylab(paste0(title)) +
    xlab("") + ylab("")
  
  return(p)
}

data <- create_violin_plot(data, features = score, score_name = "score", title = "Score")
