###
## Script to run compass
## More info: https://yoseflab.github.io/Compass/install.html
# https://yoseflab.github.io/Compass/Compass-Postprocessing-Tutorial.html

prep.compass <- function(mat, dir, agg=NULL, name=NULL, filter_low_exp=FALSE){
  
  ## check parameters and data

  ## Add error to mandate gene as rownames
  
  ## Signature from CPM RNASeq matrix
  require(DropletUtils, quietly = T)
  require(Matrix, quietly = T)
  require(tibble, quietly = T)
  
  if(filter_low_exp==TRUE){
    ## Transform to cpm and keeping genes expressed in at least 5% of samples
    exc = 0.05*ncol(mat)
    mat = mat[rowSums(mat > 0) > exc,]  
    mat <- edgeR::cpm(mat)
    tsv = mat
    mat = as(mat, "sparseMatrix")
  }else{
    mat <- edgeR::cpm(mat)
    tsv = mat
    mat = as(mat, "sparseMatrix")
  }
  

  # create output directory
  subdir = paste0("/compass_", name) #, '_', gsub("\\w+ (\\w+) (\\d+) (\\d+):(\\d+):\\d+ (\\S+)", "\\2_\\1_\\3_\\4", date()))    

  final.dir = paste0(dir,subdir)
  dir.create(path = final.dir, showWarnings = F)
  setwd(final.dir)
  
  print("Matrix dimension after filtering:")
  print(dim(mat))
  
  # create config file
  parameters="--num-threads 54 --num-processes 50 --calc-metabolites"
  
  if(!is.null(agg)){
    parameters = paste(parameters, '--microcluster-size', agg, '--lambda 0', '--species homo_sapiens') # see manual
  }
  
  command=paste(parameters, paste('--output-dir', final.dir))
  input=paste0(final.dir,"/linear_gene_expression_matrix")
  
  # Save matrix
  if(ncol(mat) >100 | !is.null(agg)){
    write10xCounts(x = mat, path = input, overwrite = F, type = "sparse")
    system(paste("awk '{print $1}'", paste0(input,"/genes.tsv"), '> tmp && mv tmp', paste0(input,"/genes.tsv")))

    rows <- rownames(tsv)
    cols <- colnames(tsv)
    tsv = data.table::data.table(tsv)
    tsv$SYMBOL = rows
    tsv = tsv %>% relocate(SYMBOL)
    tsv = data.table::setnames(x = tsv, new = c('SYMBOL',cols))

    write.table(tsv, file = paste0(input,'.tsv'), quote = F, sep = "\t", col.names = T, row.names = F) # needed for postprocessing

    input=paste(paste0(input,"/matrix.mtx"), paste0(input,"/genes.tsv"), paste0(input,"/barcodes.tsv"))
    compass.data = paste("compass --data-mtx", input)
  }else{
    compass.data = paste("compass --data", input)
    write.table(mat, file = input, quote = F, sep = "\t", col.names = T, row.names = F)
  }
  
  command = paste(compass.data, command)
  command = paste0("command='",command,"'")
  write.table(command, file = paste0(final.dir, "/config.txt"), quote = F, col.names = F, row.names = F)
  print("To launch compass run:")
  print(paste("sbatch Compass.slurm", paste0(final.dir, "/config.txt")))

  outdir=paste0("outdir='", final.dir, "'")
  write.table(outdir, file = paste0(final.dir, "/config.txt"), quote = F, col.names = F, row.names = F, append = T)
}


## Run script at
## Load conda environment
# . /tools/Anaconda2/etc/profile.d/conda.sh
# conda activate compass
library(compassR)
library(tidyverse)

pos.compass <- function(dir, sample.info, cell_id_col_name){
  require(compassR)
  # prepare cell_metadata.csv
  # receive df with sample and 'grouṕ' / 'condition'
  sample.info %>%
    write.csv(paste0(dir,"/cell_metadata.csv"), row.names = F, quote = F)
  
  compass_settings <- CompassSettings$new(
    user_data_directory = dir,
    cell_id_col_name = cell_id_col_name,
    gene_id_col_name = "HGNC.symbol"
  )
  return(compass_settings)
}

## PLOT: function to compare pathways between two samples
plot.Compass.path <- function(wilcoxon_results, title=NULL, return_data=FALSE){
  # Create reaction ID no direction
  wilcoxon_results = wilcoxon_results %>%
    mutate(reaction_no_direction=gsub('(\\S+)_\\S+','\\1',metadata_r_id))
  
  cohens_d_by_subsystem <- wilcoxon_results %>%
    #left_join(dplyr::select(compass_data$reaction_partitions, "reaction_id", "reaction_no_direction"), by = "reaction_id") %>%
    inner_join(labeled_reactions_filtered, by = "reaction_no_direction") %>%
    # Keep only "confident reactions", as defined in our paper.
    dplyr::filter(!is.na(EC_number)) %>%
    dplyr::filter(confidence == "0" | confidence == "4") %>%
    # Keep only "interesting subsystems", as defined in our paper.
    dplyr::filter(!(subsystem == "Miscellaneous" | subsystem == "Unassigned")) %>%
    dplyr::filter(!(startsWith(subsystem, "Transport") | startsWith(subsystem, "Exchange"))) %>%
    # Keep only subsystems of non-negligible size.
    group_by(subsystem) %>%  dplyr::filter(n() > 5) %>%  ungroup() %>%
    # Order subsystems in a manner that will lend itself to a visually aesthetic plot.
    mutate(cohens_d=as.numeric(cohens_d),
           subsystem_priority = factor(subsystem))#%>%
 #            fct_reorder2(cohens_d, adjusted_pval, .fun = function(cohens_d, adjusted_p_value) {
  #             abs(median(cohens_d[adjusted_pval < 0.1]))}, .desc = FALSE))
  
  if(return_data==TRUE){
    return(cohens_d_by_subsystem)
  }else{
    ggplot(cohens_d_by_subsystem,  aes(x = subsystem_priority,
                                       y = cohens_d,
                                       color = if_else(cohens_d > 0, "up_regulated", "down_regulated"),
                                       alpha = if_else(adjusted_pval < 0.1, "significant", "insignificant")
    )) +
      ggtitle(label = title, subtitle = "Up- and Down-Regulated Reactions Cross Pathway Boundaries") +
      xlab("") + ylab("Cohen's d") +
      scale_color_manual(
        values = c(up_regulated = "#ca0020", down_regulated = "#0571b0"),
        guide = FALSE
      ) +
      scale_alpha_manual(
        name = "",
        values = c(significant = 1, insignificant = 0.25),
        labels = c(significant = "BH-adjusted p-value < 0.1", insignificant = "Not significant")
      ) +
      coord_flip() +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_bw() +
      theme(legend.position = "bottom", legend.direction = "horizontal")
  }
  
}

## PLOT: function to perform VolcanoPlot

plot.Compass.volcano <- function(wilcoxon_results, facets=c("Glycolysis", "TCA cycle", "Fatty acid oxidation", "Amino acid metabolism")){
  require(compassR, quietly = TRUE)
  require(ggrepel, quietly = TRUE)
  require(tidyverse, quietly = TRUE)
  
  compass_scores_by_cell_type <- wilcoxon_results %>%
    left_join(dplyr::select(compass_data$reaction_partitions, "reaction_id", "reaction_no_direction"),by = "reaction_id") %>%
    left_join(compass_data$reaction_metadata, by = "reaction_no_direction") %>%
    # Keep only "confident reactions", as defined in our paper.
    dplyr::filter(!is.na(EC_number)) %>%
    dplyr::filter(confidence == "0" | confidence == "4") %>%
    # Exclude non-mitochondrially localized reactions from TCA.
    mutate(subsystem = case_when(
      reaction_id == "SPMDOX_pos" ~ "Arginine and Proline Metabolism",
      subsystem == "Citric acid cycle" & !grepl("[m]", formula, fixed = TRUE) ~ "Other",
      TRUE ~ subsystem)) %>%
    # Assign reactions to the appropriate subsystem.
    mutate(subsystem_priority = factor(subsystem) %>%
        fct_recode("Glycolysis" = "Glycolysis/gluconeogenesis",
          "TCA cycle" = "Citric acid cycle") %>%
        fct_collapse("Amino acid metabolism" = c(
          "Alanine and aspartate metabolism",
          "Arginine and Proline Metabolism",
          "beta-Alanine metabolism",
          "Cysteine Metabolism",
          "D-alanine metabolism",
          "Folate metabolism",
          "Glutamate metabolism",
          "Glycine, serine, alanine and threonine metabolism",
          "Histidine metabolism",
          "Lysine metabolism",
          "Methionine and cysteine metabolism",
          "Taurine and hypotaurine metabolism",
          "Tryptophan metabolism",
          "Tyrosine metabolism",
          "Urea cycle",
          "Valine, leucine, and isoleucine metabolism")) %>%
        fct_other(keep = facets) %>%
        fct_relevel(facets)) %>%
    # Keep only the subsystems for which we want to plot a facet.
    dplyr::filter(subsystem_priority != "Other") %>%
    # Lower-bound the adjusted p-value.
    mutate(adjusted_p_value = if_else(
      subsystem_priority == "Amino acid metabolism" & adjusted_p_value <= 1e-12,
      1e-12,
      adjusted_p_value)) %>%
    # Assign descriptive labels to various reactions.
    mutate(label = case_when(
      reaction_id == "PGM_neg" ~ "phosphoglycerate mutase (PGAM)",
      reaction_id == "LDH_L_neg" ~ "lactate dehydrogenase",
      reaction_id == "PDHm_pos" ~ "pyruvate dehydrogenase (PDH)",
      reaction_id == "TPI_neg" ~ "triosephosphate isomerase (DHAP forming)",
      reaction_id == "FACOAL1821_neg" ~ "long-chain fatty-acid-CoA ligase",
      reaction_id == "r1257_pos" ~ "long-chain fatty-acid-CoA ligase",
      reaction_id == "FACOAL1831_neg" ~ "long-chain fatty-acid-CoA ligase",
      reaction_id == "CSNATr_neg" ~ "carnitine O-acetyltransferase",
      reaction_id == "C160CPT1_pos" ~ "carnitine O-palmitoyltransferase",
      reaction_id == "ACONTm_pos" ~ "aconitate hydratase",
      reaction_id == "SUCOASm_pos" ~ "succinate-CoA ligase",
      reaction_id == "AKGDm_pos" ~ "alpha-ketoglutarate dehydrogenase",
      reaction_id == "SUCD1m_pos" ~ "succinate dehydrogenase",
      reaction_id == "ICDHyrm_pos" ~ "isocitrate dehydrogenase",
      reaction_id == "CK_pos" ~ "creatine\nkinase",
      reaction_id == "PGCD_pos" ~ "phosphoglycerate dehydrogenase",
      reaction_id == "ARGSS_pos" ~ "arginosuccinate synthase",
      reaction_id == "r0281_neg" ~ "putrescine diamine oxidase",
      reaction_id == "SPMDOX_pos" ~ "spermidine dehydrogenase (spermidine -> GABA)",
      reaction_id == "ARGDCm_pos" ~ "arginine decarboxylase",
      reaction_id == "AGMTm_pos" ~ "agmatinase",
      reaction_id == "GHMT2r_pos" ~ "serine hydroxymethyltransferase",
      reaction_id == "AHC_pos" ~ "adenosylhomocysteinase",
      reaction_id == "METAT_pos" ~ "methionine adenosyltransferase",
      reaction_id == "METS_pos" ~ "methionine\nsynthase",
      reaction_id == "ARGN_pos" ~ "arginase",
      TRUE ~ ""))
  
  ggplot(compass_scores_by_cell_type,
    aes(x = cohens_d, y = -log10(adjusted_p_value), color = subsystem_priority)) +
    ggtitle("Differential COMPASS Scores for Th17p vs. Th17n Cells") +
    xlab("Cohen's d") + ylab("-log(BH-adjusted p-value)") +
    xlim(-2.2, 2.2) +
    facet_wrap(vars(subsystem_priority), scales = "free_y", ncol = 2) +
    scale_color_manual(values = c(
      "Glycolysis" = "#662D8C",
      "TCA cycle" = "#B87013",
      "Fatty acid oxidation" = "#0B0D9D",
      "Amino acid metabolism" = "#B82130"
    )) +
    guides(color = FALSE) +
    geom_point(size = 1, alpha = 0.5) +
    geom_hline(yintercept = 1, linetype="dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype="dashed", color = "blue") +
    geom_text_repel(aes(label = label),min.segment.length = 0.1,point.padding = 0.5,
      size = 2,seed = 7) +theme_bw()
}


