#!/usr/bin/env python
### Integration
###-----------------------------------------------------------------------------###
DOC="Subsets"
AUTHOR="Gabriela Rapozo"
###-----------------------------------------------------------------------------###

import scanpy as sc
import matplotlib.pyplot as plt
import infercnvpy as cnv
import pandas as pd
from matplotlib import rcParams

adata = sc.read('/results/integration.h5ad')

adata = adata.raw.to_adata()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

cnv.io.genomic_position_from_gtf("/GRCh38/gtf/Homo_sapiens.GRCh38.83.gtf", adata=adata, gtf_gene_id='gene_name', inplace=True)
adata.var['chromosome'] = ['chr'+str(i) for i in adata.var['chromosome']]
adata.var['chromosome'].value_counts()

exclude_chromosomes = ['chrnan', 'chrMT', 'chrX', 'chrY', 'chrKI270728.1', 'chrKI270734.1', 'chrGL000194.1', 'chrKI270726.1', 
                       'chrKI270721.1', 'chrGL000195.1', 'chrKI270731.1', 'chrKI270711.1', 'chrGL000218.1', 'chrGL000219.1', 'chrGL000009.2']

reference_cat = ["T/NK Cells", "Mononuclear Phagocytes", 'Follicular B Cells', 'Plasma B Cells',
                 'Plasmacytoid Dendritic Cells', 'Mast Cells']

# Inferir CNV
cnv.tl.infercnv(
    adata,
    reference_key="celltype",
    exclude_chromosomes=exclude_chromosomes,
    reference_cat=reference_cat,
    window_size=250
)

cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
cnv.tl.cnv_score(adata)

adata_subset.write('/results/integration_cnv.h5ad')
