#!/usr/bin/env python
### Integration
###-----------------------------------------------------------------------------###
DOC="Integration"
AUTHOR="Giovanna Maklouf"
###-----------------------------------------------------------------------------###

# load packages
import scvi
import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
import glob


indir = Path("data/paper/07_classifyCells/output/")
outdir = Path("/data/paper/08_integration/")
filenames = os.listdir(indir)
filenames = [file for file in os.listdir(indir) if os.path.splitext(file)[1] == '.h5ad']
filenames 

adatas = [sc.read(str(indir)+ "/" + i) for i in filenames]
adata = adatas[0].concatenate(adatas[1:], join = 'outer', fill_value = 0)
print(adata.shape)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(adata.shape)

adata.layers["counts"] = adata.X.copy()
adata.raw = adata
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata,
    flavor = "seurat_v3",
    n_top_genes = 3000,
    batch_key = "batch_key", 
    layer = 'counts',
    subset = True
)

scvi.model.SCVI.setup_anndata(adata,
                        layer = "counts",
                        batch_key = "batch_key")
vae = scvi.model.SCVI(
  adata,
  n_layers = 3,
  n_latent = 50,
  dropout_rate = 0.1
)
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()
print('scvi')

sc.tl.pca(adata, svd_solver = 'arpack')
print('pca')
sc.pp.neighbors(adata, use_rep = "X_scVI")
print('neighbors')

resolutions = np.arange(0.3, 2, 0.1)
for resolution in resolutions:
    key_added = f"leiden_{resolution:.1f}"
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)

sc.tl.umap(adata)
print('umap')

adata.write(str(outdir) + "/" + "integration.h5ad")
