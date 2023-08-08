#!/usr/bin/env python
### Subsets
###-----------------------------------------------------------------------------###
DOC="Subsets"
AUTHOR="Giovanna Maklouf"
###-----------------------------------------------------------------------------###

import scvi
import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
import glob


def FilterBatchByValue(adata, threshold):
    batch = pd.DataFrame(adata.obs.batch_key.value_counts())
    batch = batch.reset_index()
    batch_rows = batch[batch['count'] > threshold].batch_key
    adata = adata[adata.obs['batch_key'].isin(batch_rows)]
    return adata


adata = sc.read("/data/integration.h5ad")

for celltype in adata.obs['celltype'].unique():
    
    adata_subset = adata[adata.obs['celltype'] == celltype].copy()
    adata_subset = adata_subset.raw.to_adata()
    adata_subset.obs["batch_key"] = adata_subset.obs["batch_key"].astype(str) 
    adata_subset = FilterBatchByValue(adata_subset, 10)

    adata_subset.layers["counts"] = adata_subset.X.copy()
    adata_subset.raw = adata_subset  # manter a dimensão completa segura
    sc.pp.normalize_total(adata_subset, target_sum=1e4)
    sc.pp.log1p(adata_subset)
    sc.pp.highly_variable_genes(
        adata_subset,
        flavor="seurat",
        n_top_genes=3000,
        batch_key="batch_key",
        subset=True
    )

    scvi.model.SCVI.setup_anndata(adata_subset, layer="counts", batch_key="batch_key")

    vae = scvi.model.SCVI(adata_subset, dispersion="gene-batch", dropout_rate=0.1)

    vae.train()

    adata_subset.obsm["X_scVI"] = vae.get_latent_representation()

    sc.tl.pca(adata_subset, svd_solver='arpack')
    sc.pp.neighbors(adata_subset, use_rep="X_scVI")

    # Loop para iterar sobre os valores de resolução
    resolutions = np.arange(0.3, 2, 0.1)
    for resolution in resolutions:
        key_added = f"leiden_{resolution:.1f}"  # Definindo o nome da chave com base no valor de resolução
        sc.tl.leiden(adata_subset, resolution=resolution, key_added=key_added)

    sc.tl.umap(adata_subset)

    # Salvar os resultados
    output_file = f"/data/paper/09_subsets/anndata/{celltype.lower().replace(' ', '_')}_subset.h5ad"
    adata_subset.write(output_file)
