import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

# need to limit scaling to hvgs and copy te adata before filtering

# --------- Load Data --------------
adata = sc.read_h5ad("data/processed/brain_preprocessed.h5ad")

# ------ PCA -------------
# we reduce these thousand of genes to components
# what does the package does underneath?
# should be done manually? 
sc.tl.pca(adata, svd_solver='arpack')
# pca scatter plot
# scanpy hides tick labels by default to keep plots clean for multi-panel layouts.
sc.pl.pca(
    adata,
    annotate_var_explained=True,
    title="PCA Scanpy",
)

# --------- Variance ratio plot (helps choose n_pcs) --------
# pc number and the % variance explained
# which PCs explain little additional variance -> we will use 30 PCs
sc.pl.pca_variance_ratio(adata, log=True)
# these 2 plots are to check if PCs are strongly correlated with mitochondrial %
sc.pl.pca(adata, color="pct_counts_mt")
sc.pl.pca(
    adata,
    color=["pct_counts_mt", "n_genes_by_counts"],
    dimensions=[(0, 1), (2, 3)],
    ncols=2,
    size=20,
)

# --- Neighborhood graph -------
# Use 30 PCs because from PC1 to PC10 is a strong drop and then from PC10 to PC20 is gradual decline and until PC30 we still have a non zero signal
# choosing less might make neighbor graphs incomplete
# each cell is connected to its nearest neighbors
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=10)

# --- UMAP --------
sc.tl.umap(adata)
sc.pl.umap(adata, color=["n_genes_by_counts", "pct_counts_mt"])
# ------t-SNE------------
sc.tl.tsne(adata, n_pcs=30)
sc.pl.tsne(adata, color=["n_genes_by_counts"])