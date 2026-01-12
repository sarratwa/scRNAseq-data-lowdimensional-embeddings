import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import time
import psutil
import os

process = psutil.Process(os.getpid())

# measuring rss: Resident Set Size: how much physical system RAM the operating system has allocated to the Python process
def cpu_mem_mb():
    return process.memory_info().rss / 1024**2

# need to limit scaling to hvgs and copy te adata before filtering

# --------- Load Data --------------
SAMPLE_SIZE = 10000 # 3000 or 10000
adata = sc.read_h5ad(f"data/processed/brain_{SAMPLE_SIZE}_preprocessed.h5ad")

#before pca
mem_before = cpu_mem_mb()
# start timer
start = time.perf_counter()

# ------ PCA -------------
# we reduce these thousand of genes to components
# what does the package does underneath?
# should be done manually? 
sc.tl.pca(
    adata,
    svd_solver="arpack",
    n_comps=30,
    mask_var="highly_variable"
)

# stop timer
pca_time = time.perf_counter() - start
mem_after = cpu_mem_mb()
# results
print(f"Classical PCA time: {pca_time:.2f}s")
print(f"CPU RAM increase: {mem_after - mem_before:.1f} MB")

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
sc.pl.pca_variance_ratio(adata, log=False)
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