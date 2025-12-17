import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

# need to limit scaling to hvgs and copy te adata before filtering

# -------------------------------
# 1. Load data
# -------------------------------
adata = sc.read_h5ad("data/processed/brain_3000_sample.h5ad")

# -------------------------------
# 2. Annotate mitochondrial genes
# -------------------------------
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# -------------------------------
# 3. Compute QC metrics
# -------------------------------
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    inplace=True
)

# -------------------------------
# 4. Filter low-quality cells
# -------------------------------
adata = adata[adata.obs["n_genes_by_counts"] > 200, :]
adata = adata[adata.obs["pct_counts_mt"] < 5, :]

# Filter genes detected in at least 3 cells
sc.pp.filter_genes(adata, min_cells=3)

# -------------------------------
# 5. Normalization + log-transform
# -------------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# -------------------------------
# 6. HVG detection (classical Seurat)
# -------------------------------
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat",
    n_top_genes=2000,
    inplace=True
)

# Extract mean + dispersion
means = adata.var["means"]
dispersions = adata.var["dispersions"]
hvg_mask = adata.var["highly_variable"]

# -------------------------------
# 7. Custom HVG PLOT 
# -------------------------------
plt.figure(figsize=(7,6))

# non-HVGs in grey
plt.scatter(
    means[~hvg_mask],
    dispersions[~hvg_mask],
    s=5,
    color="lightgrey",
    label="Other genes"
)

# HVGs in red
plt.scatter(
    means[hvg_mask],
    dispersions[hvg_mask],
    s=5,
    color="red",
    label="Highly variable"
)

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Mean Expression (log scale)")
plt.ylabel("Dispersion (log scale)")
plt.title("Highly Variable Gene Selection")
plt.legend(frameon=False)
plt.tight_layout()
plt.show()

# -------------------------------
# Scale HVGs (needed for PCA)
# -------------------------------
sc.pp.scale(adata, max_value=10)

# -------------------------------
# PCA
# -------------------------------
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color=None)

# Variance ratio plot (helps choose n_pcs)
sc.pl.pca_variance_ratio(adata, log=True)

sc.pl.pca(adata, color="pct_counts_mt")
sc.pl.pca(
    adata,
    color=["pct_counts_mt", "n_genes_by_counts"],
    dimensions=[(0, 1), (2, 3)],
    ncols=2,
    size=20,
)

# -------------------------------
# Neighborhood graph
# -------------------------------
# Use 30 PCs unless variance plot suggests fewer
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=10)

# -------------------------------
# UMAP
# -------------------------------
sc.tl.umap(adata)
sc.pl.umap(adata, color=["n_genes_by_counts", "pct_counts_mt"])

sc.tl.tsne(adata, n_pcs=30)
sc.pl.tsne(adata, color=["n_genes_by_counts"])

'''
# -------- Load data --------
adata = sc.read_h5ad("/home/sarra/brain_subsamples/brain_3000_sample.h5ad")

# optional: set verbosity & plotting style
sc.settings.verbosity = 3  # show warnings/info
sc.set_figure_params(figsize=(6, 4), facecolor='white')

# assume you already loaded your data:
# adata = sc.read_h5ad("...")

# 1. Annotate mitochondrial genes (human: "MT-")
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# 2. Calculate QC metrics: total counts, number of genes, percent mito, etc.
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],             # you can also add ribosomal ("ribo") or other gene sets if relevant
    percent_top=None,           # or e.g. [50,100,200] to compute library complexity stats  
    log1p=False,
    inplace=True
)

# 3. Plot QC metrics to inspect data and choose thresholds
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True
)
sc.pl.scatter(
    adata,
    x="total_counts",
    y="n_genes_by_counts",
    color="pct_counts_mt"
)
plt.show()

# 4. Filter cells and genes based on QC
# thresholds here are examples — adapt based on your data / QC plots
min_genes = 200       # cells must express at least 200 genes
max_genes = None      # optionally set an upper limit if needed
min_counts = None     # optionally set minimum counts
max_counts = None     # optionally set maximum total counts
max_pct_mito = 5      # exclude cells with >5% mitochondrial counts

sc.pp.filter_cells(adata, min_genes=min_genes)
if min_counts is not None:
    sc.pp.filter_cells(adata, min_counts=min_counts)
if max_counts is not None:
    sc.pp.filter_cells(adata, max_counts=max_counts)

adata = adata[adata.obs.pct_counts_mt < max_pct_mito, :]

# filter genes: keep only genes detected in a minimum number of cells
min_cells = 3
sc.pp.filter_genes(adata, min_cells=min_cells)

print(f"After QC filter: {adata.n_obs} cells, {adata.n_vars} genes")

# 5. Normalize counts per cell (library-size normalization) + log transform
sc.pp.normalize_total(adata, target_sum=1e4)  # e.g. scale each cell to 10,000 counts
sc.pp.log1p(adata)

# ---- HVG computation (classical method) ----
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat",    # This produces the classical HVG curve
    n_top_genes=2000,
    inplace=True
)

# ---- HVG Plot ----
sc.pl.highly_variable_genes(
    adata,
    show=True
)
'''
