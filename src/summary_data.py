import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -------- Load data --------
# only after running load_data.py and obtaining the h5ad file of the data
adata = sc.read_h5ad("data/raw/brain_3000_sample.h5ad")

print("============================================")
print("DATASET SUMMARY")
print("============================================")
print(f"Cells (n_obs):    {adata.n_obs:,}")
print(f"Genes (n_vars):   {adata.n_vars:,}")
print(f"Matrix format:    {type(adata.X).__name__}")
print(f"Layers available: {list(adata.layers.keys()) if hasattr(adata, 'layers') else 'None'}")
print(f"Metadata columns: {list(adata.obs.columns)}")
print(f"Gene annotation columns: {list(adata.var.columns)}")
print("============================================\n")

# -------- Tissue distribution (QC) --------
if "tissue" in adata.obs.columns:
    tissue_counts = adata.obs["tissue"].value_counts()

    print("Tissue distribution (full table):")
    print(tissue_counts)

else:
    tissue_counts = None


# -------- Top genes by mean expression --------
# Convert mean to 1D array
gene_means = np.asarray(adata.X.mean(axis=0)).ravel()

# Top 10 genes
top_gene_indices = gene_means.argsort()[-10:][::-1]
top_genes = adata.var.index[top_gene_indices].tolist()

print("\nTop 10 detected genes (by mean expression):")
for g in top_genes:
    print(" -", g)


# -------- Plot: Tissue distribution (Top N + Other) --------
if tissue_counts is not None:
    N_TOP = 30  # number of tissues to show

    top_tissues = tissue_counts.head(N_TOP)
    other_count = tissue_counts.iloc[N_TOP:].sum()

    plot_counts = top_tissues.copy()
    plot_counts["Other"] = other_count

    plt.figure(figsize=(6, 5))
    plot_counts.sort_values().plot(
        kind="barh",
        color="steelblue"
    )

    plt.xlabel("Number of Cells")
    plt.title(f"Cell distribution by tissue (top {N_TOP})")
    plt.tight_layout()
    plt.show()


# -------- Table view --------
summary_table = pd.DataFrame({
    "Metric": ["Cells", "Genes", "Matrix Type", "Metadata Columns"],
    "Value": [adata.n_obs, adata.n_vars, type(adata.X).__name__, len(adata.obs.columns)]
})

print("\n============================================")
print("Summary Table")
print("============================================")
print(summary_table)
