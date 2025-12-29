import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

# need to limit scaling to hvgs and copy te adata before filtering

# --------- Load Data --------------
SAMPLE_SIZE = 10_000 # 3_000 or 10_000
adata = sc.read_h5ad(f"data/raw/brain_{SAMPLE_SIZE}_sample.h5ad")

# ---- Annotate mitochondrial genes --------------
# focus on mitochondrial genes mt to indicate cell stress or death -> huge percentage indicates a low-quality cell
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# ------ Compute QC metrics ---------------
# these qc statistics are essential for filtering and they add 
# adata.obs
# ------- n_genes_by_counts (how many genes are detectedin each cell)
# ------- total_counts (total UMIs(Unique Molecular Identifier) per cell )
# ------- pct_counts_mt (% of counts from mitochondrial genes)
# adata.var
# -------- n_cells_by_counts (how many cells express this gene)

# redirect scanpy to the actual adata var names
adata = adata.copy()  # important if this is a view

adata.var_names = adata.var["feature_name"].astype(str)
adata.var_names_make_unique()

# need to rename the column of feature name because it raised a conflict later in dr: ValueError: DataFrame.index.name ('feature_name') is also used by a column whose values are different. This is not supported. Please make sure the values are the same, or use a different name.
adata.var.rename(columns={"feature_name": "gene_symbol"}, inplace=True)

# First let s mark the mitochondrial genes, "MT-" for human
adata.var["mt"] = adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    inplace=True
)

# Violin plots: univariate QC metrics
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True
)

# does this dataset not have %mt?
# print(adata.var_names[:10])
# print(adata.var.head())
# my mistake was using the scanpy adata.var_names but the gene names existed in adata.var[feature_name]

# ------ Filter low-quality cells -------------------
# cells with less than 200 genes are likely empty droplets
# cells with more than 5% mt RNA are dying 
adata = adata[adata.obs["n_genes_by_counts"] > 200, :]
adata = adata[adata.obs["pct_counts_mt"] < 5, :]

#----------- Filter genes detected in at least 3 cells ---------- 
# if these genes are expressed in fewer than 3 cells then they could be rare genes and they add noise 
sc.pp.filter_genes(adata, min_cells=3)

# ----- Normalization + log-transform ---------------------
# each cell is scaled to have 10k total counts to remove difference in sequencing depth
# target sum here 1r4 = 10k is just a random number in order to scale to the same total
# different normalisation -> incremental pca should not care about the feature names -> possibly some gene have the same name ? 
# which genes symbole are there twice. 
# which renaming caused the duplication
# go back to the duplication reason 
# coloumn names should be unique and this should not bother torchdr nor should it care about it
sc.pp.normalize_total(adata, target_sum=1e4)
# the log transformation compresses extreme value -> make distribution more Gaussian
sc.pp.log1p(adata)

# ------ HVG detection (classical Seurat) ------------
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat",
    n_top_genes=2000,
    inplace=True
)

# -------  Extract mean + dispersion ----------
# select top 2000 genes with highest normalized dispersion
means = adata.var["means"]
dispersions = adata.var["dispersions"]
hvg_mask = adata.var["highly_variable"]

# ----- Custom HVG PLOT --------------------
plt.figure(figsize=(7,6))

# --------- non-HVGs in grey ------------
plt.scatter(
    means[~hvg_mask],
    dispersions[~hvg_mask],
    s=5,
    color="lightgrey",
    label="Other genes"
)

# --------- HVGs in red ----------
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

# ------- Scale HVGs (needed for PCA)-----------------
# For scaling we make sure each gene has mean = 0 and variance = 1
sc.pp.scale(adata, max_value=10)
# what if i limit scaling to hvgs?
# sc.pp.scale(adata[:, adata.var.highly_variable])

# saving the processed data to use ind dr
adata.write(f"data/processed/brain_{SAMPLE_SIZE}_preprocessed.h5ad")