import scanpy as sc
import numpy as np
import torch
import matplotlib.pyplot as plt
from torchdr import IncrementalPCA

# ---------- Load preprocessed data -----------
adata = sc.read_h5ad("data/processed/brain_preprocessed.h5ad")

# Make sure we only use HVGs
adata = adata[:, adata.var.highly_variable].copy()

X = adata.X

# Convert sparse matrix to dense if needed
if not isinstance(X, np.ndarray):
    X = X.toarray()

# Convert to torch tensor
X_torch = torch.tensor(X, dtype=torch.float32)


# ---------- Incremental PCA (torchdr) --------------------
n_components = 30
batch_size = 512

ipca = IncrementalPCA(
    n_components=n_components,
    batch_size=batch_size,
    device="cpu",   # change to "cuda" if available
)

# Fit in batches
ipca.fit(X_torch)

# Transform full dataset
X_ipca = ipca.transform(X_torch).cpu().numpy()

# Store results in AnnData (Scanpy-compatible)
adata.obsm["X_ipca_torchdr"] = X_ipca


# --------- Variance explained --------------------
explained_var = ipca.explained_variance_ratio_.cpu().numpy()

pcs = np.arange(1, n_components + 1)

plt.figure(figsize=(6, 4))
plt.yscale("log")
plt.plot(
    np.arange(1, n_components + 1),
    explained_var,
    marker="o"
)
for i, var in enumerate(explained_var):
    plt.text(
        pcs[i],
        var,
        f"PC{pcs[i]}",
        fontsize=8,
        ha="center",
        va="bottom"
    )
plt.xlabel("Principal Component")
plt.ylabel("Explained Variance Ratio")
plt.title("Incremental PCA (torchdr)")
plt.tight_layout()
plt.show()


# ---------- Scatter plots (PC1 vs PC2) --------------------
plt.figure(figsize=(6, 5))
plt.scatter(
    X_ipca[:, 0],
    X_ipca[:, 1],
    s=5,
    alpha=0.6
)
plt.xlabel("IPCA 1")
plt.ylabel("IPCA 2")
plt.title("Incremental PCA (torchdr)")
plt.tight_layout()
plt.show()


# -------------- Optional: compare to Scanpy PCA -----------------
if "X_pca" in adata.obsm:
    X_pca = adata.obsm["X_pca"]

    plt.figure(figsize=(10, 4))

    plt.subplot(1, 2, 1)
    plt.scatter(X_pca[:, 0], X_pca[:, 1], s=5)
    plt.title("Scanpy PCA")

    plt.subplot(1, 2, 2)
    plt.scatter(X_ipca[:, 0], X_ipca[:, 1], s=5)
    plt.title("torchdr Incremental PCA")

    plt.tight_layout()
    plt.show()