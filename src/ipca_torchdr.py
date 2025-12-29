import scanpy as sc
import numpy as np
import torch
import matplotlib.pyplot as plt
from torchdr import IncrementalPCA
import time
import torch

# GPU memory instrumentation
def gpu_mem(prefix=""):
    if torch.cuda.is_available():
        allocated = torch.cuda.memory_allocated() / 1024**2
        reserved = torch.cuda.memory_reserved() / 1024**2
        max_alloc = torch.cuda.max_memory_allocated() / 1024**2
        print(
            f"{prefix} | "
            f"allocated={allocated:.1f} MB | "
            f"reserved={reserved:.1f} MB | "
            f"max_alloc={max_alloc:.1f} MB"
        )

# ---------- Load preprocessed data -----------
SAMPLE_SIZE = 3000 # 3000 or 10000
adata = sc.read_h5ad(f"data/processed/brain_{SAMPLE_SIZE}_preprocessed.h5ad")

# Make sure we only use HVGs
adata = adata[:, adata.var.highly_variable].copy()

X = adata.X

# Convert sparse matrix to dense if needed
if not isinstance(X, np.ndarray):
    X = X.toarray()

# Convert to torch tensor
X_torch = torch.tensor(X, dtype=torch.float32)

# move data to GPU
X_torch = X_torch.cuda()

# start timer 
start_total = time.perf_counter()

# GPU memory instrumentation before torchDR
batch_times = []
torch.cuda.reset_peak_memory_stats()
gpu_mem("Before IPCA fit")

# ---------- Incremental PCA (torchdr) --------------------
# batch size
n_components = 30
# change batch sizes: 256, 512, 1024
batch_size = 256

ipca = IncrementalPCA(
    n_components=n_components,
    batch_size=batch_size,
    device="cuda",   # change to "cuda" if available
)

# Fit in batches but this hides batch times
# ipca.fit(X_torch)

# batch-timed IPCA implementation
for start_idx in range(0, X_torch.shape[0], batch_size):
    end_idx = min(start_idx + batch_size, X_torch.shape[0])
    batch = X_torch[start_idx:end_idx]

    if batch.shape[0] < ipca.n_components:
        break

    torch.cuda.synchronize()
    t0 = time.perf_counter()

    ipca.partial_fit(batch)

    torch.cuda.synchronize()
    batch_times.append(time.perf_counter() - t0)

# Transform full dataset
X_ipca = ipca.transform(X_torch)

# Move result back to CPU for numpy / scanpy
X_ipca = X_ipca.cpu().numpy()

# stop timer
ipca_time = time.perf_counter() - start_total
print(f"Incremental PCA time: {ipca_time:.2f}s")

# GPU memory instrumentation after torchDR
gpu_mem("After IPCA fit")

# result
print(f"Total IPCA time: {ipca_time:.2f}s")
print(f"Mean batch time: {np.mean(batch_times)*1000:.2f} ms")
print(f"Std batch time: {np.std(batch_times)*1000:.2f} ms")
print(f"Number of batches: {len(batch_times)}")

gpu_mem("After IPCA fit")

# Store results in AnnData (Scanpy-compatible)
adata.obsm["X_ipca_torchdr"] = X_ipca


# --------- Variance explained --------------------
explained_var = ipca.explained_variance_ratio_.cpu().numpy()

pcs = np.arange(1, n_components + 1)

plt.figure(figsize=(6, 4))
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
pc1_pct = explained_var[0] * 100
pc2_pct = explained_var[1] * 100

plt.xlabel(f"IPCA 1 ({pc1_pct:.2f}%)")
plt.ylabel(f"IPCA 2 ({pc2_pct:.2f}%)")
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