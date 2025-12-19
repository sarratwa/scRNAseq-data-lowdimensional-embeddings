import scanpy as sc
import torch
import numpy as np
import torchdr

# checking if torchdr is installed
# print(torch.cuda.is_available())
# print(torchdr.__version__)

adata = sc.read_h5ad("data/processed/brain_preprocessed.h5ad")
# check if data as loaded 
# print(adata.obs)
# verify if pca exists
# print(adata.obsm.keys())
# print("PCA exists before saving:", "X_pca" in adata.obsm)
# print("obsm keys before saving:", adata.obsm.keys())
# fixed this by renaming the feature_name coloumn

#-------- Make sure PCA exists ---------
if "X_pca" not in adata.obsm:
    print("PCA missing — computing PCA in dr.py")
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")

# --------- Prepare PCA imput ------------
# i should not be working with the pca in the preprocessing file
# the data should recalculate the same output of the preprocessing -> transformed and normalised and not pca
# this was running a pca on a pca -> change this
X = adata.obsm["X_pca"][:, :30].astype(np.float32)

device = "cuda" if torch.cuda.is_available() else "cpu"
X_t = torch.from_numpy(X).to(device)

# to check if i have IncrementalPCA installed
# print(dir(torchdr))
# print([name for name in dir(torchdr) if "PCA" in name])

# --------- Prepare data matrix --------------------
# Use log-normalized expression
# data orientation? plot cells and not pca for genes
# 
X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
X = X.astype(np.float32)

device = "cuda" if torch.cuda.is_available() else "cpu"
X_t = torch.from_numpy(X).to(device)

# ------ TorchDR Incremental PCA ----------------------
ipca = torchdr.IncrementalPCA(
    n_components=50,
    batch_size=256
)

Z_t = ipca.fit_transform(X_t)

# -------  Store result in AnnData -------------------
adata.obsm["X_ipca"] = Z_t.cpu().numpy()

# ------ Save --------------------
adata.write("data/processed/brain_ipca.h5ad")

print("Incremental PCA shape:", adata.obsm["X_ipca"].shape)
# Incremental PCA shape: (2818, 50)
# compressed succesfully


'''

# seems to need FAISS: Facebook AI Similarity Search ? 
# __init__() got an unexpected keyword argument 'use_faiss'
# ------------- Torchdr embedding ----------
reducer = torchdr.UMAP(
    n_components=2,
    use_faiss=False  # 🚨 CRITICAL
)

Z_t = reducer.fit_transform(X_t)
adata.obsm["X_torchdr"] = Z_t.cpu().numpy()

# --------plot ----------
sc.pl.embedding(
    adata,
    basis="torchdr",
    color=["cell_type", "pct_counts_mt"],
    legend_loc="on data"
)
'''
