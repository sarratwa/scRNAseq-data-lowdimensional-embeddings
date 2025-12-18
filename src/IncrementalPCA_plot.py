import matplotlib.pyplot as plt
import scanpy as sc

adata = sc.read_h5ad("data/processed/brain_ipca.h5ad")

Z = adata.obsm["X_ipca"]

plt.figure(figsize=(6, 5))
plt.scatter(Z[:, 0], Z[:, 1], s=5)
plt.xlabel("IPCA 1")
plt.ylabel("IPCA 2")
plt.title("Incremental PCA (PC1 vs PC2)")
plt.show()

# ipca doesnt have the same structure as standard pca