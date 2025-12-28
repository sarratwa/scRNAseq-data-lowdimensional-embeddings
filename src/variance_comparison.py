import scanpy as sc
import numpy as np
import torch
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from torchdr import IncrementalPCA

# Load data
adata = sc.read_h5ad("data/processed/brain_3000_preprocessed.h5ad")

X = adata.X
if not isinstance(X, np.ndarray):
    X = X.toarray()

# Global scaling
scaler = StandardScaler(with_mean=True, with_std=True)
X_scaled = scaler.fit_transform(X)

# Classical PCA
pca = PCA(n_components=30)
pca.fit(X_scaled)
pca_var = pca.explained_variance_ratio_

# Incremental PCA 
X_torch = torch.tensor(X_scaled, dtype=torch.float32, device="cuda")

ipca = IncrementalPCA(
    n_components=30,
    batch_size=256,
    device="cuda"
)

for i in range(0, X_torch.shape[0], 256):
    batch = X_torch[i:i+256]
    if batch.shape[0] >= ipca.n_components:
        ipca.partial_fit(batch)

ipca_var = ipca.explained_variance_ratio_.cpu().numpy()

#  Plot cumulative variance 
plt.figure(figsize=(5, 4))
plt.plot(np.cumsum(pca_var), label="PCA", marker="o")
plt.plot(np.cumsum(ipca_var), label="IPCA", marker="s")
plt.xlabel("Number of principal components")
plt.ylabel("Cumulative explained variance")
plt.title("Cumulative explained variance: PCA vs IPCA")
plt.legend()
plt.tight_layout()
plt.show()
