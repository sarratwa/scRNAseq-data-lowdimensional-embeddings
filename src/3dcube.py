import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import scanpy as sc

# get the data
SAMPLE_SIZE = 3_000 # change to 10_000 or 3_000
adata = sc.read_h5ad(f"data/raw/brain_{SAMPLE_SIZE}_sample.h5ad")

# Pick a small slice: 20 cells × 20 genes
n_cells = 200
n_genes = 200

cell_idx = np.random.choice(adata.n_obs, n_cells, replace=False)
gene_idx = np.random.choice(adata.n_vars, n_genes, replace=False)

# print(adata.obs)

matrix_block = adata.X[cell_idx][:, gene_idx].toarray()

# 2D Grid (cells × genes)

fig, ax = plt.subplots(figsize=(6, 5))

# Draw grid
for x in range(n_cells + 1):
    ax.axvline(x, color="black", linewidth=0.5)
for y in range(n_genes + 1):
    ax.axhline(y, color="black", linewidth=0.5)

ax.set_xlim(0, n_cells)
ax.set_ylim(0, n_genes)
ax.invert_yaxis()

# Axis labels
ax.set_xlabel("Cells")
ax.set_ylabel("Genes")

# No ticks
ax.set_xticks([])
ax.set_yticks([])

plt.title(f"2D Cells × Genes Grid (from_{SAMPLE_SIZE}_dataset)")
plt.tight_layout()
plt.show()

# Conceptual 3D Cube
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='3d')

# Cube corners
r = [0, 1]
vertices = np.array([[x, y, z] for x in r for y in r for z in r])

# Cube faces
faces = [
    [vertices[j] for j in [0,1,3,2]],
    [vertices[j] for j in [4,5,7,6]],
    [vertices[j] for j in [0,1,5,4]],
    [vertices[j] for j in [2,3,7,6]],
    [vertices[j] for j in [0,2,6,4]],
    [vertices[j] for j in [1,3,7,5]],
]

# Draw cube
ax.add_collection3d(Poly3DCollection(faces, alpha=0.30, edgecolor='black'))

# Remove numeric tick labels only
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

# BUT keep grid & ticks (just invisible labels)
ax.set_xticks([0, 0.5, 1])
ax.set_yticks([0, 0.5, 1])
ax.set_zticks([0, 0.5, 1])

# Axis labels
ax.set_xlabel("Cells")
ax.set_ylabel("Genes")
ax.set_zlabel("Feature dimension")

plt.title("Conceptual 3D Cube")
plt.tight_layout()
plt.show()


# Expression Landscape (3D Surfacea)

X, Y = np.meshgrid(range(n_cells), range(n_genes))
Z = matrix_block  # expression slice

fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, Y, Z, cmap="viridis", edgecolor="none", alpha=0.9)

ax.set_xlabel("Cells")
ax.set_ylabel("Genes")
ax.set_zlabel("Expression level")

plt.title(f"Expression Landscape (_{SAMPLE_SIZE}_ Data Slice)")
plt.tight_layout()
plt.show()