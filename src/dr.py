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
print(adata.obsm.keys())

