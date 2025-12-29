# Low-Dimensional Embeddings for scRNA-seq Data

This project explores preprocessing and low-dimensional embeddings for single-cell RNA sequencing (scRNA-seq) data, with a focus on scalability and interpretability. The analysis is conducted using datasets from the CELLxGENE Census and standard single-cell analysis tools such as Scanpy, with exploratory integration of TorchDR for scalable dimensionality reduction.



## Project Scope and Motivation

Single-cell RNA-seq datasets are characterized by extremely high dimensionality, sparsity, and technical noise. Dimensionality reduction is therefore a crucial step for downstream analysis such as visualization and clustering. However, embedding results are highly sensitive to preprocessing decisions and algorithmic parameters.

This project investigates:
- how preprocessing and quality control affect embeddings,
- how classical methods (PCA, UMAP) behave on scRNA-seq data,
- and how scalable approaches (e.g. TorchDR, incremental PCA) can be integrated for larger datasets.

## Analysis Pipeline

The analysis follows a standard scRNA-seq workflow:

1. Data acquisition from CELLxGENE Census  
2. Exploratory data analysis and metadata inspection  
3. Quality control (cell- and gene-level filtering)  
4. Normalization and log-transformation  
5. Highly Variable Gene (HVG) selection  
6. Principal Component Analysis (PCA)  
7. Neighborhood graph construction  
8. Non-linear embeddings (UMAP)  
9. Exploratory integration of TorchDR for scalable dimensionality reduction  

## Repository Structure

```
scRNAseq-data-lowdimensional-embeddings/
├── README.md
├── environment.yml
├── .gitignore
├── data/                  
│   ├── raw/
│   └── processed/
├── src/
│   ├── __init__.py
│   ├── load_data.py
│   ├── dataset_summary.py
│   ├── preprocessing.py
│   ├── pca.py
│   ├── ipca_torchdr.py
│   └── variance_comparison.py
├── figures/
│   ├── qc_metrics.png
│   ├── PCA_pct_counts_mt.png
│   ├── PCA1v2_pct_count_mt_PCA3v4_n_genes.png
│   ├── HVG.png
│   ├── celldistributionbytissue.png
│   ├── pca/diagnostics/variance_ratio.png
│   ├── pca/projections/pca_scanpy.png
│   ├── ipca/diagnostics/varianceplot_ipca.png
│   ├── ipca/projections/ipca_torchdr.png
│   ├── Cumulative_explained_variance_PCAvsIPCA.png
│   ├── UMAP.png
│   └── t-SNE.png
└── docs/
    ├── progress_log.md
    └── meeting_notes.md
```

## Current Status

- [x] Project structure and environment setup  
- [x] Data loading and metadata inspection  
- [X] Exploratory data analysis and QC refinement 
- [X] PCA-based dimensionality reduction and interpretation
- [X] TorchDR integration
- [ ] Large-scale benchmarking (planned) 

All active development takes place on the `dev` branch. The `main` branch is kept as a clean, stable reference.

## Data

- Due to size constraints, raw and processed scRNA-seq data files (`.h5ad`) are not included in this repository.
- The directory structure is preserved, and all data can be re-generated using the provided scripts via the CELLxGENE Census API.

## Reproducibility and Environment

All experiments were conducted on a local machine using WSL2 (Ubuntu) with the following hardware configuration:

- CPU: Intel Core i7-10870H
- RAM: 32 GB
- GPU: NVIDIA GeForce RTX 3050, 4 GB VRAM
- CUDA: Version 13.0 (available via WSL2)
- OS: Windows 11 with Ubuntu (WSL2)

> [!NOTE]
> GPU acceleration was available and used where supported (TorchDR). However, the limited VRAM (4 GB) constrained the maximum dataset size, motivating the use of subsampling and incremental methods.

## Environment Setup

This project uses a Conda-based Python environment.

Main dependencies include:
- Python 3.11
- scanpy
- cellxgene-census
- numpy, pandas, matplotlib
- PyTorch (for TorchDR experiments)

To recreate the environment:

```bash
conda env create -f environment.yml
conda activate census
```

## Data Source

This project uses the CELLxGENE Census as the primary data source for large-scale single-cell RNA-seq datasets.
- The Census API provides programmatic access to curated scRNA-seq data stored in a cloud-backed TileDB-SOMA format.
- Please follow the [official installation instructions](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_installation.html)