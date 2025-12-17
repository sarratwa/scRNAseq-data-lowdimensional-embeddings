# Progress Log

This document tracks the technical and conceptual progress of the project.
Entries are ordered chronologically and summarize work done, issues encountered,
and next steps.

## After 13-10-2025

### Work done
- Began data acquisition using the CELLxGENE Census API.
- Followed the official tutorial to install `cellxgene-census` via pip.
- Identified that the Census API requires Python 3.11.
- Created a dedicated conda environment (`census`) with Python 3.11.
- Investigated installation failures related to native dependencies.
- Installed Windows Subsystem for Linux (WSL) with Ubuntu to enable Linux-based tooling.
- Installed Miniconda inside WSL and recreated the `census` environment.
- Successfully installed the Census API inside WSL.
- Verified installation by importing `cellxgene_census` and checking the version.

### Issues encountered
- Installation failed on Windows due to missing native build dependencies.
- The `tiledbsoma` package failed to build because required compiled C/C++ libraries (`tiledbsoma.dll`) are not available for Windows.
- Attempted installation via `conda-forge` failed, as `tiledbsoma` is not published for `win-64`.
- Initial Census queries caused kernel crashes due to excessive memory usage.
- Census API requires an explicit dataset version date to ensure data consistency and reproducibility.

### Resolutions
- Switched from native Windows to WSL (Linux environment) to resolve dependency issues.
- Installed `tiledbsoma` successfully via pip within WSL.
- Specified a Census snapshot version using `open_soma(census_version="2025.01.30")`.
- Reduced query size by limiting queries to a specific tissue (e.g. lung) to avoid RAM exhaustion.

### Next steps
- Implement more selective queries to control memory usage.
- Load queried data into AnnData format using `get_anndata`.
- Inspect downloaded metadata (`obs` and `var`) to understand dataset structure.
- Integrate data loading into the project codebase (`load_data.py`).

## After 27-10-2025

### Work done
- 

### Observations / issues
- 

### Next steps
- 

## After 27-11-2025

### Work done
- account for variance mean dispersion and fix plot and filter out cells
- focus on the DR section

## After 04-12-2025

### Work done
- research into:
    - PCA: preserve global variance ...
    - UMAP / t-SNE: preserve local neighborhoods ...
    - TorchDR: same math goal as PCA but is incremental and GPU-based

### Next steps
- studying the influence of these methods on the apperance of the structure, interpreting the result and verfiying its reproducibilty


## After 11-12-2025

### Work done
- 

### Observations / issues
- 

### Next steps
- 
