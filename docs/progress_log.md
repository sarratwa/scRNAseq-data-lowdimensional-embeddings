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
- Continued development within WSL (Ubuntu) due to prior compatibility issues on native Windows.
- Successfully configured a Conda environment (census, Python 3.11) and installed the CELLxGENE Census API. 
- Verified correct data access by loading Census data into an AnnData (adata) object, confirming: correct cell count and gene count
- availability of adata.obs (cell metadata) and adata.var (gene metadata)
- Explored available metadata fields (e.g. assay, tissue, cell type, disease annotations).
- Successfully generated and stored a mouse lung subset (.h5ad) using WSL as a technical validation step for large data handling and file persistence.
- Attempted migration of the workflow to Google Colab to explore scalability and cloud-based execution.
- Tested API-based dataset selection in Colab using CELLxGENE dataset IDs, which allowed partial access without immediate crashes.

### Observations / issues
- Google Colab memory limitations:
    - Kernel crashes occurred when accessing large metadata tables or attempting to load full expression matrices.
    - RAM usage spiked even during metadata inspection (e.g. adata.obs.columns).
- Rejected approaches based on feedback:
    - Pre-downloaded .h5ad files (e.g. mouse data) were rejected as the project focus is on human data.
    - Manual downloads of human datasets were rejected in favor of fully API-driven workflows.
    - Dataset-ID-based loading in Colab was considered only partially automated and therefore insufficient.
- Infrastructure constraints:
    - Current local machine does not support large-scale data loading or GPU acceleration.
    - Cloud-based environments available so far were unstable for this dataset size.
- Overall, these limitations prevented stable end-to-end execution on large datasets in Colab at this stage.

### Next steps
- Continue development locally in WSL using small subsets (3k–5k cells) to:
    - complete the full pipeline (loading → preprocessing → PCA → embeddings)
    - ensure methodological correctness and conceptual understanding
- Document hardware specifications (CPU, RAM, GPU availability) to clearly contextualize computational constraints.
- Implement preprocessing steps in Scanpy:
    - cell filtering
    - gene filtering
    - highly variable gene selection
- Generate baseline dimensionality reduction results (PCA, UMAP) on small subsets.
- Prepare the pipeline for scalability testing, with the intention to:
    - migrate the same code to GPU-enabled infrastructure (e.g. KI-Werkstatt)
    - experiment with TorchDR, incremental PCA, and batch-based processing once stable.
- Discuss with the supervisor:
    - appropriate infrastructure for large-scale runs
    - realistic GPU usage expectations for TorchDR in this project.

## After 06-11-2025
### Work done
- Continued working locally in WSL (Ubuntu) to ensure stability and full control over memory usage.
- Implemented an automated, API-driven data loading pipeline using the CELLxGENE Census API.
- Wrote a dedicated data-loading script that:
    - connects to the Census using cellxgene_census.open_soma
    - queries human data only
    - filters cells by metadata (tissue_general == "brain" and is_primary_data == True)
- Retrieved cell-level metadata first (fast operation) and converted it to a pandas DataFrame to:
    - inspect available cells
    - avoid loading the full expression matrix prematurely
- Subsetting Strategy
- Adopted a random subsetting strategy to control computational load:
    - target subset size: 3,000 cells
    - sampling performed with a fixed random seed (random_state=42) for reproducibility
- Subsetting was applied at the cell level using soma_joinid values.
- Rationale for subsetting:
    - local hardware limitations
    - need to complete the full pipeline once before scaling
    - align with supervisor feedback to prioritize methodological correctness
- The 3k subset is explicitly treated as a development and validation dataset, not a final-scale experiment
- Expression Data Retrieval
    - Used human.axis_query(...) with an explicit obs_query to:
        - download expression data only for the selected subset of cells
        - avoid loading unnecessary data into memory
    - Converted the query result directly to an AnnData object using to_anndata, with:
    - X_name="raw" for raw counts
    - selected metadata columns (cell_type, tissue, donor_id, assay, disease, sex)
- Saved the resulting AnnData object as:
    - brain_3000_sample.h5ad
- This produced a self-contained, reusable dataset suitable for downstream preprocessing and dimensionality reduction.
- Initial Dataset Inspection and EDA
    - Loaded the generated .h5ad file using Scanpy.
    - Performed initial dataset inspection:
        - number of cells (adata.n_obs)
        - number of genes (adata.n_vars)
        - matrix type and available layers
        - available metadata columns in adata.obs and adata.var
    - Computed and inspected:
        - tissue distribution across the sampled cells
        - top genes by mean expression
        - Generated basic exploratory outputs:
            - printed summaries for logging and screenshots
            - bar plot of tissue distribution
            - summary table describing dataset dimensions and structure

### Observations / issues
- This local, subset-based workflow:
    - is fully automated and API-driven
    - avoids manual downloads
    - avoids unstable cloud environments
- The goal at this stage is pipeline validation, not biological inference.
- The current setup provides a stable baseline for:
    - preprocessing (normalization, filtering, HVG selection)
    - PCA and classical dimensionality reduction
    - later integration of TorchDR and scalability experiments

### Next steps
- Apply Scanpy preprocessing steps to the 3k subset:
    - normalization
    - log-transformation
    - gene and cell filtering
    - highly variable gene selection
- Generate and validate:
    - mean–variance (dispersion) plots
    - PCA embeddings
- Use this subset as a reference pipeline, then:
    - migrate the same code to larger subsets
    - test TorchDR, batching, and GPU-enabled execution once infrastructure allows

## After 13-11-2025
### Work done
- Writing the zwischenbericht according to feedback 

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
