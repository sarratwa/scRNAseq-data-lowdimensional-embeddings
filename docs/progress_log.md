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
- Loaded a pre-extracted brain subset (brain_3000_sample.h5ad) using sc.read_h5ad
- Worked on a reduced dataset to enable fast iteration and debugging of preprocessing steps
- Annotated mitochondrial genes based on gene name prefix (MT-)
- Computed standard QC metrics using sc.pp.calculate_qc_metrics
- Filtered low-quality cells based on:
    - minimum number of detected genes (n_genes_by_counts > 200)
    - mitochondrial read fraction (pct_counts_mt < 5)
- Filtered genes detected in fewer than 3 cells to remove uninformative features
-> These steps follow standard scRNA-seq QC practices and aim to remove technical noise before downstream analysis.
- Applied total-count normalization (target_sum = 1e4)
- Performed log-transformation using sc.pp.log1p
- Identified highly variable genes using Scanpy’s Seurat-style method
- Extracted mean expression and dispersion values from adata.var
- Created a custom mean–dispersion plot:
- non-HVGs shown in gray
- HVGs highlighted in red
- Used log–log scaling to better visualize the mean–variance relationship
- Computed PCA using the ARPACK solver
- Generated:
    - PCA scatter plots
    - variance ratio plots (log-scaled)
- Used these plots to guide the selection of the number of principal components
- Visualized PCA embeddings colored by QC metrics (pct_counts_mt, n_genes_by_counts) to identify potential technical structure

### Observations / issues
- The preprocessing pipeline is implemented following Scanpy best practices and tutorials.
- While the code runs correctly and produces expected outputs, my current understanding of what constitutes a “good” versus “problematic” HVG, PCA, or UMAP plot is still developing.
- At this stage, I rely on:
    - official Scanpy tutorials
    - consistency with recommended workflows
- Further work is needed to:
    - better interpret variance structures
    - identify preprocessing artifacts visually
    - justify filtering thresholds more confidently

### Next steps
- Understand why PCA and HVG selection behave as they do.
- Justify all filtering decisions.
- Distinguish clearly between:
- biological signal
- technical artifacts
- Use TorchDR and incremental PCA as scalability tools, not black boxes.
- Demonstrate interpretability, not just implementation.

## After 15-12-2025

### Work done
- Saved the fully processed dataset to data/processed/brain_preprocessed.h5ad.
- Evaluated correlations between principal components and QC metrics (mitochondrial percentage, gene counts).
- Constructed neighborhood graphs and computed UMAP and t-SNE embeddings for exploratory analysis.
- TorchDR Preparation
- Implemented initial TorchDR workflow (dr.py) to compare:
    - Standard PCA
    - Incremental PCA (TorchDR)
- Prepare a draft of the structure for the final report.

### Observations / issues
- Identified an incorrect PCA-on-PCA workflow: During initial TorchDR integration, PCA was mistakenly applied to an already PCA-transformed representation (X_pca). This resulted in a PCA-on-PCA workflow, which is conceptually incorrect, as dimensionality reduction methods must operate on the normalized gene expression space rather than on previously reduced embeddings. The pipeline was subsequently corrected to ensure that both classical PCA (Scanpy) and Incremental PCA (TorchDR) are applied independently to the same preprocessed expression matrix.
- Incremental PCA and classical PCA are vastly different. Go over code. 

### Next steps
- Refactor the pipeline to ensure TorchDR operates directly on the log-normalized expression matrix.
- TorchDR Incremental PCA integrated and under active refinement.
- establish Classical PCA baseline.
- Get Preprocessing pipeline to stable and reproducible state.
- clean comparison between Scanpy PCA and TorchDR Incremental PCA under controlled data sizes.
- Start writing the final report. 

## After 18-12-2025 (Last meeting took place: deadline is on the 29th)

### Work done
- Draft2 for paper
- Draft3 for paper
- total runtime, peak memory CPU or GPU cuda for pca and ipca and batch size sensitivity for ipca. explained variance ratio? 
- Data scaling: 3k cells done, now downloading 10k cells, this time sampling runtime, CPU, RAM and final matrix size (cells x genes). this can become a table

### Observations / issues
- concerning draft3 must change/improve:
    - Explicit Research Question: This study investigates whether incremental PCA can reproduce the variance structure of classical PCA while reducing memory requirements under realistic computational constraints. (Introduction)
    - Explained Variance Results: comparison PCA vs IPCA / a cumulative explained variance plot (PCA vs IPCA), or a small table: PC1–PC10 variance (PCA vs IPCA)
    - Add memory usage section
    - Batch sizes that was used in ipca must be explained and described (found in ipca_torchdr.py). Answer questions what batch size did we use and why? 
    - Add Quantitative PCA Similarity Metric. example: PC1 in PCA vs IPCA
    - Runtime comparasion: runtime PCA vs IPCA (in result section)
    - Add UMAP and t-SNE figure
    - Showing HVGs is not enough, must explain how many HVGs and why were they retained.(preprocessing pipeline)
    - Add table for preprocessing summary and for PCA vs IPCA comparison
    - Computational envi description can be improved and add RAM, CPU and GPU specs.
    - explain how TorchDR is the same math as PCA but different approach: TorchDR does not introduce a new dimensionality reduction objective but implements PCA under a different computational model.
    - Prof talked about protein coding gene restriction, see what can be done about that
    - talk about the effect of batches
    - Explain why PCA is not about plots
    - Cite the references correctly
- After running peak memory usagem found that IPCA is slower than PCA (2,52s vs 1.75s). The goal now is to find out when does PCA runs out of memory and becomes slower than IPCA. The idea here is that PCA memory usage increases with the amount of cells whereas IPCA memory usage stays relatively the same. From the testing on the 3k cell dataset, we ve found that smaller batch sizes have been (batchsize:256 , num of batches: 11, runtime: 2.52s / batchsize:512 , num of batches: 6, runtime: 4.10s / batchsize: 1024, num of batches: 3, runtime: 6.43s). A plot would be intersting to see: batch size vs runtime.  
- Datascaling first attempt failed (30min): Timeout was reached. S3: Failed to read S3 object. 
- possible solution if second attempt fails: reducing concurency for TileDB-SOMA: allowing TileDB-SOMA to use fewer threads (limit the number of simultaneous remote storage requests). "sm.num_reader_threads": 2
- second data scaling failed: third datascaling with possible solution implemented.
- third data scaling failed: curlCode: 28, Timeout was reached. The connection stayed open too long S3 or the network closed it due to a timeout. So reducing concurrency did change much. At ~10k cells, Census downloads are unreliable over standard network connections.
- Conclusion after 3 failed attempts: the connection stays open too long S3 or this network drops it. the curlcode means The server did not finish sending the data before the connection expired. Trying chunking: downloading the same data, but in several small, independent pieces instead of one big request. https://stackoverflow.com/questions/65130143/how-to-processes-the-extremely-large-dataset-into-chunks-in-python-pandas-whi 
- chunking also failed: curlCode: 28, Timeout was reached. I attempted to scale data acquisition to 10,000 cells using the CellxGene Census. However, repeated network timeouts occurred due to remote streaming limitations, even when using chunked downloads. Therefore, I focused my benchmarking on a stable 3,000-cell dataset and analyzed memory and runtime behavior of PCA versus incremental GPU-based methods under controlled conditions.
- retry datascaling for 10k from a home network: 
    - Unexpectedly, a full download of 10,000 brain cells succeeded when executed from a home network environment, whereas previous attempts on other networks repeatedly failed. The same codebase and query parameters were used, suggesting that the success was not due to a code change but rather to differences in network conditions.
- specifiy : Gene expression matrices were obtained from the RNA measurement of the CellxGene Census and in preprocessing add: Analyses were restricted to protein-coding genes and identical preprocessing steps were applied across all dataset sizes. Only after implementing it. 
- PCA vs IPCA memory performance. 

### Observation PCA vs IPCA runtime
- PCA runtime:
    - Run1: 3,2s
    - Run2: 1.2s
    - Run3: 1,6s
    - Run4: 1,4s
    -> First run is slow because CPU cash is empty: The first PCA run exhibits higher runtime due to cold-start overheads such as memory allocation and library initialization. Subsequent runs benefit from cache reuse and show reduced execution times.
    - Low CPU RAM increase (~28–29 MB)
- IPCA runtime: 
    - First run is faster than normal PCA: 2,23s
    - Run 2, 3 and 4 are ~2.0 s
    - GPU memory:
        - allocated ~35 MB
        - reserved ~98 MB
        - max_alloc ~65 MB
    - Batch timing consistent (~170–195 ms)
- This has a number of implications: 
    - for small datasets 3k, classical PCA is faster or equal to IPCA
    - PCA sensitive to cache state
    - PCA has lower overhead for small N
    - IPCA has higher constant overhead due to batching, repeated partial SVD updates, GPU memory reservation
    - IPCA is scalable because: 
        - batch processing dominates
        - GPU execution is more deterministic
        - less reliance on CPU cache state
    - Memory behavior is fundamentally different   
        - CPU RAM: PCA:~28 MB / IPCA:minimal
        - GPU RAM: PCA:none / IPCA:~35–98 MB GPU
- For moderate scRNA-seq datasets (~3k cells), classical PCA remains efficient. Incremental PCA becomes advantageous primarily in memory-constrained or large-scale settings.
- reattempting data scaling 19:37

### Observation on the effect of batche sizes in IPCA
We evaluated batch sizes {256, 512, 1024} (2⁸, 2⁹, 2¹⁰) to probe the tradeoff between per-batch computation cost and the overhead of performing more incremental updates, using standard power-of-two batch sizes commonly used in GPU workloads.
- Batch size 256:
    - Total time ≈ 2 s
    - Mean batch time ≈ 633–649 ms
    - Batches = 11
    - GPU memory: allocated ≈ 22 MB, reserved ≈ 98 MB, max_alloc ≈ 65 MB
- Batch size 512:
    - Total time ≈ 3.84–3.98 s
    - Mean batch time ≈ 633–649 ms
    - Batches = 6
    - GPU memory: allocated ≈ 34.6 MB, reserved ≈ 104 MB, max_alloc ≈ 83.5 MB
- Batch size 1024: 
    - Total time ≈ 5.90–5.98 s
    - Total time ≈ 5.90–5.98 s
    - batches = 3
    - GPU memory increased: allocated ≈ 42.4 MB, reserved ≈ 170 MB, max_alloc ≈ 137 MB

### Observation on performance of IPCA on a 3k vs 10k cells dataset
- Classical PCA:
    - Classical PCA time: 4.25s
    - CPU RAM increase: 44.0 MB
- IPCA: batchsize=1024
    - Incremental PCA time: 30.57s
    - Number of batches: 10
    - Mean batch time: 3049.09 ms
    - After IPCA fit | allocated=84.1 MB | reserved=292.0 MB | max_alloc=201.9 MB
- IPCA: batchsize=256
    - After IPCA fit | allocated=83.5 MB | reserved=252.0 MB | max_alloc=201.3 MB
    - Total IPCA time: 6.65s
    - Mean batch time: 177.65 ms
    - Number of batches: 37
- The successful download enabled benchmarking of both classical PCA and incremental PCA on a 10k-cell dataset. This larger dataset revealed clear scaling effects: classical PCA remained faster as long as the dataset fit in CPU memory, while incremental PCA exhibited substantially higher runtime due to batching and GPU overhead, but with predictable and bounded memory behavior. These results reinforce the conclusion that incremental PCA is primarily advantageous in memory-constrained or larger-scale settings, rather than as a drop-in faster alternative for moderate dataset sizes.