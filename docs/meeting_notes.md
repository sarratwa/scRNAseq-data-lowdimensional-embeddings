# Meeting Notes

This document summarizes meetings with the supervisor.
Each entry records feedback, decisions, and agreed action items.

---

## Meeting – 13-10-2025

**Focus of meeting**
- Introduction to the research thema

**Supervisor feedback**
- Dimensionality reduction is an early but critical step for:
    - batching
    - clustering
    - exploratory analysis
- The goal is to understand how to represent high-dimensional scRNA-seq data in low dimensions (e.g. 20k genes → 2D).
- Dataset source:
    - CELLxGENE
- Dataset characteristics:
    - scRNA-seq data
    - annotated cells
    - disease-related samples:
        - epilepsy
        - Alzheimer’s
        - brain disease distribution
- Possible tissues:
    - brain (focus)
    - liver or other tissues (mentioned as alternatives)
- Data access:
    - use APIs for the dataset
    - download data programmatically
- Data format:
    - annotations are available
    - cell metadata exists
- Tasks:
    - load the data
    - understand where the data comes from
    - inspect dataset metadata
- Perform exploratory statistics before modeling
- Goals:
    - understand dataset composition
    - inspect distributions
- Analyses:
    - cell type categories and counts
    - descriptive statistics
    - distribution of cells across categories
    - gene expression summaries
- Visualizations:
    - bar plots (e.g. brain vs liver, cell types)
    - frequency plots
- Tools:
    - Python
    - pandas DataFrames
- Questions:
    - what annotations exist?
    - are annotations reasonable?
    - how many cells per type?
    - how many genes?
- Gene expression analysis:
    - genes in common
    - frequency of genes in the expression matrix
- Inspect:
    - how often genes appear across cells
    - sparsity patterns
- Data handling:
    - represent gene expression in tabular form (pandas)
- Goal:
    - understand structure of the gene–cell matrix
- Filter the dataset:
    - filter cells
    - filter genes
- Reasons for filtering:
    - memory constraints
    - computational feasibility
- Subsetting:
    - create a submatrix (e.g. brain cells only)
- Focus:
    - reduce dataset size before advanced analysis
- Preparation steps:
    - get obs (cell metadata)
    - get var (gene metadata)
    - ensure data is “analysis-ready”
- 2 analysis pathways:
- Path A: Scanpy (classical pipeline)
    - Use Scanpy for:
        - clustering
        - classical dimensionality reduction
    - Strategy:
        - initially ignore annotations
        - focus on intrinsic data structure
    - Questions:
        - can cells be grouped?
        - do clusters correspond to known cell types?
    -   Additional steps:
        - batch correction
        - data batching
- Path B: TorchDR (deep learning / scalable DR)
    - Use TorchDR for:
        - dimensionality reduction
        - scalability experiments
    - Motivation:
        - memory limitations on GPU
        - large datasets cannot always be processed at once
    - Core question:
        - how do we represent ~20k features in 2D?
    - Methods:
        - PCA
        - incremental PCA
    - Characteristics:
        - PCA transforms data into fewer dimensions
        - fewer x/y parameters
    - Limitations:
        - memory problems on GPU
    - Comparison goals:
        - speed
        - incremental processing
        - scalability
- Dimensionality reduction concepts
    - PCA:
        - Principal Component Analysis
        - linear transformation
        - reduces dimensionality
    - Incremental PCA:
        - processes data in batches
        - suitable for large datasets
    - Dimensionality reduction is:
        - a prerequisite for clustering
        - a visualization tool
        - not the final biological result
- Subsetting and sampling strategy:
    - TorchDR applied to:
        - a subset of the dataset (e.g. 20k cells)
    - Tasks:
        - describe what the data subset contains
        - learn how to select the subset
    - Considerations:
        - representativeness
        - computational limits
- Infrastructure and env
    - Tools and platforms:
        - local machine
        - Google Colab
        - JupyterHub
        - KI Werkstatt (HTW)
    - Tasks:
        - install required libraries
        - check Scanpy tutorials
    - Environment readiness:
        - prepare for GPU usage
        - handle memory constraints
- Visualization & outputs
    - Expected outputs:
        - dimensionality reduction plots
        - clustering visualizations
    - Visual goals:
        - clear 2D representations
        - compare methods
    - Notes:
        - clusters may or may not align with annotations
        - visualization is exploratory, not confirmatory
- Compare:
    - Scanpy-based DR
    - TorchDR-based DR
- Metrics:
    - speed
    - scalability
    - feasibility
- Testing:
    - fast vs incremental methods
    - performance under larger data sizes
- Overall Workflow:
    - Load data via CELLxGENE API
    - Explore metadata and basic statistics
    - Filter cells and genes
    - Subset data for feasibility
    - Run classical Scanpy pipeline
    - Run TorchDR / incremental PCA
    - Visualize and compare embeddings
    - Reflect on scalability and limitations

## Meeting – 27-10-2025

**Focus of meeting**
- Preprocessing

**What I presented**
- load_data.py and explained the switch from native Windows to WSL

**Supervisor feedback**
- Switch to google collab or KI wersktatt
- 1.Preprocessing goals
    - Primary objective: filter and reduce the dataset to a form suitable for downstream dimensionality reduction
    - Emphasis on:
        - computational feasibility (memory, RAM)
        - meaningful signal extraction
        - avoiding uninformative features
- 2.Gene-level filtering
    - Perform filtering at the gene level
    - Restrict analysis to protein-coding genes
        - Motivation:
            - reduce feature space
            - improve interpretability
            - limit noise
    - Typical scale:
        - ~50,000 genes initially
        - reduced to protein-coding subset

- 3.Cell-level focus and biological context
    - Focus on human brain cells
    - Disease context includes:
        - Alzheimer’s disease
        - Epilepsy
    - Interest in:
        - specific cell subtypes (e.g. neuronal cells)
    - Subselection of cells based on:
        - tissue type
        - disease relevance
        - analysis goals

- 4.Highly Variable Genes (HVGs)
    - Further reduce the gene set by selecting highly variable genes
    - Rationale:
        - retain genes that vary meaningfully across cells
        - remove uninformative genes that show little change
    - Use Scanpy’s HVG selection:
        - based on mean–variance (dispersion) relationships
    - Outcome:
        - stabilize downstream PCA
        - reduce noise

- 5.Read counts and bioinformatics considerations
    - Read counts are a central measure in scRNA-seq data
    - Bioinformatics preprocessing maps:
        - raw read counts
        - gene expression profiles
    - Need to understand:
        - how read depth influences variance
        - how filtering affects biological signal

- 6.Dimensionality reduction through feature reduction
    - Reduce the number of columns (genes) in the data matrix
    - Example:
        - from ~100,000 × 50,000 to ~100,000 × 5,000
    - This step is necessary to:
        - enable PCA
        - reduce memory usage
        - improve runtime
- 7.Scanpy-based preprocessing strategy
    - Study Scanpy preprocessing tutorials in detail
    - Use Scanpy to:
        - identify highly variable genes
        - normalize data
        - prepare data for PCA
    - Ensure correct preprocessing order before dimensionality reduction
- 8.Subsampling for computational feasibility
    - Use random subsets of cells when needed
    - Motivation:
        - limited computing power
        - large-scale datasets
    - Subsampling is used as:
        - a practical compromise
        - a testbed for methods
- 9.Exploratory plots and diagnostics
    - Generate plots to validate preprocessing choices:
        - histograms of gene expression
        - distribution of highly variable genes
        - variance-related diagnostics
    - These plots serve as: 
        - sanity checks
        - evidence of understanding preprocessing effects

- 10.Zwischenbericht (intermediate report)
    - The intermediate report should:
        - illustrate preprocessing decisions
        - demonstrate understanding of the data
        - include relevant plots (e.g. histograms)
    - Target length:
        - approximately 1 page
    - Focus:
        - explanation, not results
        - clarity over completeness

## Meeting – 06-11-2025

**Focus of meeting**
- Feedback on current progress

**What I presented**
- mouse data download code
- project set up locally
- attempts at migrating to google collab

**Supervisor feedback**
- Overall assessment
    - Working with mouse data was explicitly discouraged:
        - The project focus must remain on human scRNA-seq data.
        - Mouse datasets may only be used for isolated technical tests, not for core analysis.
    - Manual dataset downloads (e.g. downloading .h5ad files via the browser) were discouraged:
        - The workflow should remain API-driven and reproducible.
    - Loading datasets by manually copying dataset IDs into code was considered suboptimal:
        - While technically functional, this approach was described as “half-automated”.
        - Preference was given to workflows that programmatically query and select datasets.
- Feedback on data aquisition strategy
    - Google Colab was suggested as a possible environment for scalability testing, but:
        - The observed memory crashes were acknowledged.
        - It was emphasized that infrastructure limitations should be explicitly documented, not worked around silently.
    - The importance of GPU-aware thinking was stressed:
        - Even if GPU resources are not currently available, the pipeline should be designed with GPU execution in mind.
        - TorchDR was highlighted as relevant primarily in the context of batching and GPU-based scalability.
    - KI-Werkstatt was suggested as a more suitable environment for future large-scale experiments.
- Feedback on infrastructure and scalability
- Expectations for next iteration

## Meeting – 13-11-2025

**Focus of meeting**
- show progress
- ask about the prefered structure for the zwischenbericht

**What I presented**
- the subset of 3k cells and subsetting criteria 
- the local, subset-based workflow
- preliminary dataset inspection

**Supervisor feedback**
- Initial 3k subset could be used initially.
- Preprocessing and variance diagnostics are th next step: scanpy tutorials
- The structure of the zwischenbericht should be as follow: abstract + introduction + env + methods and handling progress + exploratory stat analysis and visualization and conclusion. This will not be the final bericht s structure as this does not cover dr yet which is the meat of the problem. This covers current progress; working API access, a sampled dataset, summary tables, one figure and some first statistics. Must include a research question, a clear decription of methods and a concrete interpertation of final results for this stage.

## Meeting – 27-11-2025

**Focus of meeting**
- zwischenbericht feedback
- possible final bericht structure

**Supervisor feedback**
- 1.The intro has to be longer: 
    - better contextualize contextualize recent technological advances in single-cell sequencing. In particular, the development of high-throughput sequencing technologies has enabled the generation of extremely large scRNA-seq datasets.
    - Special emphasis should be placed on the role of dataset size and metadata. The availability of rich metadata (e.g. cell type annotations, disease status, donor information) introduces both opportunities for analysis and challenges related to data handling, scalability, and interpretation.
    - It was suggested to explicitly discuss why large-scale datasets and extensive metadata constitute a computational challenge, particularly in the context of dimensionality reduction and visualization.
    - Additionally, the distribution of cell types and biological subtypes (e.g. healthy vs diseased samples) should be examined to provide biological context and to assess potential imbalance in the data.
    - Finally, while modern data formats and APIs (e.g. h5ad files and CELLxGENE Census APIs) significantly simplify access to large datasets compared to a few years ago, it is still important to acknowledge the historical difficulty of working with scRNA-seq data and to position current tools as enabling technologies rather than trivial solutions.
- 2.Work env: 
    - check out KI werkstatt(prof recomendation) / Kunkel / Google collab / JupyterHub
    - mention software env and not only hardware
- 3.exploration of the metadata:
    - subset/ distribution/ diseases ...
    - Exploratory analysis should focus on inspecting available metadata and mapping selected metadata variables onto low-dimensional embeddings. In particular, metadata such as disease status and cell subtype annotations should be visualized on the dimensionality reduction results to assess whether biologically meaningful structure emerges.
    - After normalization and preprocessing, embeddings are expected to show partial separation of cells based on metadata (e.g. healthy vs diseased samples), which can be visualized using color-coded overlays. These visualizations serve as qualitative sanity checks rather than definitive biological conclusions.
    - In addition, exploratory plots should summarize the distribution of cell subtypes and the proportion of healthy versus diseased cells to identify potential class imbalance that may influence downstream analysis.
    - Show an informative plot about the metadata
    - Include an informative plot about each of the methods you used
    - DR should hold the biggest space
    - refer to you github repo for more details but the final bericht should be heavily focused on DR
    - Decribe pipelline and methods in your readme
- DR: 
    - Exploratory analysis of metadata and its visualization on low-dimensional embeddings is treated as an initial sanity-check step in the analysis pipeline. While useful for understanding dataset composition, this step is not the main focus of the study.
    - The core objective is the systematic investigation of dimensionality reduction methods for scRNA-seq data. Following standard Scanpy tutorials, PCA is used as the primary linear dimensionality reduction step and as the basis for downstream embedding methods.
    - The central analysis compares traditional dimensionality reduction approaches (PCA, UMAP, t-SNE) with Torch-based methods, with particular emphasis on scalability. As dataset size increases, both memory and runtime constraints may limit the applicability of classical methods, motivating the use of TorchDR and batch-based or GPU-accelerated approaches.
    - Where available, precomputed embeddings and annotations (e.g. UMAP or t-SNE coordinates provided via the CELLxGENE API) are used for reference and comparison, avoiding unnecessary recomputation. TorchDR embeddings are then generated on the same underlying data representation to enable a fair methodological comparison. 
    - The central figures of the analysis consist of PCA, UMAP, and t-SNE embeddings obtained through classical methods, alongside corresponding TorchDR embeddings, allowing direct visual and qualitative comparison under increasing data scale.
- For the plots in the zwischenbericht:
    - account for variance 
    - The current variance-related plots require revision, as the relationship between mean expression and dispersion is not correctly represented. A clearer understanding of how these diagnostic plots should appear is necessary in order to identify preprocessing issues.
    - It was recommended to apply more stringent filtering, both at the cell and gene level, to reduce noise and improve the quality of downstream visualizations. In particular, lowly expressed genes should be removed, as they contribute limited biological signal while inflating variance estimates.
    - Gene filtering thresholds should be explicitly examined and justified (e.g. retaining genes expressed in at least a fixed percentage of cells, such as 20%). This analysis should clarify what information is retained or lost through filtering and how these choices affect variance structure and subsequent dimensionality reduction results.

**Key decisions / clarifications**
- rewrite the intro according to new feedback
- redirect the attention to the DR when writing the bericht

## Meeting – 04-12-2025

**Focus of meeting**
- zwischenbericht feedback
- missing graphs and sections

**Supervisor feedback**
- In the intro: the final report must showcase an understandung of batches and their relevance (batch correction)
- In the intro: it is important to describe why single cell sequencing, cellxgene dataset and the metadata are relevant. 
- For reproducibility purposes create a Github repo with the code and progress log.
- Visualize how different dimensionality reduction methods transform the same scRNA-seq data (PCA, UMAP, t-SNE, TorchDR). 
    - Same input data
    - DIfferent DR methods
    - compare the resulting embeddings
- Visualize different DR methods -> what kind of influence does this method has on the data? 
    - each DR method imposes assumptions on the data and changes what structures becomes visible. We will be showing method induced bias.
    - compare by visualizing their embeddings on the same dataset, and analyze how each method influences the apparent structure and interpretation of the data
- Restrict the analysis to protein-coding genes to reduce noise and focus on genes with well-characterized functional roles.
- Examine the distribution of healthy and diseased samples within the dataset to assess class balance.

**Key decisions / clarifications**
- Review articles from cellxgene that have researched this topic for the intro.
- Create a Github repo with the code.
- PCA plot
- UMAP plot 
- torchDR PCA plot
- Answer the question: why do these plots look different? What does each method emphasize or distort? What conclusion are method-dependent? How much do we see is biology, and how much is method?

## Meeting – 11-12-2025

**Focus of meeting**
- Dimensionality Reduction, PCA, and QC Interpretation

**What I presented**
- Implemented a complete preprocessing and initial dimensionality reduction pipeline for a subsampled scRNA-seq dataset using Scanpy. 
- The focus was on correct quality control, normalization, highly variable gene (HVG) selection, and exploratory dimensionality reduction (PCA, UMAP, t-SNE). 
- At this stage, the workflow closely follows official Scanpy preprocessing tutorials to ensure methodological correctness, while interpretative confidence is still being developed.

**Supervisor feedback**
- 1.PCA fundamentals and interpretation
    - PCA projects high-dimensional data onto orthogonal axes (principal components), typically visualized as x- and y-axes in 2D plots.
    - PCA axes do not necessarily start at zero, and this is expected behavior.
    - PCA requires proper scaling and normalization; otherwise, misleading structures may appear.
    - The first principal component (PC1) explains the largest proportion of variance in the data.
    - Each subsequent component explains a decreasing fraction of the remaining variance.
- 2.Explained variance and component selection
    - The proportion of explained variance per principal component must be explicitly examined.
    - A cumulative explained variance plot should be used to determine:
        - how many principal components are required
        - what percentage of variance should be retained
    - Common practice:
        - inspect the first ~30–50 components
        - select the number of PCs that explain a sufficiently large proportion of variance (e.g. 70–90%), depending on the downstream task
    - The number of PCs retained influences:
        - clustering quality
        - neighborhood graph construction
        - downstream dimensionality reduction (UMAP, t-SNE, TorchDR)
- 3.Scaling considerations and dataset size
    - Scaling behavior is not primarily determined by the number of cells alone.
    - Computational complexity depends on the ratio between number of samples (cells) and number of features (genes).
    - Increasing the number of cells does not necessarily increase the number of features.
    - Highly Variable Gene (HVG) selection is therefore critical:
        - it stabilizes PCA
        - it keeps the feature space fixed as the dataset grows
    - This distinction is essential when reasoning about scalability.
- 4.PCA computation and TorchDR integration
    - PCA must be computed correctly before applying TorchDR.
    - TorchDR can be used:
        - on top of PCA embeddings
        - with incremental PCA for large datasets
    - Incremental PCA is especially relevant when:
        - memory limits are reached
        - GPU resources are constrained
    - TorchDR is motivated by:
        - scalability concerns
        - batching
        - GPU acceleration
- 5.Cell filtering and zero-expression issues
    - Cells with zero or near-zero gene expression should be removed prior to PCA.
    - Such cells can:
        - distort variance estimates
        - appear as artificial outliers in embeddings
    - Filtering decisions must be justified and documented.
- 6.Gene variability and HVG selection
    - Gene variability should be assessed using dispersion-based measures.
    - Highly Variable Genes (HVGs):
        - show meaningful variation across cells
        - often define biological cell types
    - HVG selection plots must be interpreted carefully:
        - the goal is not aesthetics but variance structure
    - Improper HVG selection leads to unstable PCA and misleading embeddings.
- 7.Outliers, dead cells, and biological vs technical interpretation
    - PCA and downstream embeddings may reveal apparent outliers.
    - These outliers can correspond to:
        - dying or stressed cells
        - technical artifacts
        - biologically interesting rare populations
    - Important clarification:
        - healthy donors can still contain dying cells
        - this is normal and expected
    - The presence of outliers alone does not justify removal without further investigation.
- 8.Mitochondrial gene content and QC thresholds
    - QC plots for mitochondrial gene expression are critical.
    - Cells with very high mitochondrial percentages (e.g. >50%) are often:
        - dying
        - stressed
    - However:
        - rigid thresholds should be applied cautiously
        - biological context matters
    - The decision to remove such cells must be explicitly reasoned, not automatic.
- 9.Total counts and doublets
    - Cells with unusually high total read counts (e.g. 80,000 vs typical 10,000) require inspection.
    - High read counts may indicate:
        - doublets (two cells captured together)
        - technical artifacts during droplet-based sequencing
    - Doublets arise when:
        - more than one cell is encapsulated in a single droplet
    - These artifacts manifest as:
        - dense “bubbles” or isolated points in embeddings
        - gray or unannotated outliers 
- 10.Droplet sequencing and experimental artifacts
    - Droplet-based scRNA-seq operates at the nanogram scale.
    - Sample preparation steps (e.g. buffer solutions, RNA extraction) are error-prone.
    - Errors introduced before or during sequencing propagate into:
        - gene expression measurements
        - downstream analyses
    - These technological artifacts should be explicitly acknowledged.
- 11.Importance for introduction and discussion
    - The impact of sequencing technology and preprocessing artifacts should be highlighted in the introduction.
    - Emphasize:
        - how fragile the measurement process is
        - why careful QC and preprocessing are essential
    - These considerations justify:
        - cautious interpretation of embeddings
        - focus on methodology rather than biological claims
- 12.Visualization and storytelling
    - PCA and DR plots are not just technical outputs; they tell a story.
    - For each plot, you should be able to explain:
        - what structure is visible
        - what could be biological
        - what could be technical
    - Gray or uncolored points should be questioned:
        - why do they appear?
        - what metadata do they lack?
    - These discussions sit at the intersection of:
        - biology
        - data science
        - experimental technology

## Meeting – 18-12-2025

**Focus of meeting**
- TorchDR integration and PCA handling
- Normalization & feature names
- Final report structure 

**What I presented**
- Updated preprocessing pipeline
- Correlations between principal components and QC metrics (mitochondrial percentage, gene counts).
- Neighborhood graphs UMAP and t-SNE embeddings for exploratory analysis.
- Initial TorchDR workflow
- A draft of the structure for the final report.

**Supervisor feedback**
- In dr.py, PCA was applied on top of an already computed PCA (X_pca).
- This resulted in a PCA-on-PCA workflow, which is conceptually wrong.
- TorchDR (including IncrementalPCA) must operate on the log-normalized gene expression matrix, not on PCA embeddings.
- Correct principle
    - Preprocessing step:
        - Normalization
        - Log-transformation
        - HVG selection
        - scaling
    - Dimensionality reduction step:
        - Classical PCA (Scanpy)
        - Incremental PCA (TorchDR)
    - These are alternative DR paths, not sequential ones.
- Decoupled preprocessing and DR:
    - preprocessing.py ends with a clean, processed adata
    - dr.py is responsible only for DR methods
    - TorchDR now receives:
        - adata.X (log-normalized, optionally HVG-filtered)
        - not adata.obsm["X_pca"]
- Normalization & feature names
    - Incremental PCA does not care about feature (gene) names, only about the numeric matrix.
    - Duplicate gene symbols are not a mathematical issue for PCA, but they are a data-structure issue in AnnData.
    - feature_name was used both as:
        - adata.var.index
        - and a column in adata.var
        - With conflicting values.
    - resolution was adata.var.rename which resolved the AnnData contraints
    - TorchDR itself was not the root cause, but exposed the inconsistency later in the pipeline.
    - Why did this conflict surface during DR and not earlier?
        - Code needs minor revision to ensure:
        - var_names are finalized before downstream steps.
- Scope and scaling: 
    - Do not scale beyond ~10k cells at this stage.
    - Focus is:
        - Conceptual correctness
        - Comparison between:
            - Classical PCA (Scanpy)
            - Incremental PCA (TorchDR)
        - Large-scale benchmarking can be discussed conceptually rather than executed fully.
- Must add: track memory usage for most steps in the piepline.
- env.yml must be cleaned up
- Final Report structure:
    - Structure is sound but the focus should be:
        - Methodological reasoning
            - Why each preprocessing and dimensionality reduction step exists
            - What assumptions each method makes
        - Pipeline decisions
            - Why PCA is used
            - Why incremental PCA is introduced
            - Why TorchDR is relevant for scalability
        - Limitations and trade-offs
            - Memory vs accuracy
            - Batch processing vs full-matrix methods
            - CPU vs GPU execution
    - Introduction and background can be merged
    - Biological interpretation should have a supportive role: 
        - Used only to validate or illustrate methodological choices
        - Not used to make novel biological claims
        - Cell types, disease labels, and QC metrics serve as sanity checks, not conclusions
    - Introduction: Do not oversell biological discoveries. Focus on:
        - Explosion of single-cell sequencing technologies
        - Resulting data volume and metadata complexity
        - Why dimensionality reduction is a computational necessity, not a visualization trick
        - Keep in intro:
            - Sequencing technologies → massive datasets
            - CELLxGENE Census → standardized access to large-scale data
            - Core problem: data no longer fits into memory → algorithms must adapt
            - Data availability (CELLxGENE)
            - Data standards (AnnData / .h5ad)
            - Computational challenge (scale, memory, batch processing)
    - github repo: env.yml must be changed because it s a dump rn
    - Stick to one dataset
        - Do not try to impress with size
        - TorchDR should be dataset-agnostic
        - Scalability argument is theoretical + empirical, not look how big
    - This section is good but needs clear justification per step:
        - QC thresholds (why 200 genes, why 5% mt)
        - Normalization (why total-count normalization)
        - Log transform (distributional reasons)
        - HVG selection (why dimensionality reduction depends on variance)
    - All of the plots for hvgs and qc metrics are more introductory, the main focus is pca and incremental pca, and the main  plots is the comparasion between them.
        - UMAP is not a result
        - UMAP is context / motivation
        - PCA is the real workhorse
    - Show:
        - Incremental PCA
        - Batch-wise computation
        - GPU memory constraints
        - Same mathematical goal, different computational strategy
    - In discussion show: 
        - PCA is not about plots
        - PCA is a foundational preprocessing step for:
        - kNN graphs
        - clustering
        - downstream embeddings
- Presentation suggestion: 5 slides:
    - PCA is not about plots
    - PCA is a foundational preprocessing step for:
    - kNN graphs
    - clustering
    - downstream embeddings
    - extra slides can exist as a backup or answer for further questions.
