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
- 

**What I presented**
- 

**Supervisor feedback**
- 

**Key decisions / clarifications**
- 

**Action items**
- [ ] 
- [ ] 

## Meeting – 18-12-2025

**Focus of meeting**
- calculating the PCA and torchDR

**What I presented**
- 

**Supervisor feedback**
- 

**Key decisions / clarifications**
- 

**Action items**
- [ ] 
- [ ] 