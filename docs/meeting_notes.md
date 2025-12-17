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

## Meeting – 06-11-2025

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

## Meeting – 13-11-2025

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

## Meeting – 20-11-2025

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