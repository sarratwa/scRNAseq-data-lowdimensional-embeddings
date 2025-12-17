# Meeting Notes

This document summarizes meetings with the supervisor.
Each entry records feedback, decisions, and agreed action items.

---

## Meeting – 13-10-2025

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

---

## Meeting – 20-10-2025

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

## Meeting – 04-12-2025

**Focus of meeting**
- zwischenbericht feedback
- possible final bericht structure

**What I presented**
- 

**Supervisor feedback**
- 

**Key decisions / clarifications**
- 

**Action items**
- [ ] 
- [ ] 

## Meeting – 11-12-2025

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