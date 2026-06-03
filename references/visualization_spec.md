# Visualization specification

## Purpose

This module defines required visualization outputs for scVelo + CellRank2 fate mapping. The goal is to support interpretation, QC, and publication-style figure assembly.

## Visualization levels

Every analysis must output three levels of figures:

1. QC figures
2. Inference figures
3. Biological interpretation figures

Velocity stream alone is insufficient. It must be accompanied by latent time, velocity confidence, phase portraits, CellRank terminal states, fate probabilities, and driver gene trends.

## Required figure groups

### 1. Input and structure figures

Purpose: show that the selected data structure is suitable for trajectory analysis.

Required:

- UMAP colored by cell type or cluster.
- UMAP colored by sample, condition, and timepoint if available.
- UMAP or violin plots for total counts, detected genes, percent mitochondrial genes.
- UMAP or violin plots for spliced and unspliced counts if available.

### 2. RNA velocity figures

Required:

- velocity stream plot
- velocity grid plot or velocity arrow plot
- latent time UMAP
- velocity confidence UMAP
- velocity length UMAP

Recommended commands:

```python
scv.pl.velocity_embedding_stream(adata, basis="umap", color=["celltype", "condition"])
scv.pl.velocity_embedding_grid(adata, basis="umap", color="celltype")
scv.pl.velocity_embedding(adata, basis="umap", arrow_length=3, arrow_size=2)
scv.pl.scatter(adata, color=["velocity_length", "velocity_confidence", "latent_time"])
```

### 3. Phase portrait figures

Required:

- phase portraits for top velocity genes
- phase portraits for top CellRank driver candidates

Purpose:

- Validate whether spliced/unspliced kinetics support the inferred direction.
- Identify genes with poor kinetic fit.

Recommended command:

```python
scv.pl.velocity(adata, var_names=top_velocity_genes, basis="umap")
```

### 4. CellRank fate figures

Required:

- terminal states on UMAP
- initial states on UMAP
- macrostates on UMAP if available
- fate probability UMAP for each terminal lineage
- coarse transition matrix
- fate probability distribution across cell types, clusters, conditions, or timepoints

Recommended commands:

```python
g.plot_macrostates(which="terminal", basis="umap")
g.plot_macrostates(which="initial", basis="umap")
g.plot_fate_probabilities(basis="umap")
g.plot_coarse_T()
```

### 5. Driver gene figures

Required:

- lineage driver ranking table
- top driver dotplot or matrixplot
- top driver heatmap
- gene trends for selected driver candidates

Recommended commands:

```python
drivers = g.compute_lineage_drivers()
cr.pl.gene_trends(adata, model=model, genes=top_genes, time_key="latent_time", data_key="Ms")
```

## Publication-style main figure plan

Suggested layout:

```text
Figure X. scVelo + CellRank2 reveals lineage fate structure and candidate fate-associated drivers

A. UMAP of selected lineage subset colored by cell type / condition / timepoint
B. Velocity stream on UMAP
C. Latent time or velocity pseudotime on UMAP
D. CellRank terminal states and fate probabilities
E. Fate probability distribution across clusters / timepoints / conditions
F. Top lineage driver gene heatmap or dotplot
G. Gene trends of selected candidate drivers along latent time
H. Summary model of inferred transition and candidate regulatory programs
```

## Extended data figure plan

```text
Extended Data Fig. X

A. spliced / unspliced count QC
B. velocity confidence
C. phase portraits of top velocity genes
D. kernel-weight sensitivity: VelocityKernel vs ConnectivityKernel
E. GPCCA n_states sensitivity
F. terminal state robustness across parameter settings
G. driver gene overlap with DEG / markers / pathway genes
H. low-confidence or excluded lineages
```

## Caption rules

Each figure caption must state:

- what the figure shows
- what it supports
- what it does not prove
- whether the result is high-confidence or exploratory

Avoid causal claims in captions unless supported by perturbation or orthogonal validation.
