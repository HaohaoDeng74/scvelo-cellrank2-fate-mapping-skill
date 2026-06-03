# Literature-guided interpretation

## Purpose

This document guides biological interpretation of scVelo + CellRank2 outputs using method papers, CNS-style application patterns, and best-practice constraints.

## Methodological foundation

### RNA velocity

La Manno et al. introduced RNA velocity as a way to infer future transcriptional states by distinguishing unspliced and spliced mRNA abundances. The method estimates the time derivative of gene expression and can predict future cell states over a short timescale.

Interpretation principle:

- RNA velocity should be interpreted as short-timescale transcriptional directionality, not as direct long-term lineage proof.
- The predicted direction must be checked against biology, markers, sampling design, and phase portraits.

### scVelo

scVelo generalizes RNA velocity using a likelihood-based dynamical model of splicing kinetics. It is particularly useful for transient cell states, development, and perturbation responses.

Interpretation principle:

- Prefer dynamical mode when feasible.
- Use stochastic mode as a fallback or exploratory analysis.
- Always inspect velocity confidence and phase portraits.

### CellRank and CellRank2

CellRank integrates directionality and cell-state similarity into a Markov-chain framework to infer initial states, terminal states, fate probabilities, and lineage-associated driver genes. CellRank2 extends this idea to multiview single-cell data.

Interpretation principle:

- Fate probabilities are probabilistic absorption tendencies toward terminal states, not direct lineage labels.
- Driver genes are ranked by association with fate probabilities, not causality.
- Kernel choice must match data modality and biological question.

## CNS-style application patterns

### Pattern 1. Directionality must be biologically anchored

Do not present velocity stream alone. Anchor it using at least two of the following:

- known lineage markers
- sampled timepoints
- developmental or perturbation stage
- clonal tracing or genetic barcoding
- spatial localization
- flow cytometry or histology
- public dataset validation

### Pattern 2. Fate probability is more informative than pseudotime alone

For branching systems, report:

- terminal states
- initial states
- fate probability distribution
- branch-biased intermediate populations
- branch-specific driver gene candidates
- gene trends along each fate

### Pattern 3. Driver genes require layered evidence

A high-confidence lineage-associated driver candidate should satisfy:

1. Strong CellRank driver score or correlation with fate probability.
2. Dynamic expression trend along latent time or lineage-specific pseudotime.
3. Cell-state or branch specificity.
4. DEG, marker, pathway, TF/regulon, or literature support.
5. Preferably orthogonal validation or perturbation evidence.

### Pattern 4. Trajectory inference is a hypothesis generator

Velocity and CellRank outputs support hypotheses about dynamic transitions. They do not prove causality without experimental evidence.

### Pattern 5. Subset before inference

For complex immune tissues, do not infer all-cell fate maps directly. Recommended order:

1. global annotation
2. lineage subset
3. recompute neighbors/UMAP if needed
4. velocity and CellRank within subset
5. integrate results back into global atlas

## Infection/immunology-specific notes

In CLP spleen, infection, sepsis, transplant, or inflammatory datasets, consider the following axis-specific logic:

### Erythroid / EMH axis

Check whether inferred direction agrees with erythroid maturation or stress erythropoiesis markers such as `Gata1`, `Klf1`, `Alas2`, `Hbb`, `Hba`, `Mki67`, and `Top2a`.

### Myeloid activation axis

Check whether inferred direction agrees with inflammatory, antimicrobial, phagocytic, or repair-associated programs such as `Ly6c2`, `S100a8`, `S100a9`, `Lcn2`, `Il1b`, `Ccl3`, `Ccl4`, `Mrc1`, `C1qa`, `C1qb`, and `C1qc`.

### B cell to plasma cell axis

Check whether inferred direction agrees with `Ms4a1`, `Cd79a`, `Cd74`, `Mzb1`, `Jchain`, `Xbp1`, `Prdm1`, and immunoglobulin genes.

### T cell activation/exhaustion/Treg axis

Check whether inferred direction agrees with `Cd3d`, `Cd4`, `Cd8a`, `Il7r`, `Ccr7`, `Mki67`, `Pdcd1`, `Lag3`, `Havcr2`, `Foxp3`, `Il2ra`, and `Ctla4`.

## Required interpretation template

For each lineage, write:

```markdown
### Inferred transition
Initial state: ...
Terminal state: ...
Intermediate state: ...

### Computational evidence
- Velocity direction:
- Latent time:
- CellRank fate probability:
- Driver gene candidates:
- Gene trends:

### Biological evidence
- Known markers:
- Sample timepoint consistency:
- DEG/pathway support:
- Literature support:
- Optional spatial/flow validation:

### Confidence grade
Grade A/B/C/D

### Interpretation
This result supports / suggests / weakly suggests ...
```
