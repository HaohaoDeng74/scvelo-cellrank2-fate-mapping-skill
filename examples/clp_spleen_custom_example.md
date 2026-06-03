# CLP spleen custom application example

## Context

For CLP spleen single-cell data, do not run global all-cell RNA velocity as a primary biological conclusion. The spleen contains unrelated immune and hematopoietic lineages. Global stream plots may reflect annotation geometry, batch, stress, or abundance differences rather than true fate transitions.

## Recommended lineage subsets

### 1. Erythroid / EMH axis

Potential question:

- Does CLP induce or remodel extramedullary hematopoiesis?
- Is there a maturation or stress erythropoiesis trajectory?

Markers to inspect:

- `Gata1`, `Klf1`, `Alas2`, `Hbb-bs`, `Hba-a1`, `Mki67`, `Top2a`

### 2. Monocyte–macrophage axis

Potential question:

- Are inflammatory monocytes transitioning toward macrophage-like or repair-like states?
- Does CLP alter fate bias among myeloid states?

Markers to inspect:

- `Ly6c2`, `S100a8`, `S100a9`, `Lcn2`, `Il1b`, `Ccl3`, `Ccl4`, `Mrc1`, `C1qa`, `C1qb`, `C1qc`

### 3. B cell to plasma cell axis

Potential question:

- Does CLP reshape antibody-producing cell differentiation?

Markers to inspect:

- `Ms4a1`, `Cd79a`, `Cd74`, `Mzb1`, `Jchain`, `Xbp1`, `Prdm1`, immunoglobulin genes

### 4. T cell activation / exhaustion / Treg axis

Potential question:

- Does CLP induce persistent T cell activation, suppression, exhaustion, or Treg-biased fate?

Markers to inspect:

- `Cd3d`, `Cd4`, `Cd8a`, `Il7r`, `Ccr7`, `Mki67`, `Pdcd1`, `Lag3`, `Havcr2`, `Foxp3`, `Il2ra`, `Ctla4`

## Suggested report interpretation

Use conservative wording:

- `CLP-associated myeloid cells showed a fate-probability gradient toward ...`
- `The inferred transition is consistent with inflammatory-to-phagocytic remodeling ...`
- `Candidate fate-associated genes included ...`

Avoid:

- `CLP causes monocytes to become ...` unless experimentally validated.
- `Gene X drives fate decision` unless perturbation or strong orthogonal evidence exists.
