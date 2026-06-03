# scVelo + CellRank2 fate mapping report

## 1. Analysis overview

Dataset: `{dataset_name}`  
Subset: `{subset_name}`  
Cells analyzed: `{n_cells}`  
Genes analyzed: `{n_genes}`  
Velocity mode: `{velocity_mode}`  
CellRank kernel: `{kernel_description}`  
Estimator: `{estimator}`  
Date: `{analysis_date}`

## 2. Input QC

Summarize whether the input data are appropriate for velocity and fate mapping.

- Cell annotation column:
- Sample/condition/timepoint metadata:
- Spliced/unspliced layers:
- UMAP/neighbors:
- Selected subset size:
- Major concerns:

## 3. Velocity results

### 3.1 Velocity direction

Describe velocity stream and whether it agrees with known biology.

### 3.2 Latent time

Describe latent time pattern and whether it agrees with markers or sampled timepoints.

### 3.3 Velocity QC

Describe velocity confidence, velocity length, and phase portraits.

## 4. CellRank fate inference

### 4.1 Initial states

Initial state(s):

### 4.2 Terminal states

Terminal state(s):

### 4.3 Fate probabilities

Summarize fate probability structure across clusters, conditions, and timepoints.

## 5. Lineage-associated driver candidates

For each lineage:

| Lineage | Top candidate genes | Evidence summary | Confidence |
|---|---|---|---|
| `{lineage}` | `{genes}` | `{evidence}` | `{grade}` |

## 6. Gene trends

Summarize whether top candidate genes show coherent dynamic expression trends.

## 7. Biological interpretation

### Lineage conclusion template

#### Inferred transition

Initial state:  
Terminal state:  
Intermediate state:  

#### Computational evidence

- Velocity direction:
- Latent time:
- Fate probability:
- Driver gene candidates:
- Gene trends:

#### Biological evidence

- Known markers:
- Sample timepoint consistency:
- DEG/pathway support:
- Literature support:
- Optional spatial/flow validation:

#### Confidence grade

Grade A/B/C/D

#### Interpretation

This result supports / suggests / weakly suggests ...

## 8. Limitations

List major limitations and low-confidence areas.

## 9. Recommended next analyses

- DEG overlap
- pathway enrichment
- TF/regulon inference
- spatial validation
- flow cytometry or perturbation validation
- public dataset validation
