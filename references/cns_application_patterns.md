# CNS-style application patterns for fate mapping

## Principle

High-impact papers rarely rely on a single trajectory plot. They build a layered evidence chain: dynamic process, directionality, fate bias, driver candidates, and validation.

## Pattern 1. Dynamic process first

Before running CellRank2, define the biological process:

- differentiation
- activation
- exhaustion
- reprogramming
- regeneration
- inflammation-to-resolution transition
- stress hematopoiesis

If no plausible dynamic process exists, label the analysis exploratory.

## Pattern 2. Directionality needs external anchors

Use at least two anchors:

- known markers
- sampling timepoint
- lineage tracing
- perturbation
- spatial position
- histology
- flow cytometry
- public dataset

## Pattern 3. Fate probability is a main result

For branching processes, fate probability is usually more interpretable than pseudotime alone. Show:

- fate probability UMAP
- fate probability by cluster
- fate probability by condition/timepoint
- branch-biased intermediate states

## Pattern 4. Driver genes are candidates

Use layered support:

- CellRank correlation with fate probability
- expression trend
- branch specificity
- pathway enrichment
- TF/regulon evidence
- known literature
- perturbation if available

## Pattern 5. Main figure and extended data separation

Main figure:

- biological story
- directionality
- fate probability
- driver candidates

Extended data:

- QC
- phase portraits
- parameter sensitivity
- excluded low-confidence lineages

## Pattern 6. Appropriate wording

Use:

- `supports`
- `suggests`
- `is consistent with`
- `candidate driver`
- `fate-associated gene`

Avoid:

- `proves`
- `causes`
- `determines`
- `master regulator`
- `causal driver`
