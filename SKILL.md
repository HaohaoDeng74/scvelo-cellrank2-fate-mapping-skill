---
name: scvelo-cellrank2-fate-mapping
description: RNA velocity and CellRank2 fate-mapping workflow for single-cell RNA-seq, including scVelo preprocessing, velocity QC, CellRank2 kernel selection, GPCCA fate probability inference, terminal/initial state review, lineage-associated driver candidate analysis, gene trend visualization, evidence grading, and publication-style reporting. Use when asked to run or design scVelo, RNA velocity, CellRank/CellRank2, fate probability, terminal state, initial state, latent time, or lineage driver analysis.
---

# scVelo + CellRank2 Fate Mapping Skill

## Purpose

Use this skill to run or design a reproducible scVelo + CellRank2 fate-mapping workflow from AnnData. The skill prioritizes: input validity, velocity assumptions, CellRank kernel choice, fate probability inference, lineage-associated driver candidate analysis, visualization QC, and conservative biological interpretation.

## Trigger conditions

Use this skill when the user asks for any of the following:

- RNA velocity or scVelo analysis
- CellRank or CellRank2 analysis
- cell fate, lineage trajectory, fate probability, initial state, terminal state, or latent time analysis
- lineage driver gene or fate-associated gene discovery
- publication-style visualization of dynamic single-cell fate inference
- infection/immunology cell-state transition analysis using scRNA-seq

## Bundled resources

Read additional files only when needed.

- `references/qc_checklist.md`: read when validating AnnData or deciding whether velocity is interpretable.
- `references/failure_modes.md`: read when velocity confidence is low, terminal states are implausible, batch effects dominate, or results conflict with known biology.
- `references/visualization_spec.md`: read when generating analysis figures, publication panels, or supplementary QC figures.
- `references/literature_guided_interpretation.md`: read when interpreting fate decisions and driver genes.
- `references/cns_application_patterns.md`: read when writing a paper-style narrative or figure plan.
- `references/api_compatibility.md`: read when CellRank/scVelo versions differ or API calls fail.
- `templates/report_template.md`: use for Markdown report generation.
- `templates/figure_plan_template.md`: use for main/extended figure planning.
- `templates/codex_prompt.md`: use when delegating the task to Codex/agent.

## Required inputs

Minimum input:

- AnnData `.h5ad` file.
- Cell or cluster annotation in `adata.obs`, e.g. `celltype`, `major_celltype`, or `leiden`.
- Low-dimensional embedding, preferably `adata.obsm["X_umap"]`, or permission to recompute PCA/neighbors/UMAP.

Preferred velocity input:

- `adata.layers["spliced"]`
- `adata.layers["unspliced"]`

If `spliced` and `unspliced` are absent, do not run RNA velocity. Use the non-velocity CellRank route and label the output as non-velocity fate mapping.

## Execution workflow

Use `config_template.yaml` as the starting point.

```bash
cp config_template.yaml config.yaml
python scripts/00_validate_input.py --config config.yaml
python scripts/01_run_scvelo.py --config config.yaml
python scripts/02_velocity_qc.py --config config.yaml
python scripts/03_run_cellrank2.py --config config.yaml
python scripts/04_driver_gene_analysis.py --config config.yaml
python scripts/05_gene_trend_and_pathway.py --config config.yaml
python scripts/06_generate_report.py --config config.yaml
python scripts/smoke_test.py --skill-root .
```

All scripts write to `run_manifest.json`. Each step must record status, warnings, errors, input files, output files, parameters, and software versions where available.

## Decision rules

### Velocity route

Use the velocity route only when both `spliced` and `unspliced` layers exist and have compatible dimensions. Report unspliced/spliced count distributions and the unspliced fraction. If unspliced fraction is extremely low, mark velocity confidence as high risk and downgrade biological conclusions.

### Non-velocity route

If spliced/unspliced layers are unavailable, run a non-velocity CellRank route. Prefer:

1. `PseudotimeKernel` when a valid pseudotime column is configured.
2. `RealTimeKernel` when a valid experimental time column is configured and ordered.
3. `ConnectivityKernel` when no directional prior is available.

Clearly report that fate direction is not RNA-velocity-supported.

### Lineage-specific analysis

For immune or infected tissues, do not interpret global all-cell velocity across unrelated lineages. Subset to coherent compartments first, then recompute neighbors if configured.

Recommended subsets:

- Erythroid / extramedullary hematopoiesis lineage
- Monocyte–macrophage lineage
- Neutrophil / granulocytic lineage
- B cell to plasmablast/plasma cell lineage
- T cell activation/exhaustion/Treg axis
- Dendritic cell activation/differentiation axis

## Core computational outputs

Required output directories:

```text
00_input_qc/
01_velocity/
02_velocity_qc/
03_cellrank_fate/
04_driver_genes/
05_gene_trends/
06_report/
07_publication_figures/
```

Required machine-readable outputs:

- `run_manifest.json`
- `00_input_qc/input_validation.json`
- `01_velocity/adata_scvelo.h5ad` when velocity route is available
- `03_cellrank_fate/adata_cellrank2.h5ad`
- `03_cellrank_fate/fate_probabilities.csv` when available
- `03_cellrank_fate/terminal_states.csv` or a warning if unavailable
- `03_cellrank_fate/n_states_sensitivity.csv` when configured
- `04_driver_genes/lineage_driver_candidates.csv`
- `05_gene_trends/gene_trends_summary.csv`
- `06_report/scvelo_cellrank2_report.md`
- `07_publication_figures/figure_plan.md`

## Biological interpretation rules

- Velocity stream alone is insufficient evidence.
- Terminal states must be checked against cell type labels, timepoints, markers, and sample composition.
- Lineage driver outputs are fate-probability-associated candidate genes, not causal drivers. Prefer native CellRank `compute_lineage_drivers()` when available; otherwise use the explicitly labeled `correlation_fallback` table.
- Flag mitochondrial, ribosomal, hemoglobin, cell-cycle, stress, and ambient RNA genes as high-risk candidates.
- Report confidence as A/B/C/D using computational and biological evidence. If velocity falls back from dynamical to stochastic, grade cannot exceed B. If velocity is absent, fate direction cannot exceed C unless supported by real-time, lineage tracing, perturbation, or strong biological priors.

## Evidence grading

The automatic report is a computational audit. Publication-grade biological interpretation requires marker, timepoint, batch/sample, DEG/pathway, literature, and orthogonal validation review.


Grade A: velocity direction, latent time/pseudotime, fate probabilities, terminal states, marker dynamics, sample timepoints, and literature/orthogonal evidence are concordant.

Grade B: computational evidence is mostly concordant, but experimental or literature support is incomplete.

Grade C: trajectory or fate probability is suggestive but directionality or driver interpretation is uncertain.

Grade D: velocity assumptions fail, cells are unrelated, phase portraits are poor, terminal states are unstable, or results contradict known biology.

## Failure handling

Never silently continue after a critical failure. Write the failure to `run_manifest.json`, emit a clear warning, and route to the safest available analysis branch. If a result is exploratory, state that explicitly in the report.
