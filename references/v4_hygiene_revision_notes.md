# v4 hygiene and release-candidate revision notes

This revision addresses the v3 re-audit focus on packaging hygiene, method-label accuracy, and redundant artifacts.

## Changes

- Removed all `__pycache__/` directories and `*.pyc` bytecode files before packaging.
- Updated `SKILL.md` so driver outputs are described as fate-probability-associated candidate genes. Native CellRank `compute_lineage_drivers()` and `correlation_fallback` are now distinguished explicitly.
- Updated `driver_gene_summary.json` generation so `interpretation` is dynamic according to the observed driver-gene methods.
- Updated report language to distinguish native CellRank driver candidates from fallback expression/fate-probability correlation candidates.
- Expanded `agents/openai.yaml` with `display_name`, `short_description`, and `default_prompt`.
- Added cell-level fate-assignment stability fields to `n_states_sensitivity.csv` when dependencies support them: `assignment_match_fraction` and `assignment_ami_vs_main`.
- Added optional fate assignment composition tables by cell type/cluster/condition/timepoint for each sensitivity setting.

## Validation in current environment

- Syntax validation: PASS using `py_compile` with temporary bytecode output.
- Lightweight smoke test: PASS.
- Minimal workflow test: SKIPPED in this container because `anndata`, `scanpy`, and `cellrank` are not installed.

## Remaining required validation before formal publication use

- Run `workflow_test_minimal.py` in an environment with `anndata`, `scanpy`, and `cellrank`.
- Run one small velocity-enabled dataset with `scvelo` and `cellrank` installed to confirm `VelocityKernel`, GPCCA export, native `compute_lineage_drivers()`, and figure generation.
- Treat the generated report as a computational audit; publication-grade biological interpretation still requires marker, timepoint, batch/sample, DEG/pathway, literature, and orthogonal validation review.
