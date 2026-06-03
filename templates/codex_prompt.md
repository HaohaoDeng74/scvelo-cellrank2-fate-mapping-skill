# Codex prompt for scVelo + CellRank2 fate mapping

You are working inside a skill package named `scvelo-cellrank2-fate-mapping`.

Follow these rules strictly:

1. Read `SKILL.md` first.
2. Read `references/qc_checklist.md` before input validation.
3. Read `references/failure_modes.md` if velocity confidence is low, terminal states look implausible, or batch/timepoint confounding appears.
4. Read `references/visualization_spec.md` before generating final figures.
5. Read `references/literature_guided_interpretation.md` before writing biological conclusions.

Execution order:

```bash
python scripts/00_validate_input.py --config config.yaml
python scripts/01_run_scvelo.py --config config.yaml
python scripts/02_velocity_qc.py --config config.yaml
python scripts/03_run_cellrank2.py --config config.yaml
python scripts/04_driver_gene_analysis.py --config config.yaml
python scripts/05_gene_trend_and_pathway.py --config config.yaml
python scripts/06_generate_report.py --config config.yaml
python scripts/smoke_test.py --skill-root .
```

Hard constraints:

- Do not run RNA velocity without `spliced` and `unspliced` layers.
- If velocity inputs are absent, use the non-velocity CellRank route and label it explicitly.
- Do not interpret CellRank driver genes as causal genes.
- Always write or update `run_manifest.json`.
- Never leave `{placeholder}` strings in the final report.
- Flag mitochondrial, ribosomal, hemoglobin, cell-cycle, stress and ambient RNA genes as high-risk driver candidates.
- Use lineage-specific subsets for complex immune tissues.
- Report confidence as A/B/C/D and state why.
