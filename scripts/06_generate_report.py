#!/usr/bin/env python
"""Generate Markdown report and publication figure plan from machine-readable outputs."""

from __future__ import annotations

import argparse
import re
from datetime import date
from pathlib import Path

import pandas as pd

from utils import load_config, output_dir, read_json, update_manifest, write_json


def safe(value, default="not available"):
    if value is None or value == "":
        return default
    return value


def load_csv(path):
    path = Path(path)
    if path.exists():
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.DataFrame()
    return pd.DataFrame()


def build_report(config):
    out = output_dir(config)
    manifest = read_json(out / "run_manifest.json", default={})
    input_qc = read_json(out / "00_input_qc" / "input_validation.json", default={})
    velocity_status = read_json(out / "01_velocity" / "scvelo_status.json", default={})
    velocity_qc = read_json(out / "02_velocity_qc" / "velocity_qc_summary.json", default={})
    cellrank_status = read_json(out / "03_cellrank_fate" / "cellrank_status.json", default={})
    driver_summary = read_json(out / "04_driver_genes" / "driver_gene_summary.json", default={})
    trend_summary = read_json(out / "05_gene_trends" / "gene_trend_summary.json", default={})

    drivers = load_csv(out / "04_driver_genes" / "top_lineage_driver_candidates.csv")
    sens = load_csv(out / "03_cellrank_fate" / "n_states_sensitivity.csv")

    route = cellrank_status.get("route", "not available")
    directionality = cellrank_status.get("directionality_supported_by", "not available")
    velocity_mode = velocity_status.get("velocity_mode_used", "not computed")
    median_conf = velocity_qc.get("median_velocity_confidence")

    if directionality == "connectivity_only":
        max_grade = "D"
        max_grade_note = "Connectivity-only route has no intrinsic directionality; report as state-structure mapping unless external time, lineage tracing, or perturbation data supplies direction."
    elif route.startswith("non_velocity"):
        max_grade = "C"
        max_grade_note = "Non-velocity route: fate direction is exploratory unless supported by external time, lineage tracing, or perturbation data."
    elif velocity_mode == "stochastic":
        max_grade = "B"
        max_grade_note = "Velocity fell back to stochastic mode; conclusions should not exceed moderate confidence without orthogonal support."
    else:
        max_grade = "A"
        max_grade_note = "Velocity route completed; confidence still depends on QC, phase portraits, marker direction, and biological plausibility."

    def cap_grade(grade):
        order = ["A", "B", "C", "D"]
        return order[max(order.index(grade), order.index(max_grade))]

    top_rows = []
    if not drivers.empty and {"lineage", "gene"}.issubset(drivers.columns):
        for lineage, sub in drivers.groupby("lineage"):
            genes = ", ".join(map(str, sub.head(10)["gene"].tolist()))
            risks = int(sub.get("high_risk_gene", pd.Series(False, index=sub.index)).astype(bool).sum()) if "high_risk_gene" in sub.columns else 0
            methods = ", ".join(sorted(map(str, sub.get("method", pd.Series(["unknown"])).dropna().unique().tolist())))
            evidence = f"{len(sub)} candidates; {risks} technical-risk genes; method={methods}."
            raw_grade = "B" if input_qc.get("velocity_possible") and not risks else "C"
            top_rows.append(f"| {lineage} | {genes or 'not available'} | {evidence} | {cap_grade(raw_grade)} |")
    if not top_rows:
        top_rows.append("| not available | not available | driver module did not produce candidate genes | D |")
    warnings = []
    for _, step in (manifest.get("steps", {}) or {}).items():
        warnings.extend(step.get("warnings", []) or [])
    warnings_md = "\n".join([f"- {w}" for w in warnings]) if warnings else "- No warnings recorded."

    sensitivity_md = "No n_states sensitivity table available."
    if not sens.empty:
        sensitivity_md = sens.to_markdown(index=False)

    text = f"""# scVelo + CellRank2 fate mapping report

## 1. Analysis overview

Dataset: {safe(config.get('project', {}).get('name'))}  
Subset: {safe(config.get('subset', {}).get('subset_name'))}  
Analysis date: {date.today()}  
Cells analyzed: {safe(input_qc.get('n_cells'))}  
Genes analyzed: {safe(input_qc.get('n_genes'))}  
Velocity mode: {safe(velocity_mode)}  
CellRank route: {safe(route)}  
Directionality supported by: {safe(directionality)}  
CellRank kernel: {safe(cellrank_status.get('kernel_description'))}  
Estimator: {safe(config.get('cellrank', {}).get('estimator', 'GPCCA'))}

## 2. Input QC

Input status: {safe(input_qc.get('status'))}  
Velocity possible: {safe(input_qc.get('velocity_possible'))}  
Spliced layer: {safe(input_qc.get('has_spliced'))}  
Unspliced layer: {safe(input_qc.get('has_unspliced'))}  
Unspliced fraction: {safe(input_qc.get('unspliced_fraction'))}  
Subset cells: {safe((input_qc.get('subset') or {}).get('n_cells'))}  
UMAP available: {'X_umap' in (input_qc.get('obsm') or [])}

## 3. Velocity results

scVelo status: {safe(velocity_status.get('status'))}  
Number of velocity genes: {safe(velocity_status.get('n_velocity_genes'))}  
Median velocity confidence: {safe(velocity_qc.get('median_velocity_confidence'))}  
Median velocity length: {safe(velocity_qc.get('median_velocity_length'))}

Interpretation guardrail: {max_grade_note}

## 4. CellRank fate inference

CellRank status: {safe(cellrank_status.get('status'))}  
Route: {safe(route)}  
Directionality supported by: {safe(directionality)}  
Kernel description: {safe(cellrank_status.get('kernel_description'))}  
Fate probabilities exported: {safe(cellrank_status.get('fate_probabilities_exported'))}  
Terminal states exported: {safe(cellrank_status.get('terminal_states_exported'))}  
Initial states exported: {safe(cellrank_status.get('initial_states_exported'))}  
Coarse transition matrix exported: {safe(cellrank_status.get('coarse_transition_matrix_exported'))}

### n_states sensitivity

{sensitivity_md}

## 5. Lineage-associated driver candidates

Driver summary: {safe(driver_summary.get('interpretation'))}

| Lineage | Top candidate genes | Evidence summary | Confidence |
|---|---|---|---|
{chr(10).join(top_rows)}

## 6. Gene trends and pathway support

Trend rows: {safe(trend_summary.get('n_trend_rows'))}  
Pathway enrichment rows: {safe(trend_summary.get('n_enrichment_rows'))}

## 7. Warnings and limitations

{warnings_md}

## 8. Interpretation template

For each lineage, report:

- inferred initial/intermediate/terminal states;
- whether velocity or non-velocity direction supports the transition;
- fate probability gradient across clusters, samples, and timepoints;
- top fate-associated candidate genes;
- whether candidate genes show coherent dynamic trends;
- whether candidates are plausible markers/regulators rather than technical genes;
- confidence grade A/B/C/D.

Evidence grade ceiling for this run: {max_grade}.  

Lineage driver outputs are fate-associated candidate genes and should not be described as causal drivers without perturbation, lineage tracing, or orthogonal validation. If the method is `cellrank_native_compute_lineage_drivers`, describe it as native CellRank fate-associated driver candidates. If the method is `correlation_fallback`, describe it as a fallback expression/fate-probability correlation candidate list rather than native CellRank driver analysis.
"""
    # Guard against template placeholders.
    text = re.sub(r"\{[^{}\n]+\}", "not available", text)
    return text


def build_figure_plan(config):
    out = output_dir(config)
    text = f"""# Publication figure plan

## Main figure

Figure X. scVelo + CellRank2 reveals lineage fate structure and candidate fate-associated drivers.

A. UMAP of selected lineage subset colored by cell type, condition and timepoint.  
B. Velocity stream on UMAP, or non-velocity kernel route diagram if spliced/unspliced are absent.  
C. Latent time / pseudotime / experimental time map.  
D. CellRank terminal and initial states.  
E. Fate probability maps shown separately for each lineage.  
F. Fate probability distribution across clusters, timepoints and conditions.  
G. Top lineage-associated candidate driver genes, with technical-risk genes flagged.  
H. Gene trend summaries for selected candidates.  
I. Model diagram with explicit confidence grade and validation needs.

## Extended data figure

A. Input spliced/unspliced QC and unspliced fraction.  
B. Velocity confidence and velocity length.  
C. Phase portraits for top velocity genes.  
D. Kernel selection and weight summary.  
E. GPCCA n_states sensitivity table.  
F. Driver gene risk-category annotation.  
G. Gene trend bin summaries.  
H. Failure modes and excluded lineages.

## Output locations

- Report: `{out / '06_report' / 'scvelo_cellrank2_report.md'}`
- Driver genes: `{out / '04_driver_genes' / 'lineage_driver_candidates.csv'}`
- Gene trends: `{out / '05_gene_trends' / 'gene_trends_summary.csv'}`
"""
    return text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    step = "06_generate_report"
    outdir = output_dir(config) / "06_report"
    figdir = output_dir(config) / "07_publication_figures"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    report = build_report(config)
    report_path = outdir / "scvelo_cellrank2_report.md"
    report_path.write_text(report, encoding="utf-8")

    fig_plan = build_figure_plan(config)
    fig_path = figdir / "figure_plan.md"
    fig_path.write_text(fig_plan, encoding="utf-8")

    placeholders = re.findall(r"\{[^{}\n]+\}", report)
    status = "PASS" if not placeholders else "WARN"
    summary = {"status": status, "unfilled_placeholders": placeholders}
    write_json(summary, outdir / "report_status.json")
    update_manifest(config, step, status, outputs=[str(report_path), str(fig_path), str(outdir / "report_status.json")], warnings=[f"Unfilled placeholders: {placeholders}"] if placeholders else [])


if __name__ == "__main__":
    main()
