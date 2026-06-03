#!/usr/bin/env python
"""Generate basic gene trend summaries and optional pathway support tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from utils import fail_step, get_matrix, load_config, output_dir, update_manifest, write_json


def to_dense_2d(x):
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


def load_driver_table(config, outdir):
    path = output_dir(config) / "04_driver_genes" / "top_lineage_driver_candidates.csv"
    if not path.exists():
        path = output_dir(config) / "04_driver_genes" / "lineage_driver_candidates.csv"
    if not path.exists():
        raise RuntimeError("Driver gene table not found. Run 04_driver_gene_analysis.py first.")
    return pd.read_csv(path)


def trend_summary(adata, drivers, config, warnings):
    trend_cfg = config.get("gene_trends", {})
    drv_cfg = config.get("driver_genes", {})
    time_key = trend_cfg.get("time_key", "latent_time")
    if time_key not in adata.obs.columns:
        if "latent_time" in adata.obs.columns:
            time_key = "latent_time"
        else:
            warnings.append("No latent_time/pseudotime column found; trend summary skipped.")
            return pd.DataFrame()

    layer = drv_cfg.get("expression_layer", None)
    X = to_dense_2d(get_matrix(adata, layer=layer, use_raw=drv_cfg.get("use_raw", False))).astype(float, copy=False)
    genes = np.asarray(adata.raw.var_names if drv_cfg.get("use_raw", False) and getattr(adata, "raw", None) is not None else adata.var_names)
    gene_index = {g: i for i, g in enumerate(map(str, genes))}

    t = adata.obs[time_key].to_numpy(dtype=float)
    valid = np.isfinite(t)
    n_bins = int(trend_cfg.get("n_bins", 20))
    if valid.sum() < max(20, n_bins):
        warnings.append("Too few finite time values for trend summary.")
        return pd.DataFrame()

    quantiles = np.linspace(0, 1, n_bins + 1)
    edges = np.unique(np.nanquantile(t[valid], quantiles))
    if len(edges) < 3:
        warnings.append("Time values have too few unique quantile bins.")
        return pd.DataFrame()
    bins = np.digitize(t, edges[1:-1], right=True)

    rows = []
    max_genes_per_lineage = int(config.get("driver_genes", {}).get("trend_genes_per_lineage", 8))
    if "passes_default_filter" in drivers.columns:
        drivers = drivers[drivers["passes_default_filter"].astype(bool)]
    drivers = drivers.sort_values(["lineage", "abs_correlation"], ascending=[True, False]).groupby("lineage").head(max_genes_per_lineage)

    for _, r in drivers.iterrows():
        gene = str(r["gene"])
        if gene not in gene_index:
            continue
        expr = X[:, gene_index[gene]]
        for b in sorted(set(bins[valid])):
            mask = valid & (bins == b)
            rows.append({
                "lineage": str(r["lineage"]),
                "gene": gene,
                "time_key": time_key,
                "bin": int(b),
                "bin_time_min": float(np.nanmin(t[mask])),
                "bin_time_max": float(np.nanmax(t[mask])),
                "mean_expression": float(np.nanmean(expr[mask])),
                "n_cells": int(mask.sum()),
            })
    return pd.DataFrame(rows)


def parse_gmt(path):
    gene_sets = {}
    if not path:
        return gene_sets
    path = Path(path)
    if not path.exists():
        return gene_sets
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            gene_sets[parts[0]] = set(parts[2:])
    return gene_sets


def enrichment(drivers, adata, config, warnings):
    gmt = config.get("gene_trends", {}).get("gene_sets_gmt")
    gene_sets = parse_gmt(gmt)
    if not gene_sets:
        warnings.append("No valid GMT gene sets supplied; pathway enrichment skipped.")
        return pd.DataFrame()

    try:
        from scipy.stats import fisher_exact
    except Exception:
        warnings.append("scipy unavailable; pathway enrichment skipped.")
        return pd.DataFrame()

    background = set(map(str, adata.var_names))
    rows = []
    selected = drivers[drivers.get("passes_default_filter", True).astype(bool)] if "passes_default_filter" in drivers.columns else drivers
    for lineage, sub in selected.groupby("lineage"):
        genes = set(map(str, sub["gene"])) & background
        for term, term_genes_all in gene_sets.items():
            term_genes = set(map(str, term_genes_all)) & background
            a = len(genes & term_genes)
            b = len(genes - term_genes)
            c = len(term_genes - genes)
            d = len(background - genes - term_genes)
            if a == 0:
                continue
            _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({
                "lineage": str(lineage),
                "term": term,
                "overlap": a,
                "driver_genes": ";".join(sorted(genes & term_genes)),
                "p_value": float(p),
                "background_size": len(background),
            })
    out = pd.DataFrame(rows)
    if len(out):
        from utils import bh_adjust
        out["q_value"] = bh_adjust(out["p_value"].to_numpy())
        out = out.sort_values(["lineage", "q_value", "p_value"])
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    step = "05_gene_trend_and_pathway"
    outdir = output_dir(config) / "05_gene_trends"
    outdir.mkdir(parents=True, exist_ok=True)
    warnings = []
    outputs = []
    status = "PASS"

    try:
        input_h5ad = args.input or str(output_dir(config) / "03_cellrank_fate" / "adata_cellrank2.h5ad")
        adata = ad.read_h5ad(input_h5ad)
        drivers = load_driver_table(config, outdir)

        trends = trend_summary(adata, drivers, config, warnings)
        trends.to_csv(outdir / "gene_trends_summary.csv", index=False)
        outputs.append(str(outdir / "gene_trends_summary.csv"))
        if trends.empty:
            status = "WARN"

        enr = enrichment(drivers, adata, config, warnings)
        enr.to_csv(outdir / "pathway_enrichment.csv", index=False)
        outputs.append(str(outdir / "pathway_enrichment.csv"))

        summary = {
            "status": status,
            "n_trend_rows": int(len(trends)),
            "n_enrichment_rows": int(len(enr)),
            "warnings": warnings,
        }
        write_json(summary, outdir / "gene_trend_summary.json")
        outputs.append(str(outdir / "gene_trend_summary.json"))

        update_manifest(config, step, status if not warnings else "WARN", inputs=[str(input_h5ad)], outputs=outputs, warnings=warnings, summary=summary)
    except Exception as exc:
        fail_step(config, step, exc, warnings=warnings)
        raise


if __name__ == "__main__":
    main()
