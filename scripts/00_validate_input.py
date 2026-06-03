#!/usr/bin/env python
"""Validate input AnnData for scVelo + CellRank2 fate mapping."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np

from utils import dense_sum, fail_step, layer_total, load_config, output_dir, update_manifest, write_json


def subset_if_needed(adata, config, report):
    subset_cfg = config.get("subset", {})
    if not subset_cfg.get("enabled", False):
        report["subset"] = {"enabled": False, "n_cells": int(adata.n_obs)}
        return adata
    obs_key = subset_cfg.get("obs_key")
    include_values = subset_cfg.get("include_values", [])
    if obs_key not in adata.obs.columns:
        report["errors"].append(f"Subset obs_key `{obs_key}` not found in adata.obs.")
        report["subset"] = {"enabled": True, "obs_key": obs_key, "n_cells": 0}
        return adata[:0].copy()
    mask = adata.obs[obs_key].isin(include_values)
    adata_sub = adata[mask].copy()
    n = int(adata_sub.n_obs)
    report["subset"] = {
        "enabled": True,
        "obs_key": obs_key,
        "include_values": list(map(str, include_values)),
        "n_cells": n,
    }
    if n < int(subset_cfg.get("min_cells_stop", 100)):
        report["errors"].append(f"Subset contains {n} cells, below min_cells_stop.")
    elif n < int(subset_cfg.get("min_cells_warning", 300)):
        report["warnings"].append(f"Subset contains {n} cells, below min_cells_warning.")
    return adata_sub


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    step = "00_validate_input"
    outdir = output_dir(config) / "00_input_qc"
    outdir.mkdir(parents=True, exist_ok=True)

    report = {
        "input_h5ad": str(config["project"]["input_h5ad"]),
        "exists": False,
        "status": "PASS",
        "errors": [],
        "warnings": [],
    }

    try:
        input_h5ad = Path(config["project"]["input_h5ad"])
        report["exists"] = input_h5ad.exists()
        if not input_h5ad.exists():
            report["errors"].append(f"Input file not found: {input_h5ad}")
            report["status"] = "FAIL"
            write_json(report, outdir / "input_validation.json")
            update_manifest(config, step, "FAIL", outputs=[str(outdir / "input_validation.json")], errors=report["errors"])
            raise SystemExit(1)

        adata = ad.read_h5ad(input_h5ad)
        report["n_cells_raw"] = int(adata.n_obs)
        report["n_genes_raw"] = int(adata.n_vars)
        report["obs_columns"] = list(map(str, adata.obs.columns))
        report["layers"] = list(map(str, adata.layers.keys()))
        report["obsm"] = list(map(str, adata.obsm.keys()))
        report["obsp"] = list(map(str, adata.obsp.keys()))

        # Basic metadata checks.
        meta = config.get("metadata", {})
        for key_name in ["celltype_key", "cluster_key", "condition_key", "timepoint_key", "sample_key", "batch_key"]:
            obs_key = meta.get(key_name)
            if obs_key and obs_key not in adata.obs.columns:
                report["warnings"].append(f"Configured {key_name} `{obs_key}` not found in adata.obs.")
            elif obs_key:
                report[f"{key_name}_n_categories"] = int(adata.obs[obs_key].astype(str).nunique())

        # Layer dimension checks.
        for layer in ["spliced", "unspliced"]:
            if layer in adata.layers and adata.layers[layer].shape != adata.X.shape:
                report["errors"].append(f"Layer `{layer}` shape {adata.layers[layer].shape} does not match X shape {adata.X.shape}.")

        has_spliced = "spliced" in adata.layers
        has_unspliced = "unspliced" in adata.layers
        report["has_spliced"] = bool(has_spliced)
        report["has_unspliced"] = bool(has_unspliced)
        report["velocity_possible"] = bool(has_spliced and has_unspliced)

        if has_spliced:
            report["spliced_total"] = layer_total(adata, "spliced")
        if has_unspliced:
            report["unspliced_total"] = layer_total(adata, "unspliced")
        if has_spliced and has_unspliced:
            s = max(float(report["spliced_total"]), 0.0)
            u = max(float(report["unspliced_total"]), 0.0)
            frac = u / (u + s) if (u + s) > 0 else float("nan")
            report["unspliced_fraction"] = frac
            low = config.get("qc", {}).get("unspliced_fraction_warning_low", 0.05)
            high = config.get("qc", {}).get("unspliced_fraction_warning_high", 0.40)
            if np.isfinite(frac) and frac < low:
                report["warnings"].append(f"Unspliced fraction {frac:.4f} is below warning threshold {low}.")
            if np.isfinite(frac) and frac > high:
                report["warnings"].append(f"Unspliced fraction {frac:.4f} is above warning threshold {high}; check intronic contamination/preprocessing.")
        else:
            report["warnings"].append("Missing spliced/unspliced layers. RNA velocity branch will be skipped.")
            report["recommended_route"] = "non_velocity"

        if "X_umap" not in adata.obsm:
            report["warnings"].append("X_umap not found; UMAP should be computed or recomputed before visualization.")
        if "neighbors" not in adata.uns and "connectivities" not in adata.obsp:
            report["warnings"].append("No neighbors graph detected; CellRank ConnectivityKernel requires neighbors/connectivities or recomputation.")

        adata_sub = subset_if_needed(adata, config, report)
        report["n_cells"] = int(adata_sub.n_obs)
        report["n_genes"] = int(adata_sub.n_vars)

        if report["errors"]:
            report["status"] = "FAIL"
        elif report["warnings"]:
            report["status"] = "WARN"
        else:
            report["status"] = "PASS"

        write_json(report, outdir / "input_validation.json")
        update_manifest(
            config,
            step,
            report["status"],
            inputs=[str(input_h5ad)],
            outputs=[str(outdir / "input_validation.json")],
            warnings=report["warnings"],
            errors=report["errors"],
            route=report.get("recommended_route", "velocity" if report["velocity_possible"] else "non_velocity"),
        )
        print(json.dumps(report, indent=2, ensure_ascii=False))
        if report["status"] == "FAIL":
            raise SystemExit(1)
    except SystemExit:
        raise
    except Exception as exc:
        fail_step(config, step, exc)
        raise


if __name__ == "__main__":
    main()
