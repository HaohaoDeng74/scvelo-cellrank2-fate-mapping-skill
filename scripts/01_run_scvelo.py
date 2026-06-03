#!/usr/bin/env python
"""Run scVelo preprocessing, velocity, velocity graph, latent time, and confidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import scanpy as sc
import scvelo as scv

from utils import fail_step, load_config, output_dir, update_manifest, write_json


def subset_adata(adata, config):
    subset_cfg = config.get("subset", {})
    if not subset_cfg.get("enabled", False):
        return adata.copy()
    obs_key = subset_cfg.get("obs_key")
    include_values = subset_cfg.get("include_values", [])
    if obs_key not in adata.obs.columns:
        raise KeyError(f"Subset obs_key `{obs_key}` not found.")
    return adata[adata.obs[obs_key].isin(include_values)].copy()


def maybe_recompute_embedding(adata, config):
    pp = config.get("preprocessing", {})
    if not pp.get("recompute_pca_neighbors_umap", False):
        return
    sc.pp.highly_variable_genes(adata, n_top_genes=pp.get("n_hvgs", 3000), subset=False, flavor="seurat_v3")
    if "highly_variable" in adata.var.columns and int(adata.var["highly_variable"].sum()) > 50:
        adata_pca_view = adata[:, adata.var["highly_variable"]].copy()
    else:
        adata_pca_view = adata
    sc.pp.pca(adata_pca_view, n_comps=pp.get("n_pcs", 30), random_state=pp.get("random_state", 0))
    adata.obsm["X_pca"] = adata_pca_view.obsm["X_pca"]
    sc.pp.neighbors(adata, n_pcs=pp.get("n_pcs", 30), n_neighbors=pp.get("n_neighbors", 30), random_state=pp.get("random_state", 0))
    sc.tl.umap(adata, random_state=pp.get("random_state", 0))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    step = "01_run_scvelo"
    outdir = output_dir(config)
    velocity_dir = outdir / "01_velocity"
    velocity_dir.mkdir(parents=True, exist_ok=True)
    status = "PASS"
    warnings = []
    outputs = []

    try:
        input_h5ad = Path(config["project"]["input_h5ad"])
        output_h5ad = velocity_dir / "adata_scvelo.h5ad"
        if output_h5ad.exists() and not args.force:
            warnings.append(f"Output exists and --force not set; reusing {output_h5ad}.")
            update_manifest(config, step, "SKIPPED", inputs=[str(input_h5ad)], outputs=[str(output_h5ad)], warnings=warnings)
            print("\n".join(warnings))
            return

        adata = ad.read_h5ad(input_h5ad)
        adata = subset_adata(adata, config)

        if "spliced" not in adata.layers or "unspliced" not in adata.layers:
            warnings.append("Missing spliced/unspliced layers. scVelo velocity branch skipped; copied subset for non-velocity CellRank route.")
            adata.uns["fate_mapping_route"] = "non_velocity"
            adata.write_h5ad(output_h5ad)
            write_json({"route": "non_velocity", "warnings": warnings}, velocity_dir / "scvelo_status.json")
            update_manifest(config, step, "SKIPPED", inputs=[str(input_h5ad)], outputs=[str(output_h5ad)], warnings=warnings, route="non_velocity")
            return

        maybe_recompute_embedding(adata, config)
        scv_cfg = config["scvelo"]

        scv.pp.filter_and_normalize(
            adata,
            min_shared_counts=scv_cfg.get("min_shared_counts", 20),
            n_top_genes=scv_cfg.get("n_top_genes", 2000),
        )

        scv.pp.moments(
            adata,
            n_pcs=scv_cfg.get("n_pcs", 30),
            n_neighbors=scv_cfg.get("n_neighbors", 30),
        )

        mode = scv_cfg.get("mode", "dynamical")
        velocity_mode_used = mode
        try:
            if mode == "dynamical":
                scv.tl.recover_dynamics(adata, n_jobs=scv_cfg.get("n_jobs", 8))
            scv.tl.velocity(adata, mode=mode)
        except Exception as exc:
            fallback = scv_cfg.get("fallback_mode", "stochastic")
            warnings.append(f"Primary velocity mode `{mode}` failed: {exc}. Falling back to `{fallback}`.")
            scv.tl.velocity(adata, mode=fallback)
            adata.uns["velocity_mode_fallback"] = fallback
            velocity_mode_used = fallback
            status = "WARN"

        scv.tl.velocity_graph(adata)
        try:
            scv.tl.latent_time(adata)
        except Exception as exc:
            warnings.append(f"latent_time failed: {exc}")
            status = "WARN"
        try:
            scv.tl.velocity_confidence(adata)
        except Exception as exc:
            warnings.append(f"velocity_confidence failed: {exc}")
            status = "WARN"

        n_velocity_genes = int(adata.var["velocity_genes"].sum()) if "velocity_genes" in adata.var.columns else 0
        if n_velocity_genes < config.get("qc", {}).get("min_velocity_genes_warning", 200):
            warnings.append(f"Only {n_velocity_genes} velocity genes detected.")
            status = "WARN"

        adata.uns["velocity_mode_used"] = velocity_mode_used
        adata.uns["fate_mapping_route"] = "velocity"
        adata.write_h5ad(output_h5ad)
        outputs.append(str(output_h5ad))
        write_json({"status": status, "velocity_mode_used": velocity_mode_used, "n_velocity_genes": n_velocity_genes, "warnings": warnings}, velocity_dir / "scvelo_status.json")
        outputs.append(str(velocity_dir / "scvelo_status.json"))

        update_manifest(
            config,
            step,
            status,
            inputs=[str(input_h5ad)],
            outputs=outputs,
            warnings=warnings,
            route="velocity",
            velocity_mode_used=velocity_mode_used,
            n_velocity_genes=n_velocity_genes,
        )
    except Exception as exc:
        fail_step(config, step, exc, warnings=warnings)
        raise


if __name__ == "__main__":
    main()
