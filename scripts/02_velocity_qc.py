#!/usr/bin/env python
"""Generate scVelo QC and velocity plots with verified output paths."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scvelo as scv

from utils import fail_step, load_config, output_dir, update_manifest, write_json


def capture_new_files(outdir: Path, before: set[Path]) -> list[str]:
    after = {p for p in outdir.glob("**/*") if p.is_file()}
    return [str(p) for p in sorted(after - before)]


def run_plot(outdir: Path, warnings: list[str], label: str, func, *args, **kwargs) -> list[str]:
    before = {p for p in outdir.glob("**/*") if p.is_file()}
    try:
        func(*args, **kwargs)
        files = capture_new_files(outdir, before)
        if not files:
            warnings.append(f"{label} ran but no new plot file was detected; check scVelo save naming.")
        return files
    except Exception as exc:
        warnings.append(f"{label} failed: {exc}")
        return []


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    step = "02_velocity_qc"
    outdir = output_dir(config) / "02_velocity_qc"
    outdir.mkdir(parents=True, exist_ok=True)
    warnings = []
    outputs = []
    status = "PASS"

    try:
        input_h5ad = args.input or str(output_dir(config) / "01_velocity" / "adata_scvelo.h5ad")
        adata = ad.read_h5ad(input_h5ad)
        if adata.uns.get("fate_mapping_route") == "non_velocity" or "velocity" not in adata.layers:
            warnings.append("No velocity layer detected; velocity QC plotting skipped.")
            write_json({"status": "SKIPPED", "warnings": warnings}, outdir / "velocity_qc_status.json")
            update_manifest(config, step, "SKIPPED", inputs=[str(input_h5ad)], outputs=[str(outdir / "velocity_qc_status.json")], warnings=warnings)
            return

        basis = config.get("visualization", {}).get("basis", "umap")
        formats = config.get("visualization", {}).get("formats", ["png"])
        color_keys = [c for c in config.get("visualization", {}).get("color_keys", []) if c in adata.obs.columns]
        scv.settings.figdir = str(outdir)
        scv.settings.set_figure_params(dpi=config.get("visualization", {}).get("dpi", 300))

        for fmt in formats:
            outputs.extend(run_plot(outdir, warnings, "velocity_embedding_stream", scv.pl.velocity_embedding_stream, adata, basis=basis, color=color_keys or None, save=f"_stream.{fmt}", show=False))
            outputs.extend(run_plot(outdir, warnings, "velocity_embedding_grid", scv.pl.velocity_embedding_grid, adata, basis=basis, color=color_keys[0] if color_keys else None, save=f"_grid.{fmt}", show=False))
            outputs.extend(run_plot(outdir, warnings, "latent_time_scatter", scv.pl.scatter, adata, color=[c for c in ["latent_time", "velocity_length", "velocity_confidence"] if c in adata.obs.columns], save=f"_latent_velocity_qc.{fmt}", show=False))

        try:
            scv.tl.velocity_confidence(adata)
        except Exception as exc:
            warnings.append(f"velocity_confidence computation failed: {exc}")
            status = "WARN"

        summary = {
            "status": status if not warnings else "WARN",
            "median_velocity_confidence": float(adata.obs["velocity_confidence"].median()) if "velocity_confidence" in adata.obs else None,
            "median_velocity_length": float(adata.obs["velocity_length"].median()) if "velocity_length" in adata.obs else None,
            "plot_files_detected": outputs,
            "warnings": warnings,
        }
        write_json(summary, outdir / "velocity_qc_summary.json")
        outputs.append(str(outdir / "velocity_qc_summary.json"))
        adata.write_h5ad(outdir / "adata_velocity_qc.h5ad")
        outputs.append(str(outdir / "adata_velocity_qc.h5ad"))
        update_manifest(config, step, "PASS" if not warnings else "WARN", inputs=[str(input_h5ad)], outputs=outputs, warnings=warnings, summary=summary)
    except Exception as exc:
        fail_step(config, step, exc, warnings=warnings)
        raise


if __name__ == "__main__":
    main()
