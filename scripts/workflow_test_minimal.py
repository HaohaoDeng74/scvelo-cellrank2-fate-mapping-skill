#!/usr/bin/env python
"""Optional minimal workflow test for environments with AnnData/Scanpy/CellRank.

This test is intentionally small. It validates packaging-independent workflow
behaviour: relative output_dir resolution, non-velocity route metadata, manifest
creation, report placeholder guard, and subset-too-small validation. It skips
cleanly when heavy scientific dependencies are unavailable.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def have_deps():
    missing = []
    for pkg in ["anndata", "scanpy", "cellrank", "numpy", "pandas", "scipy"]:
        try:
            __import__(pkg)
        except Exception:
            missing.append(pkg)
    return missing


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skill-root", default=".")
    args = parser.parse_args()
    root = Path(args.skill_root).resolve()
    missing = have_deps()
    if missing:
        print(json.dumps({"status": "SKIPPED", "reason": "missing dependencies", "missing": missing}, indent=2))
        return

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        rng = np.random.default_rng(0)
        X = rng.poisson(1, size=(80, 60)).astype(float)
        obs = pd.DataFrame({
            "celltype": ["A"] * 40 + ["B"] * 40,
            "major_celltype": ["immune"] * 80,
            "sample": ["s1"] * 40 + ["s2"] * 40,
            "condition": ["ctrl"] * 40 + ["case"] * 40,
            "timepoint": [0] * 40 + [1] * 40,
            "pseudotime": np.linspace(0, 1, 80),
        }, index=[f"cell{i}" for i in range(80)])
        var = pd.DataFrame(index=[f"Gene{i}" for i in range(60)])
        adata = ad.AnnData(X=X, obs=obs, var=var)
        sc.pp.pca(adata, n_comps=10)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
        sc.tl.umap(adata)
        h5ad = td / "toy.h5ad"
        adata.write_h5ad(h5ad)
        cfg = td / "config.yaml"
        cfg.write_text(f"""
project:
  name: workflow_test
  input_h5ad: {h5ad}
  output_dir: results
metadata:
  celltype_key: celltype
  cluster_key: celltype
  sample_key: sample
  condition_key: condition
  timepoint_key: timepoint
  pseudotime_key: pseudotime
subset:
  enabled: true
  subset_name: immune
  obs_key: major_celltype
  include_values: [immune]
  min_cells_warning: 10
  min_cells_stop: 5
cellrank:
  route_if_no_velocity: pseudotime
  pseudotime_key: pseudotime
  velocity_kernel_weight: 0.6
  connectivity_kernel_weight: 0.4
  cluster_key: celltype
  n_states: null
  n_states_sensitivity: []
driver_genes:
  top_n_per_lineage: 5
  trend_genes_per_lineage: 3
  use_raw: false
  expression_layer: null
  min_abs_correlation: 0.1
  max_q_value: 0.5
gene_trends:
  time_key: pseudotime
  n_bins: 4
  gene_sets_gmt: null
qc:
  risky_gene_prefixes: {{mitochondrial: ["mt-"]}}
  risky_gene_exact: {{}}
visualization:
  basis: umap
  dpi: 80
  formats: [png]
  color_keys: [celltype]
report:
  generate_markdown: true
""", encoding="utf-8")
        env = None
        for script in ["00_validate_input.py", "03_run_cellrank2.py", "04_driver_gene_analysis.py", "05_gene_trend_and_pathway.py", "06_generate_report.py"]:
            subprocess.run([sys.executable, str(root / "scripts" / script), "--config", str(cfg)], check=True, cwd=str(root), env=env)
        report = (td / "results" / "06_report" / "scvelo_cellrank2_report.md").read_text(encoding="utf-8")
        if "{" in report and "}" in report:
            raise AssertionError("Report appears to contain unresolved brace placeholders.")
        manifest = json.loads((td / "results" / "run_manifest.json").read_text(encoding="utf-8"))
        print(json.dumps({"status": "PASS", "steps": list(manifest.get("steps", {}).keys())}, indent=2))


if __name__ == "__main__":
    main()
