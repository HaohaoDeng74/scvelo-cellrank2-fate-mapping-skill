#!/usr/bin/env python
"""Shared utilities for the scVelo + CellRank2 fate mapping skill."""

from __future__ import annotations

import json
import sys
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import yaml


def load_config(path: str | Path) -> dict:
    path = Path(path).resolve()
    with open(path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise ValueError(f"Config did not parse to a dictionary: {path}")
    cfg["_config_path"] = str(path)
    cfg["_config_dir"] = str(path.parent)
    return cfg


def output_dir(config: dict) -> Path:
    out = Path(config["project"]["output_dir"])
    if out.is_absolute():
        return out
    base = Path(config.get("_config_dir", ".")).resolve()
    return base / out


def read_json(path: str | Path, default: Any = None) -> Any:
    path = Path(path)
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(obj: Any, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def manifest_path(config: dict) -> Path:
    return output_dir(config) / "run_manifest.json"


def init_manifest(config: dict) -> dict:
    path = manifest_path(config)
    if path.exists():
        manifest = read_json(path, default={})
    else:
        manifest = {
            "project": config.get("project", {}).get("name", "unnamed"),
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "steps": {},
            "warnings": [],
            "errors": [],
        }
    manifest.setdefault("steps", {})
    manifest.setdefault("warnings", [])
    manifest.setdefault("errors", [])
    return manifest


def package_versions() -> dict:
    versions = {"python": sys.version.split()[0]}
    for pkg in ["anndata", "scanpy", "scvelo", "cellrank", "numpy", "pandas", "scipy"]:
        try:
            mod = __import__(pkg)
            versions[pkg] = getattr(mod, "__version__", "unknown")
        except Exception:
            versions[pkg] = "not_installed"
    return versions


def update_manifest(config: dict, step: str, status: str, **kwargs: Any) -> None:
    manifest = init_manifest(config)
    entry = {
        "status": status,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "versions": package_versions(),
    }
    entry.update(kwargs)
    manifest["steps"][step] = entry
    for w in kwargs.get("warnings", []) or []:
        manifest["warnings"].append({"step": step, "warning": w})
    if status in {"ERROR", "FAIL"}:
        manifest["errors"].append({"step": step, "error": kwargs.get("error", "unknown")})
    write_json(manifest, manifest_path(config))


def fail_step(config: dict, step: str, exc: Exception, warnings: list[str] | None = None) -> None:
    update_manifest(
        config,
        step,
        "ERROR",
        error=str(exc),
        traceback=traceback.format_exc(),
        warnings=warnings or [],
    )


def dense_sum(x, axis=None):
    """Sum dense or sparse matrices and return a numpy array/scalar."""
    y = x.sum(axis=axis)
    if hasattr(y, "A1"):
        return y.A1
    if hasattr(y, "A"):
        return np.asarray(y.A).ravel()
    return np.asarray(y).ravel() if axis is not None else float(y)


def layer_total(adata, layer: str) -> float:
    if layer not in adata.layers:
        return float("nan")
    return float(dense_sum(adata.layers[layer]).sum())


def get_matrix(adata, layer: str | None = None, use_raw: bool = False):
    if use_raw and getattr(adata, "raw", None) is not None:
        return adata.raw.X
    if layer and layer in adata.layers:
        return adata.layers[layer]
    return adata.X


def bh_adjust(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        val = ranked[i] * n / (i + 1)
        prev = min(prev, val)
        q[order[i]] = prev
    return np.minimum(q, 1.0)


def risky_gene_category(gene: str, config: dict) -> str:
    qc = config.get("qc", {})
    for category, prefixes in qc.get("risky_gene_prefixes", {}).items():
        if any(str(gene).startswith(p) for p in prefixes):
            return category
    for category, genes in qc.get("risky_gene_exact", {}).items():
        if str(gene) in set(map(str, genes)):
            return category
    return ""


def sanitize_filename(text: str) -> str:
    return "".join(c if c.isalnum() or c in "-._" else "_" for c in str(text))[:120]
