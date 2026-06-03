#!/usr/bin/env python
"""Compute lineage-associated driver candidates.

Native CellRank lineage drivers are preferred when an estimator can be loaded.
Otherwise the module falls back to an explicit expression/fate-probability
correlation analysis. Both outputs are fate-associated candidate lists, not
causal driver evidence.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from utils import bh_adjust, fail_step, get_matrix, load_config, output_dir, risky_gene_category, update_manifest, write_json


def to_dense_2d(x):
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


def find_fate_probabilities(adata, outdir):
    fp_csv = outdir.parent / "03_cellrank_fate" / "fate_probabilities.csv"
    if fp_csv.exists():
        fp = pd.read_csv(fp_csv, index_col=0)
        common = adata.obs_names.intersection(fp.index)
        fp = fp.loc[common]
        return fp, common
    for key in adata.obsm.keys():
        if "fate" in str(key).lower() and "prob" in str(key).lower():
            arr = np.asarray(adata.obsm[key])
            if arr.ndim == 2:
                cols = [f"lineage_{i}" for i in range(arr.shape[1])]
                names = adata.uns.get(f"{key}_names", None)
                if names is not None and len(names) == arr.shape[1]:
                    cols = list(map(str, names))
                return pd.DataFrame(arr, index=adata.obs_names, columns=cols), adata.obs_names
    raise RuntimeError("No fate probability matrix found. Run 03_run_cellrank2.py first.")


def compute_correlations(adata, fate_df, config):
    drv_cfg = config.get("driver_genes", {})
    layer = drv_cfg.get("expression_layer", None)
    use_raw = bool(drv_cfg.get("use_raw", False))
    X = get_matrix(adata, layer=layer, use_raw=use_raw)
    X = to_dense_2d(X)
    var_names = np.asarray(adata.raw.var_names if use_raw and getattr(adata, "raw", None) is not None else adata.var_names)

    # Align rows.
    if X.shape[0] != fate_df.shape[0]:
        raise RuntimeError(f"Expression matrix rows ({X.shape[0]}) do not match fate probabilities ({fate_df.shape[0]}).")

    rows = []
    X = X.astype(float, copy=False)
    gene_means = np.nanmean(X, axis=0)
    gene_stds = np.nanstd(X, axis=0)
    valid_genes = gene_stds > 0

    for lineage in fate_df.columns:
        y = fate_df[lineage].to_numpy(dtype=float)
        y_std = np.nanstd(y)
        if y_std == 0 or not np.isfinite(y_std):
            continue
        y_center = y - np.nanmean(y)
        for j, gene in enumerate(var_names):
            if not valid_genes[j]:
                continue
            x = X[:, j]
            x_center = x - gene_means[j]
            denom = (len(y) - 1) * gene_stds[j] * y_std
            corr = float(np.nansum(x_center * y_center) / denom) if denom else np.nan
            if not np.isfinite(corr):
                continue
            # Approximate two-sided p-value using t distribution if scipy is available.
            n = len(y)
            try:
                from scipy import stats
                t = corr * np.sqrt((n - 2) / max(1e-12, 1 - corr * corr))
                pval = float(2 * stats.t.sf(abs(t), df=n - 2))
            except Exception:
                pval = np.nan
            rows.append({
                "lineage": str(lineage),
                "gene": str(gene),
                "correlation": corr,
                "abs_correlation": abs(corr),
                "p_value": pval,
                "mean_expression": float(gene_means[j]),
                "risk_category": risky_gene_category(str(gene), config),
            })
    df = pd.DataFrame(rows)
    if len(df) and df["p_value"].notna().any():
        df["q_value"] = np.nan
        for lineage, idx in df.groupby("lineage").groups.items():
            p = df.loc[idx, "p_value"].to_numpy(dtype=float)
            mask = np.isfinite(p)
            q = np.full(len(p), np.nan)
            q[mask] = bh_adjust(p[mask])
            df.loc[idx, "q_value"] = q
    else:
        df["q_value"] = np.nan

    df["method"] = "correlation_fallback"
    df["high_risk_gene"] = df["risk_category"].astype(str).ne("")
    df["candidate_label"] = np.where(df["high_risk_gene"], "technical_risk_fate_associated", "fate_associated_candidate")
    df = df.sort_values(["lineage", "abs_correlation"], ascending=[True, False])
    return df


def try_cellrank_driver_export(adata, config, warnings):
    """Try native CellRank lineage driver analysis from a saved GPCCA estimator.

    The native path is preferred when a pickled estimator is available and the
    installed CellRank version exposes compute_lineage_drivers(). Otherwise the
    workflow falls back to an explicit correlation analysis.
    """
    estimator_path = output_dir(config) / "03_cellrank_fate" / "gpcca_estimator.pkl"
    if not estimator_path.exists():
        warnings.append("CellRank native driver analysis skipped: gpcca_estimator.pkl not found; using correlation fallback.")
        return None
    try:
        with open(estimator_path, "rb") as fh:
            estimator = pickle.load(fh)
        if not hasattr(estimator, "compute_lineage_drivers"):
            warnings.append("CellRank estimator lacks compute_lineage_drivers(); using correlation fallback.")
            return None
        native = estimator.compute_lineage_drivers()
        if native is None:
            warnings.append("CellRank compute_lineage_drivers() returned None; using correlation fallback.")
            return None
        if isinstance(native, pd.DataFrame):
            df = native.copy()
        else:
            df = pd.DataFrame(native)
        if df.empty:
            warnings.append("CellRank native driver table is empty; using correlation fallback.")
            return None
        df = df.reset_index().rename(columns={"index": "gene"})
        # Normalize common CellRank table shapes into a long table when possible.
        if "lineage" not in df.columns:
            value_cols = [c for c in df.columns if c != "gene"]
            if value_cols:
                df = df.melt(id_vars=["gene"], value_vars=value_cols, var_name="lineage", value_name="score")
        if "correlation" not in df.columns and "score" in df.columns:
            df["correlation"] = pd.to_numeric(df["score"], errors="coerce")
        if "abs_correlation" not in df.columns and "correlation" in df.columns:
            df["abs_correlation"] = df["correlation"].abs()
        if "p_value" not in df.columns:
            for c in ["pval", "pvals", "p"]:
                if c in df.columns:
                    df["p_value"] = pd.to_numeric(df[c], errors="coerce")
                    break
        if "q_value" not in df.columns:
            for c in ["qval", "qvals", "pval_adj", "pvals_adj"]:
                if c in df.columns:
                    df["q_value"] = pd.to_numeric(df[c], errors="coerce")
                    break
        if "q_value" not in df.columns:
            df["q_value"] = np.nan
        if "mean_expression" not in df.columns:
            df["mean_expression"] = np.nan
        if "lineage" not in df.columns:
            df["lineage"] = "unknown_lineage"
        df["method"] = "cellrank_native_compute_lineage_drivers"
        df["risk_category"] = df["gene"].map(lambda g: risky_gene_category(str(g), config))
        return df
    except Exception as exc:
        warnings.append(f"CellRank native driver analysis failed: {exc}; using correlation fallback.")
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    step = "04_driver_gene_analysis"
    outdir = output_dir(config) / "04_driver_genes"
    outdir.mkdir(parents=True, exist_ok=True)
    warnings = []
    outputs = []

    try:
        input_h5ad = args.input or str(output_dir(config) / "03_cellrank_fate" / "adata_cellrank2.h5ad")
        adata = ad.read_h5ad(input_h5ad)

        fate_df, common = find_fate_probabilities(adata, outdir)
        if len(common) != adata.n_obs:
            adata = adata[common].copy()

        native = try_cellrank_driver_export(adata, config, warnings)
        if native is not None:
            df = native
            if "risk_category" not in df.columns:
                df["risk_category"] = df["gene"].map(lambda g: risky_gene_category(str(g), config))
            if "abs_correlation" not in df.columns:
                score_col = "correlation" if "correlation" in df.columns else None
                df["abs_correlation"] = pd.to_numeric(df[score_col], errors="coerce").abs() if score_col else np.nan
        else:
            df = compute_correlations(adata, fate_df, config)

        if "method" not in df.columns:
            df["method"] = "unknown"
        if "q_value" not in df.columns:
            df["q_value"] = np.nan
        if "p_value" not in df.columns:
            df["p_value"] = np.nan
        df["high_risk_gene"] = df["risk_category"].astype(str).ne("")
        df["candidate_label"] = np.where(df["high_risk_gene"], "technical_risk_fate_associated", "fate_associated_candidate")

        min_abs = config.get("driver_genes", {}).get("min_abs_correlation", 0.1)
        max_q = config.get("driver_genes", {}).get("max_q_value", 0.1)
        df["passes_default_filter"] = (df["abs_correlation"] >= min_abs) & (
            df["q_value"].isna() | (df["q_value"] <= max_q)
        )

        out_csv = outdir / "lineage_driver_candidates.csv"
        df.to_csv(out_csv, index=False)
        outputs.append(str(out_csv))

        top_n = int(config.get("driver_genes", {}).get("top_n_per_lineage", 50))
        top = df[df["passes_default_filter"]].groupby("lineage", group_keys=False).head(top_n)
        top.to_csv(outdir / "top_lineage_driver_candidates.csv", index=False)
        outputs.append(str(outdir / "top_lineage_driver_candidates.csv"))

        methods = sorted(map(str, df["method"].dropna().unique().tolist())) if "method" in df.columns else []
        if "cellrank_native_compute_lineage_drivers" in methods and "correlation_fallback" in methods:
            interpretation = "Mixed native CellRank and correlation-fallback fate-associated candidate genes; not causal drivers."
        elif "cellrank_native_compute_lineage_drivers" in methods:
            interpretation = "Native CellRank compute_lineage_drivers fate-associated candidate genes; not causal drivers."
        elif "correlation_fallback" in methods:
            interpretation = "Correlation-fallback fate-associated candidate genes; not native CellRank driver analysis and not causal drivers."
        else:
            interpretation = "Fate-associated candidate genes; method unavailable; not causal drivers."
        summary = {
            "status": "PASS",
            "n_rows": int(len(df)),
            "n_lineages": int(df["lineage"].nunique()) if len(df) else 0,
            "n_high_risk_genes": int(df["high_risk_gene"].sum()) if len(df) else 0,
            "methods": methods,
            "warnings": warnings,
            "interpretation": interpretation,
        }
        write_json(summary, outdir / "driver_gene_summary.json")
        outputs.append(str(outdir / "driver_gene_summary.json"))

        update_manifest(config, step, "PASS" if not warnings else "WARN", inputs=[str(input_h5ad)], outputs=outputs, warnings=warnings, summary=summary)
    except Exception as exc:
        fail_step(config, step, exc, warnings=warnings)
        raise


if __name__ == "__main__":
    main()
