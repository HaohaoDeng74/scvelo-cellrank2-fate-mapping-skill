#!/usr/bin/env python
"""Run CellRank kernels, GPCCA, terminal/initial states, and fate probabilities.

The implementation deliberately exports results from the estimator first and
uses AnnData key discovery only as a fallback, because CellRank storage keys can
vary by version.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from utils import fail_step, load_config, output_dir, update_manifest, write_json


def get_cellrank():
    import cellrank as cr
    return cr


def build_kernel(adata, config, warnings):
    cr = get_cellrank()
    cr_cfg = config.get("cellrank", {})
    can_use_velocity = "velocity_graph" in adata.uns or "velocity" in adata.layers
    route = "velocity" if can_use_velocity else cr_cfg.get("route_if_no_velocity", "connectivity")

    if can_use_velocity:
        vk = cr.kernels.VelocityKernel(adata)
        vk.compute_transition_matrix()
        ck = cr.kernels.ConnectivityKernel(adata)
        ck.compute_transition_matrix()
        vw = float(cr_cfg.get("velocity_kernel_weight", 0.6))
        cw = float(cr_cfg.get("connectivity_kernel_weight", 0.4))
        kernel = vw * vk + cw * ck
        desc = f"{vw:.2f} * VelocityKernel + {cw:.2f} * ConnectivityKernel"
        return kernel, desc, "velocity", "rna_velocity"

    pseudo_key = cr_cfg.get("pseudotime_key") or config.get("metadata", {}).get("pseudotime_key")
    time_key = cr_cfg.get("time_key") or config.get("metadata", {}).get("timepoint_key")
    if route == "pseudotime" and pseudo_key and pseudo_key in adata.obs.columns:
        pk = cr.kernels.PseudotimeKernel(adata, time_key=pseudo_key)
        pk.compute_transition_matrix()
        return pk, f"PseudotimeKernel({pseudo_key}); no RNA velocity", "non_velocity_pseudotime", "pseudotime"
    if route == "realtime" and time_key and time_key in adata.obs.columns:
        try:
            rk = cr.kernels.RealTimeKernel(adata, time_key=time_key)
            rk.compute_transition_matrix()
            return rk, f"RealTimeKernel({time_key}); no RNA velocity", "non_velocity_realtime", "real_time"
        except Exception as exc:
            warnings.append(f"RealTimeKernel failed, using ConnectivityKernel: {exc}")
    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()
    return ck, "ConnectivityKernel only; no RNA velocity/directional prior", "non_velocity_connectivity", "connectivity_only"


def run_gpcca(kernel, adata, config, warnings, n_states_override=None):
    cr = get_cellrank()
    cr_cfg = config.get("cellrank", {})
    cluster_key = cr_cfg.get("cluster_key") or config.get("metadata", {}).get("cluster_key", "celltype")
    cluster_key = cluster_key if cluster_key in adata.obs.columns else None
    g = cr.estimators.GPCCA(kernel)
    n_states = n_states_override if n_states_override is not None else cr_cfg.get("n_states", None)

    if hasattr(g, "compute_macrostates"):
        kwargs = {"cluster_key": cluster_key} if cluster_key else {}
        if n_states is not None:
            kwargs["n_states"] = n_states
        g.compute_macrostates(**kwargs)
    else:
        kwargs = {"cluster_key": cluster_key} if cluster_key else {}
        if n_states is not None:
            kwargs["n_states"] = n_states
        g.fit(**kwargs)

    if cr_cfg.get("manual_terminal_states") and hasattr(g, "set_terminal_states"):
        try:
            g.set_terminal_states(cr_cfg["manual_terminal_states"])
        except Exception as exc:
            warnings.append(f"Manual terminal state setting failed: {exc}")
    else:
        for method in ["predict_terminal_states", "set_terminal_states"]:
            if hasattr(g, method):
                try:
                    getattr(g, method)()
                    break
                except TypeError:
                    # Some versions require explicit states for set_terminal_states.
                    continue
                except Exception as exc:
                    warnings.append(f"{method} failed: {exc}")

    if cr_cfg.get("manual_initial_states") and hasattr(g, "set_initial_states"):
        try:
            g.set_initial_states(cr_cfg["manual_initial_states"])
        except Exception as exc:
            warnings.append(f"Manual initial state setting failed: {exc}")
    else:
        for method in ["predict_initial_states", "set_initial_states"]:
            if hasattr(g, method):
                try:
                    getattr(g, method)()
                    break
                except TypeError:
                    continue
                except Exception as exc:
                    warnings.append(f"{method} failed: {exc}")

    if not hasattr(g, "compute_fate_probabilities"):
        raise RuntimeError("GPCCA object lacks compute_fate_probabilities(). Check CellRank version.")
    g.compute_fate_probabilities()
    return g


def _as_dataframe(obj, index=None, prefix="value"):
    if obj is None:
        return None
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    if isinstance(obj, pd.Series):
        return obj.to_frame()
    try:
        arr = np.asarray(obj)
    except Exception:
        try:
            arr = np.asarray(obj.X)
        except Exception:
            return None
    if arr.ndim == 0:
        return pd.DataFrame({prefix: [arr.item()]})
    if arr.ndim == 1:
        return pd.DataFrame({prefix: arr}, index=index if index is not None and len(index) == arr.shape[0] else None)
    cols = [f"{prefix}_{i}" for i in range(arr.shape[1])]
    return pd.DataFrame(arr, index=index if index is not None and len(index) == arr.shape[0] else None, columns=cols)


def _lineage_names(obj, n):
    for attr in ["names", "columns", "categories"]:
        val = getattr(obj, attr, None)
        try:
            if val is not None and len(val) == n:
                return list(map(str, val))
        except Exception:
            pass
    return [f"lineage_{i}" for i in range(n)]


def _get_estimator_attr(g, names):
    for name in names:
        if hasattr(g, name):
            try:
                val = getattr(g, name)
                if val is not None:
                    return val, name
            except Exception:
                continue
    return None, None


def export_from_estimator(g, adata, outdir, warnings):
    outputs = []
    flags = {
        "fate_probabilities_exported": False,
        "terminal_states_exported": False,
        "initial_states_exported": False,
        "macrostates_exported": False,
        "coarse_transition_matrix_exported": False,
    }

    fp, attr = _get_estimator_attr(g, ["fate_probabilities", "lineage_drivers"])  # second attr unlikely; harmless fallback.
    if attr == "lineage_drivers":
        fp = None
    if fp is not None:
        df = _as_dataframe(fp, index=adata.obs_names, prefix="lineage")
        if df is not None and df.shape[0] == adata.n_obs:
            if df.shape[1] > 1:
                df.columns = _lineage_names(fp, df.shape[1])
            df.to_csv(outdir / "fate_probabilities.csv")
            outputs.append(str(outdir / "fate_probabilities.csv"))
            flags["fate_probabilities_exported"] = True

    for label, attrs, fname in [
        ("terminal", ["terminal_states", "terminal_states_memberships"], "terminal_states.csv"),
        ("initial", ["initial_states", "initial_states_memberships"], "initial_states.csv"),
        ("macro", ["macrostates", "macrostates_memberships"], "macrostates.csv"),
    ]:
        obj, attr = _get_estimator_attr(g, attrs)
        if obj is not None:
            df = _as_dataframe(obj, index=adata.obs_names, prefix=label)
            if df is not None:
                df.to_csv(outdir / fname)
                outputs.append(str(outdir / fname))
                if label == "terminal":
                    flags["terminal_states_exported"] = True
                elif label == "initial":
                    flags["initial_states_exported"] = True
                else:
                    flags["macrostates_exported"] = True

    obj, attr = _get_estimator_attr(g, ["coarse_T", "coarse_transition_matrix"])
    if obj is not None:
        df = _as_dataframe(obj, prefix="state")
        if df is not None:
            df.to_csv(outdir / "coarse_transition_matrix.csv")
            outputs.append(str(outdir / "coarse_transition_matrix.csv"))
            flags["coarse_transition_matrix_exported"] = True
    return outputs, flags


def export_from_adata(adata, outdir, warnings, flags):
    outputs = []
    if not flags.get("fate_probabilities_exported"):
        for key in list(adata.obsm.keys()):
            if "fate" in str(key).lower() and "prob" in str(key).lower():
                arr = np.asarray(adata.obsm[key])
                if arr.ndim == 2:
                    cols = [f"lineage_{i}" for i in range(arr.shape[1])]
                    names = adata.uns.get(f"{key}_names", None)
                    if names is not None and len(names) == arr.shape[1]:
                        cols = list(map(str, names))
                    pd.DataFrame(arr, index=adata.obs_names, columns=cols).to_csv(outdir / "fate_probabilities.csv")
                    outputs.append(str(outdir / "fate_probabilities.csv"))
                    flags["fate_probabilities_exported"] = True
                    break

    for state_type in ["terminal", "initial"]:
        flag_name = f"{state_type}_states_exported"
        if not flags.get(flag_name):
            candidates = [c for c in adata.obs.columns if state_type in str(c).lower()]
            if candidates:
                adata.obs[candidates].to_csv(outdir / f"{state_type}_states.csv")
                outputs.append(str(outdir / f"{state_type}_states.csv"))
                flags[flag_name] = True
            else:
                warnings.append(f"No {state_type} state columns detected in adata.obs after CellRank.")
    return outputs, flags


def fate_entropy(fp_df: pd.DataFrame | None) -> float | None:
    if fp_df is None or fp_df.empty:
        return None
    arr = fp_df.to_numpy(dtype=float)
    arr = np.clip(arr, 1e-12, 1.0)
    row_sum = arr.sum(axis=1, keepdims=True)
    arr = arr / np.where(row_sum == 0, 1, row_sum)
    ent = -(arr * np.log(arr)).sum(axis=1)
    if arr.shape[1] > 1:
        ent = ent / np.log(arr.shape[1])
    return float(np.nanmedian(ent))


def load_fp(path: Path):
    if path.exists():
        try:
            return pd.read_csv(path, index_col=0)
        except Exception:
            return None
    return None



def fate_assignment(fp_df: pd.DataFrame | None):
    if fp_df is None or fp_df.empty:
        return None
    arr = fp_df.to_numpy(dtype=float)
    if arr.ndim != 2 or arr.shape[1] == 0:
        return None
    cols = list(map(str, fp_df.columns))
    return pd.Series([cols[i] for i in np.nanargmax(arr, axis=1)], index=fp_df.index, name="max_fate_assignment")


def assignment_stability(main_fp: pd.DataFrame | None, fp: pd.DataFrame | None):
    main_assign = fate_assignment(main_fp)
    test_assign = fate_assignment(fp)
    if main_assign is None or test_assign is None:
        return {"assignment_match_fraction": None, "assignment_ami_vs_main": None}
    common = main_assign.index.intersection(test_assign.index)
    if len(common) == 0:
        return {"assignment_match_fraction": None, "assignment_ami_vs_main": None}
    a = main_assign.loc[common].astype(str).to_numpy()
    b = test_assign.loc[common].astype(str).to_numpy()
    match = float(np.mean(a == b)) if len(a) else None
    ami = None
    try:
        from sklearn.metrics import adjusted_mutual_info_score
        ami = float(adjusted_mutual_info_score(a, b))
    except Exception:
        ami = None
    return {"assignment_match_fraction": match, "assignment_ami_vs_main": ami}


def export_assignment_composition(adata, fp_df: pd.DataFrame | None, out_path: Path, config):
    assign = fate_assignment(fp_df)
    if assign is None:
        return None
    meta_cfg = config.get("metadata", {})
    group_keys = [meta_cfg.get("celltype_key"), meta_cfg.get("cluster_key"), meta_cfg.get("condition_key"), meta_cfg.get("timepoint_key")]
    group_keys = [g for g in group_keys if g and g in adata.obs.columns]
    rows = []
    obs = adata.obs.loc[assign.index.intersection(adata.obs_names)].copy()
    assign = assign.loc[obs.index]
    for key in group_keys:
        tab = pd.crosstab(obs[key].astype(str), assign.astype(str), normalize="index")
        counts = pd.crosstab(obs[key].astype(str), assign.astype(str))
        for group in tab.index:
            for fate in tab.columns:
                rows.append({
                    "group_key": key,
                    "group": str(group),
                    "fate_assignment": str(fate),
                    "fraction": float(tab.loc[group, fate]),
                    "count": int(counts.loc[group, fate]),
                })
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(out_path, index=False)
        return str(out_path)
    return None

def sensitivity_analysis(adata, config, outdir, warnings, main_fp):
    rows = []
    settings = config.get("cellrank", {}).get("n_states_sensitivity", []) or []
    if not settings:
        return rows
    main_cols = set(map(str, main_fp.columns)) if main_fp is not None else set()
    for setting in settings:
        setting_label = str(setting)
        try:
            adata_tmp = adata.copy()
            local_warnings = []
            kernel_tmp, desc, route, direction = build_kernel(adata_tmp, config, local_warnings)
            g_tmp = run_gpcca(kernel_tmp, adata_tmp, config, local_warnings, n_states_override=setting)
            try:
                adata_tmp = g_tmp.adata
            except Exception:
                pass
            tmp_dir = outdir / "sensitivity" / (setting_label.replace(" ", "").replace("[", "").replace("]", "").replace(",", "_"))
            tmp_dir.mkdir(parents=True, exist_ok=True)
            _, flags = export_from_estimator(g_tmp, adata_tmp, tmp_dir, local_warnings)
            _, flags = export_from_adata(adata_tmp, tmp_dir, local_warnings, flags)
            fp = load_fp(tmp_dir / "fate_probabilities.csv")
            cols = set(map(str, fp.columns)) if fp is not None else set()
            jaccard = (len(main_cols & cols) / len(main_cols | cols)) if main_cols and cols else None
            comp_path = export_assignment_composition(adata_tmp, fp, tmp_dir / "fate_assignment_composition.csv", config)
            stability = assignment_stability(main_fp, fp)
            rows.append({
                "n_states_setting": setting_label,
                "status": "PASS",
                "route": route,
                "directionality_supported_by": direction,
                "n_lineages": int(fp.shape[1]) if fp is not None else None,
                "fate_probability_entropy_median": fate_entropy(fp),
                "lineage_column_jaccard_vs_main": jaccard,
                **stability,
                "composition_table": comp_path,
                "warnings": "; ".join(local_warnings),
            })
        except Exception as exc:
            rows.append({"n_states_setting": setting_label, "status": "FAIL", "error": str(exc)})
    pd.DataFrame(rows).to_csv(outdir / "n_states_sensitivity.csv", index=False)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    step = "03_run_cellrank"
    outdir = output_dir(config) / "03_cellrank_fate"
    outdir.mkdir(parents=True, exist_ok=True)
    warnings = []
    outputs = []
    status = "PASS"

    try:
        input_h5ad = args.input or str(output_dir(config) / "01_velocity" / "adata_scvelo.h5ad")
        if not Path(input_h5ad).exists():
            input_h5ad = config["project"]["input_h5ad"]
            warnings.append("scVelo output not found; using original input for CellRank non-velocity route.")
        adata = ad.read_h5ad(input_h5ad)

        kernel, desc, route, directionality = build_kernel(adata, config, warnings)
        adata.uns["cellrank_kernel_description"] = desc
        adata.uns["cellrank_route"] = route
        adata.uns["directionality_supported_by"] = directionality

        g = run_gpcca(kernel, adata, config, warnings)
        try:
            with open(outdir / "gpcca_estimator.pkl", "wb") as fh:
                pickle.dump(g, fh)
            outputs.append(str(outdir / "gpcca_estimator.pkl"))
        except Exception as exc:
            warnings.append(f"Estimator pickle serialization skipped: {exc}")
        try:
            if hasattr(g, "write"):
                g.write(outdir / "gpcca_estimator.cellrank")
                outputs.append(str(outdir / "gpcca_estimator.cellrank"))
        except Exception as exc:
            warnings.append(f"Estimator native serialization skipped: {exc}")

        try:
            adata = g.adata
        except Exception:
            pass

        est_outputs, export_flags = export_from_estimator(g, adata, outdir, warnings)
        outputs.extend(est_outputs)
        adata_outputs, export_flags = export_from_adata(adata, outdir, warnings, export_flags)
        outputs.extend(adata_outputs)

        adata.write_h5ad(outdir / "adata_cellrank2.h5ad")
        outputs.append(str(outdir / "adata_cellrank2.h5ad"))

        main_fp = load_fp(outdir / "fate_probabilities.csv")
        sens_rows = sensitivity_analysis(adata, config, outdir, warnings, main_fp)
        if sens_rows:
            outputs.append(str(outdir / "n_states_sensitivity.csv"))

        cellrank_status = {
            "status": status,
            "route": route,
            "directionality_supported_by": directionality,
            "kernel_description": desc,
            "fate_probability_entropy_median": fate_entropy(main_fp),
            "warnings": warnings,
            **export_flags,
        }
        if directionality == "connectivity_only":
            cellrank_status["interpretation_limit"] = "Connectivity-only route has no intrinsic directionality; report as state-structure mapping unless external direction evidence is supplied."
        write_json(cellrank_status, outdir / "cellrank_status.json")
        outputs.append(str(outdir / "cellrank_status.json"))
        update_manifest(config, step, status if not warnings else "WARN", inputs=[str(input_h5ad)], outputs=outputs, warnings=warnings, route=route, directionality_supported_by=directionality, kernel_description=desc, export_flags=export_flags)
    except Exception as exc:
        fail_step(config, step, exc, warnings=warnings)
        raise


if __name__ == "__main__":
    main()
