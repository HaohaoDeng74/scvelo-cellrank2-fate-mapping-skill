# API compatibility notes

Target environment: Python 3.10, scanpy >=1.9, scvelo >=0.2.5, cellrank >=2.0.

CellRank APIs can differ across versions. Scripts should:
- record package versions in `run_manifest.json`;
- prefer documented kernels and estimators;
- attempt `GPCCA.compute_macrostates()` when available;
- fall back to `fit()`/`predict_terminal_states()` only when necessary;
- write machine-readable outputs from `adata.obsm`, `adata.obs`, and `adata.uns` even when estimator serialization fails.

CellRank2 supports multiview kernels, including gene-expression/connectivity, RNA velocity, pseudotime, developmental potential, experimental time and metabolically labeled data. Select the kernel according to available data and biological assumptions.
