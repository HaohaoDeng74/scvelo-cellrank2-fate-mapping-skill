# v3 revision notes

This revision addresses the v2 re-review findings.

## Accepted and implemented

- Fixed `n_states_sensitivity` to use `adata_tmp` and its local estimator outputs rather than the original `adata`.
- Added estimator-first CellRank result export for fate probabilities, terminal states, initial states, macrostates and coarse transition matrix, with AnnData key guessing only as fallback.
- Added export flags to `cellrank_status.json`.
- Added `directionality_supported_by` and explicit interpretation limits for connectivity-only non-velocity routes.
- Added native CellRank `compute_lineage_drivers()` attempt from `gpcca_estimator.pkl`, with `correlation_fallback` as an explicit fallback method.
- Fixed `smoke_test.py` so `py_compile` writes bytecode to a temporary directory rather than `__pycache__` inside the skill tree.
- Fixed relative `output_dir` resolution so relative paths are resolved against the config file directory.
- Fixed scVelo plot output tracking by globbing actual files created by scVelo.
- Added `workflow_test_minimal.py`, which runs a minimal non-velocity workflow when `anndata`, `scanpy` and `cellrank` are installed, and skips cleanly otherwise.

## Still intentionally limited

- Full end-to-end RNA velocity workflow was not executed in this environment because `anndata`, `scanpy`, `scvelo` and `cellrank` are not installed.
- CellRank API objects vary across versions. The exporter uses public/common attributes first and falls back to AnnData storage discovery.
- Native CellRank driver export depends on successful estimator serialization and compatible CellRank version; otherwise results are explicitly labeled `correlation_fallback`.
