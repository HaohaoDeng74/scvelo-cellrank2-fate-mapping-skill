# QC checklist

## Input QC

- [ ] AnnData file loads successfully.
- [ ] Cell annotation column exists.
- [ ] Sample / condition / timepoint metadata exist or their absence is documented.
- [ ] `X_umap` exists or UMAP is recomputed.
- [ ] `spliced` and `unspliced` layers exist for velocity analysis.
- [ ] Selected subset has enough cells.
- [ ] Major batch effects have been assessed.

## Velocity QC

- [ ] `filter_and_normalize` completed.
- [ ] `moments` completed.
- [ ] `recover_dynamics` completed or fallback mode was documented.
- [ ] `velocity_graph` computed.
- [ ] `latent_time` computed.
- [ ] `velocity_confidence` computed.
- [ ] Velocity stream agrees with known biological direction.
- [ ] Latent time agrees with markers or sampling timepoints.
- [ ] Phase portraits are acceptable for top genes.
- [ ] Low-confidence regions are identified.

## CellRank QC

- [ ] Kernel type is documented.
- [ ] Kernel weights are documented.
- [ ] GPCCA n_states range is documented.
- [ ] Terminal states are biologically plausible.
- [ ] Initial states are biologically plausible.
- [ ] Fate probabilities are interpretable.
- [ ] Sensitivity analysis was performed.

## Driver gene QC

- [ ] Driver genes were computed per lineage.
- [ ] Driver genes were not interpreted as causal without validation.
- [ ] Top genes show dynamic gene trends.
- [ ] Top genes overlap with marker/DEG/pathway/literature evidence where possible.
- [ ] Phase portraits were inspected for key driver candidates.

## Reporting QC

- [ ] Each lineage has evidence grade A/B/C/D.
- [ ] Low-confidence conclusions are labeled as exploratory.
- [ ] Figures have captions.
- [ ] Methods include scVelo mode, CellRank kernel, GPCCA settings, and software versions.
