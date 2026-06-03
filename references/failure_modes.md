# Failure modes and fallback strategies

## Missing spliced/unspliced layers

Problem: RNA velocity cannot be computed.

Action:

- Stop scVelo branch.
- Use CellRank2 non-velocity kernels only if pseudotime, real time, or connectivity structure is meaningful.
- Label output as non-velocity fate mapping.

## Poor velocity confidence

Problem: Velocity vectors are noisy or inconsistent.

Action:

- Subset more narrowly.
- Recompute neighbors and UMAP.
- Try stochastic mode as exploratory fallback.
- Do not present terminal states as high-confidence.

## Incoherent phase portraits

Problem: Key genes do not support fitted kinetics.

Action:

- Exclude poor-fit genes from interpretation.
- Reduce causal language.
- Use marker, DEG, and pathway support to decide whether the trajectory remains interpretable.

## Biologically implausible terminal states

Problem: GPCCA predicts terminal states that do not match known biology.

Action:

- Review cluster annotations.
- Scan `n_states`.
- Manually define terminal states only if justified.
- Use caution in driver gene interpretation.

## Global immune mixture artifact

Problem: Velocity appears to connect unrelated immune lineages.

Action:

- Do not interpret global velocity.
- Subset by major lineage.
- Recompute neighbors and velocity within each subset.

## Batch- or stress-dominated trajectory

Problem: Latent time or fate probability aligns with batch, mitochondrial genes, ribosomal genes, dissociation stress, or ambient RNA rather than biology.

Action:

- Reassess QC and integration.
- Remove technical confounders if justified.
- Mark trajectory as not interpretable if artifact persists.
