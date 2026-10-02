# Changelog

## 0.2.0

Backward compatible with 0.1.0: the legacy single-structure command, its
defaults, output filenames, and CSV schemas are unchanged. Everything below is
opt-in.

### Added
- `sscna` command with `analyze` (identical to the legacy script), `compare`,
  `wire`, and `substrate` subcommands.
- Motif-aware coordination tables (`contacts`) with exact atom-pair evidence,
  direct metal-ligand annotation, arbitrary shell depth, per-site analysis,
  carbon-seed inclusion, and named cofactor-class cutoffs.
- `--cohesive-viewer`: the all-in-one interactive viewer and its contact/atom
  tables for the legacy command.
- Pairwise comparison, family profiles, geometry-aware alignment, a
  template-anchored conservation map, and a profile conservation viewer.
- Cofactor-to-target protein-wire analysis (shared relay, proton, electron,
  PCET, and residue-level redox modes) with an interactive viewer.
- Hypothetical substrate-point chains through the coordination network, as an
  API (`modules/substrate_seeds.py`) and `sscna substrate`.
- ML-facing residue labels (`modules/ml_export.py`): labeling contract 1.0
  with golden fixtures, per-site/per-model examples, one cofactor family per
  analysis, optional solvent exclusion, assembly/symmetry provenance, and
  equivalent-site grouping for train/test splits.
- Opt-in `include_model_id=True` for the residue table.
- Benchmarks and sensitivity reports under `benchmarks/`; CI for the test suite.
- `docs/stability.md` defines stable, legacy-compatibility, and experimental
  tiers; the stable core is import-isolated from plotting code and tested for it.

### Changed
- Legacy plots no longer open windows or browser tabs, or block, when stdin or
  stdout is not a terminal (previously a headless run could hang in
  `fig.show()`). Output files are unchanged; `--show` / `--no-show` override.

### Notes
- Analysis uses the deposited asymmetric unit; biological-assembly and
  crystal-symmetry transforms are not applied.
- Crystallographic waters are ordinary shell members (motif `water`).
- In pooled multi-model mode, residues repeated across models are represented
  by their closest atom; use `site_model_mode="per-model"` to keep all models.
- Wire and substrate chains are structural hypotheses, not transfer-rate,
  binding, or reactivity calculations.
- The tracked example output for 4ub6 under `SSCNA_output/` was removed from
  the repository; generated output belongs in ignored directories.
