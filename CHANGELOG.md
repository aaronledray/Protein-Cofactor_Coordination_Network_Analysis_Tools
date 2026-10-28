# Changelog

## Unreleased

### Removed (breaking)
- The legacy single-structure script (`1_Single_Structure_..._SSCNA_v0.0.2.py`), the legacy
  comparison scripts (`2_...STSQNC`, `3_...STMQNC`), and `modules/legacy_adapter.py`.
  The `Coord_Breakdown` report schema, sidecar configuration files, interactive prompts,
  static PNG plots, the `Residues_of_Interest` mode, and the `--mode`/`--cohesive-viewer`
  options are gone.

### Changed (breaking)
- `sscna analyze` is now implemented in `modules/analyze_cli.py` on top of
  `analyze_structure()`. It writes `Coordination_Residues.csv`, `Coordination_Atoms.csv`,
  `Coord_Links.csv`, `Coord_Contacts.csv`, and the cohesive viewer
  (`coordination_network.html`) with three shells by default; `--no-plots` writes the
  tables only and `--output-dir` sets the output directory.

### Added
- Opt-in Arpeggio typed-interaction layer (`modules/interaction_typing.py`,
  `type-coordination-interactions`) and an `--external-residue-network` overlay for the
  cohesive viewer.

## 0.2.1

### Added
- Opt-in `explicit_sites` for `analyze_structure()`: caller-defined cofactor sites
  (for example heme plus copper).
- Opt-in `exclude_hydrogens` for `analyze_structure()`.

### Fixed
- mmCIF files that omit the optional `_atom_site.pdbx_PDB_ins_code` column now load.

### Notes
- Hydrogens still seed shells by default (legacy behavior unchanged), so shell sizes
  differ between models with and without deposited hydrogens unless
  `exclude_hydrogens=True`. Consider making exclusion the default in a future
  labeling-contract version.

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

### License
- Relicensed from GPLv3 to the BSD Zero Clause License (0BSD) for this release.

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
