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

### Fixed
- `conserved_interactions.load_atoms` dropped an atom entirely when its highest-occupancy conformer was not the blank or `A` conformer, because Biopython reports the altloc of that conformer for a disordered atom. It now keeps the `A` (or only) conformer of every atom, as documented. Results on single-conformer atoms are unchanged; the heme-copper oxidase comparison, which does not use this loader, is unaffected.

- `conserved_interactions.load_atoms` set `seq_index` from the order of observed residues, so residues on either side of an unobserved stretch or a skipped hetero residue (for example a selenomethionine) counted as sequence neighbors and their contact was dropped. `seq_index` now follows the residue numbering. In the heme-copper oxidase set, 13 of 253 structures have such breaks.
- `conserved_interactions.load_atoms` dropped atoms that have no `A` or blank conformer (only `B`, `C`, ...). It now keeps the lowest-lettered conformer. One oxidase structure (7XMB) lost 14 residues (109 atoms) this way.

### Added
- `scripts/survey_cu_type1.py` (`fetch` then `summarize`): reproducible survey of type 1 copper ligand distances in high-resolution PDB structures, one value per protein; its report and per-protein table are in `docs/scope_checks/02_cu_type1_survey.md` and `.csv`. This is the evidence for the per-metal, per-donor-class coordination cutoffs.
- `scripts/survey_sulfur_acceptors.py` (`fetch` then `summarize`): N-H...S geometry around Cys SG and Met SD in 2,048 atomic-resolution (1.2 A or better) PDB structures, one per 30% sequence cluster, with the amide hydrogen placed from the backbone; with hydrogens also placed on the planar side-chain donors (Trp NE1, Asn ND2, Gln NE2, Arg N, and each His ring N as if protonated) and a His tautomer-sensitivity table; report in `docs/scope_checks/03_sulfur_acceptor_survey.md` (the per-contact table `03_sulfur_acceptor_survey_contacts.csv` is written by the script and not committed). Evidence for the sulfur-acceptor windows and donors.
- Claim audit (`docs/scope_checks/10_claim_audit.md`): documentation and caption wording corrected where it said more than the evidence supports; new section "What the output means" in `docs/conserved_interactions.md`.
- `tests/test_single_structure_consistency.py`: the CNI functions agree with `sscna analyze` where their rules overlap (three structures); known differences pinned (`docs/scope_checks/09_single_structure_consistency.md`).
- `scripts/shuffle_control.py` (position-label shuffle, local shuffles and scrambled sequences; one worker by default, one oxidase vote needs about 2.2 GB and a minute) and `scripts/leave_one_out.py` (leave-one-sequence-out and per-contact margins from the per-sequence votes, no geometry rerun): reports `docs/scope_checks/07_shuffle_control.md`, `08_leave_one_out.md`.
- `scripts/plastocyanin_vote_check.py` (`fetch`, `run`): multi-species plastocyanin check of the vote (30 structures, 12 sequences aligned to 1AG6 numbering), report in `docs/scope_checks/04_plastocyanin_vote.md`; `scripts/survey_sequence_neighbors.py`: side-chain acceptors of backbone N-H bonds by sequence offset (`05_sequence_neighbors_survey.md`); `scripts/survey_neighbor_chemistry.py`: background rates of the neighbor classes by residue-type pair (`06_neighbor_chemistry_survey.md`, `06_neighbor_chemistry_background.csv`).
- `residue_contacts(..., neighbor_contacts=True, neighbor_chemistry=...)` and `conserved_contacts(..., neighbor_contacts=True, neighbor_chemistry=...)` (experimental, opt-in): contacts between residues one apart in sequence, judged per atom pair (at least 4 bonds apart through the peptide bond), with a `neighbor` column. Classes: `polar` (a hydrogen bond with a side-chain atom; the placed hydrogen must point at the acceptor; a Pro N and metal-bound atoms never count), `ionic` (default) and `hydrophobic`, `aromatic` (side-chain atoms beyond C-beta only; a ring beside the neighbor's backbone is not counted). Default output and the manuscript reproduction are unchanged. Evidence and background rates in `docs/scope_checks/05_sequence_neighbors_survey.md` and `06_neighbor_chemistry_survey.md` (+ `.csv`); on 12 plastocyanins no neighbor contact is conserved.
- `sscna conserve --single-protein` (`modules/single_protein.py`, experimental): the structures of one protein as the evidence instead of the vote: `contact_persistence` and `cofactor_persistence` (structures in which each contact is seen, closest distances, a `near_cutoff` flag, optional per-state counts), `write_log(mode="single_protein")`, and opt-in `observed_contact` / `observed_cofactor_contact` directives with claim level `observed` and tiers `observed_in_all_structures` / `observed_in_some_structures` / `observed_in_one_structure`. Not evidence of conservation; refuses a manifest with several sequence ids. Default mode and the default log are unchanged except that the manifest now records `mode`.
- `sscna conserve` (`modules/conserve_cli.py`, experimental): writes the conserved-interaction log for a set of structures from a manifest CSV and an optional position-map CSV; manuscript defaults, everything else opt-in (cofactor metals, coordination tiers, anchor audit, metal ligands, sulfur acceptor, neighbor contacts, directives). A thin layer over the API; its tests check that the log equals the API's tables.
- `modules/design_log.py` (experimental): `write_log` writes an analysis as a folder (tidy CSV per table, `structures.csv`, `manifest.json` with provenance, parameters, columns and the legends, and a readable `SUMMARY.md`); refuses to overwrite an existing log unless asked. Optional descriptive directives (`contact_directives`, `cofactor_directives`, `write_log(..., directives=...)`): each restates one table row with its per-group evidence in three tiers (conserved_in_all, conserved_in_subset, family_specific), generated from templates, claim_level descriptive, no importance or effect claims. `docs/stability.md` now lists `coordination_table`, `site_figures` and `design_log`.
- `conserved_interactions.load_metal_ligands` (opt-in, descriptive, never voted): the non-protein atoms (waters, O2, CO, azide, peroxide, modified residues, other hetero residues) within 2.8 A of the metal atoms of one structure, with a legend; `scripts/survey_metal_ligands.py` and `docs/scope_checks/12_oxidase_metal_ligands.md` survey the oxidase binuclear centres (heme iron 37%, CuB 42% occupied, 16 species).
- `conserved_contacts(..., coverage=True)` and `vote_coverage` (opt-in): the depth of the vote behind every contact (`conserved_in`, `failing_groups`, `spare_negatives_<group>`, `spare_voters_<group>`, `held_back_by`, `one_negative_away`); the conserved set and default columns are unchanged. `scripts/group_options.py` and `docs/scope_checks/11_group_options.md` record what each treatment of a small group does to the heme-copper oxidase set (group C: 490 / 721 / 523).
- `conserved_interactions.sequence_votes` and `contact_observations` expose the two halves of the vote (the per-sequence votes, and the table of candidate pairs per structure); a behaviour-preserving extraction, the 490-contact reproduction is unchanged.
- `conserved_contacts(..., sulfur_acceptor=True, metals=...)` passes the opt-in sulfur flags through to the vote; default output is unchanged.
- `conserved_interactions.class_rules`, `rules_of` and `format_rules`: the rule and cutoff behind every reported column (`residue_contacts`, `cofactor_contacts`, `water_bridges`, `relay_steps`, the conserved-contact vote), built from the parameters actually used. Every table those functions return carries them in `frame.attrs["rules"]`, and `site_figures` writes them into `checklist.md`. The threshold literals in `residue_contacts` (3.5 A polar, 4.0 A hydrophobic/aromatic/ionic) became the constants `POLAR_DISTANCE` and `CHEMISTRY_DISTANCE` with unchanged values; tests put a pair just inside and just outside each documented cutoff.
- `residue_contacts(..., sulfur_acceptor=True)` (experimental, opt-in): adds `polar_s` (a backbone N within 3.7 A of a Cys SG or 3.6 A of a Met SD, or a listed side-chain N/OH within 3.6 A, in another residue), `polar_s_bb` (backbone N only) and `polar_s_geom` (the hydrogen placed from the heavy atoms for the backbone amide, Trp NE1, Arg, Asn/Gln amide and a His ring N as if protonated, pointing at the sulfur: H...S <= 2.9 A, N-H...S >= 120 degrees; hydroxyl and Lys NZ are not tested), from the survey in `docs/scope_checks/03_sulfur_acceptor_survey.md`. Hydrogen placement helpers are in `modules/protonation.py`. A disulfide SG is never an acceptor; a donor that a metal coordinates (its `bond` tier in `coordination_table`, else 2.8 A) is deprotonated and excluded, via the optional `metals` table (`conserved_interactions.load_metals`). Either His ring N is a capable donor otherwise, since the tautomer is unknown from heavy atoms. `polar` and every other column are unchanged, and with the option off the output and the manuscript reproduction are identical. Sequence neighbors remain excluded.
- `modules/coordination_table.py` (experimental, opt-in): coordination cutoffs by metal and donor class with a `bond` and a `weak` tier and a source on every row; only type 1 copper (His N, Cys S, Met S with a 2.7 / 3.3 A convention) so far. `cofactor_contacts(..., coordination_model="by_pair")` and `coordination_table.annotate_contacts(analyze_structure(...)["contacts"])` add `coord_tier`, `coord_donor_class` and `coord_rule`; unlisted pairs follow the global cutoff and read `default`. Defaults, the existing columns and the manuscript reproduction are unchanged. `site_figures` now draws its tiers from the table.
- `modules/site_figures.py` (experimental, opt-in): four standard human-check figures for one metal cofactor site (distance ladder against every cutoff, site geometry, ligand-sulfur environment, residue contact map) plus a checklist, via `python -m modules.site_figures STRUCTURE RESNAME`. Output goes to the git-ignored `site_check_output/`.
- `conserved_interactions.backbone_only_audit` (experimental) flags backbone-only models (more than half of the non-Gly/Ala residues with at most five heavy atoms) before voting, alongside `anchor_audit`.
- Conserved-interaction analysis (`modules/conserved_interactions.py`, experimental): typed residue-residue contacts, cofactor-feature contacts, single-water bridges, hydrogen-bond link counts through waters and side chains, the one-vote-per-sequence conservation rule, and an anchor audit for misaligned correspondences; reproduces the heme-copper oxidase comparison exactly (`docs/conserved_interactions.md`).
- Water-aware wires across structure sets (`modules/conserved_wires.py`, experimental): relay graph per structure, shortest-path signatures through residues and waters, vertex-disjoint path counts, and the one-vote-per-sequence vote of reach, signatures and independent routes.
- Optional protonation helpers (`modules/protonation.py`, experimental): PROPKA pKa, `reduce` hydrogens, directional hydrogen-bond geometry, carboxylate O...O distances; `pip install .[protonation]` installs PROPKA.
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
