# Single-structure consistency: CNI functions against `sscna analyze` (scope check 4)

Test: `tests/test_single_structure_consistency.py` (7 tests). Structures: spinach plastocyanin 1AG6 (Cu), sperm whale myoglobin 1A6M (heme,
oxygen-bound) and human carbonic anhydrase II 1CA2 (Zn). `analyze_structure` runs with three shells and `first_model_only`, as the command does;
the CNI side uses `load_atoms` and `cofactor_contacts` on the same files.

## Where they agree

- **The command is the API.** `sscna analyze --no-plots` writes the four tables with the columns, row counts and values that
  `analyze_structure` returns (the command's default of three shells against the API's two is the only difference, and it is deliberate).
  One format quirk is pinned: the `links` table keeps `distance_A` as a three-decimal string (`"2.010"`), so a CSV round trip reads `2.01`.
- **The primary shell has one stated rule.** A residue is in the primary shell exactly when one of its N, O or S atoms lies within 3.6 A of any heavy
  atom of the cofactor residue (carbon atoms of the cofactor count as sources; hydrogens do not; waters and ligands are included). This
  reproduces the analysis's primary shell exactly on all three structures (5, 17 and 4 residues). The heme case needed the carbon sources: four residues
  touch only a methyl carbon of the heme (backbone O of Lys42 and Tyr103, two waters).
- **Distances.** The coordinates that `load_atoms` produces give every cofactor-to-shell distance of the analysis to five decimals.
- **Coverage.** Every protein residue in the analysis's primary shell is recorded by `cofactor_contacts` with `contact` true and a minimum distance not larger
  than the analysis's (CNI's minimum includes carbon atoms).
- **Direct coordination.** At the same cutoff (2.6 and 2.8 A tested), the protein residues `analyze_structure` marks `direct_coordination` are exactly the
  residues `cofactor_contacts` marks `coord`, on all three structures.
- **Tiers.** `coordination_table.annotate_contacts(analysis contacts)` and `cofactor_contacts(coordination_model="by_pair")` give the same bond or weak tier per
  copper ligand of 1AG6.

## Where they differ (known, and stated so)

1. **Waters and ligands.** The analysis counts a metal-bound water or small ligand as a direct ligand: the water on the zinc of carbonic anhydrase (2.09 A) and the
   oxygen molecule on the heme iron of myoglobin (1.81 A). The CNI loader keeps standard residues only, so `cofactor_contacts` cannot see either. For a log of the residues
   that coordinate a catalytic metal this is a real gap: the directive would list the three histidines of carbonic anhydrase and not the water that is the substrate
   site. It is a scope decision (waters are not residues; `water_bridges` handles waters through a separate input), not a bug. Decided 2026-11-02: they are
   reported per structure by `conserved_interactions.load_metal_ligands` (descriptive, never voted); the test `test_the_metal_ligand_table_reproduces_the_non_protein_direct_ligands_of_the_analysis`
   shows that table equals the analysis's non-protein direct ligands on the zinc and heme structures.
2. **Cutoffs.** The command's direct-coordination cutoff is 2.6 A and CNI's `coord` is 2.8 A (both pinned). At 2.88 A the axial Met92 of plastocyanin is in neither and is
   named only by the weak tier.
3. **Features must be built by hand.** `cofactor_contacts` takes a `features` table; nothing builds it from a structure and a residue name. The test needs a helper
   (`features_for`), as does the manuscript validation. A single-structure CNI run therefore needs that step, which belongs to the open
   single-protein mode question.
4. **Primary shell is not "closest atoms".** For non-metal cofactors the shell follows the stated rule above, so it contains residues that touch only the cofactor's
   carbon atoms; CNI's `contact` class (any heavy-atom pair within 4.0 A) is wider still.

No discrepancy was found in the values where the rules overlap. Mutations that drop sulfur from `coord`, widen its cutoff by 0.3 A, or shift the loader's
coordinates by 0.01 A are all caught.
