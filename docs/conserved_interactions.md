# Conserved interactions across structure sets

`modules.conserved_interactions` (experimental) asks which interactions recur in nearly every
sequence of every group of a structure set, and what chemistry they have. It generalizes the analysis
of the heme-copper oxidase comparison and uses the same rules, so results can be checked against it.

## Inputs

All functions take tidy tables; `load_atoms(path, chain, position_map)` builds the atom table from a
structure file.

| Table | Columns |
|---|---|
| `atoms` | `structure_id, position, seq_index, resname, atom, element, x, y, z` (heavy atoms; `position` is the shared label, for example an alignment column; `seq_index` is the residue's order in its own chain) |
| `structures` | `structure_id, group, sequence_id` |
| `features` | `structure_id, feature, atom, element, kind, x, y, z` with `kind` in `metal`, `oxygen` |
| `waters` | `structure_id, x, y, z` |

## Functions

| Function | What it does |
|---|---|
| `residue_contacts` | Residue-residue contacts (closest heavy atoms 1.9 to 4.0 Å, not sequence neighbors) with typed chemistry: closest-atom class (`bb_bb`, `bb_sc`, `sc_bb`, `sc_sc`), `polar` (N/O pair within 3.5 Å), `polar_sc`, `polar_bb_bb`, `hydrophobic`, `aromatic` (ring atom within 4 Å), `ionic` |
| `cofactor_contacts` | Residue contacts to named cofactor features: `coord` (metal to N/O/S within 2.8 Å), `polar`, `ionic`, `contact` |
| `water_bridges` | Residues bridged to an oxygen feature by one modeled water, and whether they also contact it directly |
| `relay_steps` | Fewest hydrogen-bond links from source atoms to each residue through waters and side-chain N/O atoms |
| `sequence_vote` | The conservation rule: one vote per sequence (positive when at least half of the sequence's structures are), then at least 90% of sequences in every group, with at least five per group; chemistry flags can be voted over contact-positive structures only (`condition`) |
| `conserved_contacts` | Contacts passing the vote with the conserved chemistry of each |
| `chemistry_summary` | Counts of conserved contacts by conserved chemistry class |
| `anchor_audit`, `misaligned_sequences` | Flag structures whose correspondence puts an anchor residue (for example the axial histidine) far from its site |

## The alignment audit

A residue correspondence that is wrong in one region makes every contact there look absent in the
affected sequences, which caps conservation below any support threshold. Run `anchor_audit` with a
residue that must lie next to a known site (the metal-binding histidine and the metal) before voting,
and decide whether to drop the flagged sequences. In the oxidase comparison eight of 63 sequences failed
this check.

## Protonation and hydrogen-bond geometry (optional tools)

`modules.protonation` wraps PROPKA (`pip install propka`, or `pip install .[protonation]`) and `reduce`
and adds geometry helpers that need neither: `run_propka`, `add_hydrogens`, `directional_hbonds` (H...A
2.5 Å or less, D-H...A 120° or more, D...A 3.5 Å or less), `virtual_histidine_hydrogens` (test a doubly
protonated histidine), and `carboxylate_distances` (a carboxylate cannot donate to a carboxylate unless
one is protonated, and an O...O contact of 2.7 Å or less is the geometry of a protonated pair).

PROPKA is empirical and built on heavy-atom structures. It does not know the redox, ligation, or charge
state of metal centers, so metal-bound histidines return large negative values and results depend on
whether cofactors are in the input. Report the sensitivity (cofactors and ions kept, hemes and copper
only, protein only) alongside any value.

## Validation

On the 253 structures (55 sequences) of the oxidase comparison, with positions from the alignment, the
functions reproduce the manuscript analysis exactly: the same 484 conserved contacts (no extra, none
missing), every chemistry class identical for every contact, the same conserved cofactor contacts
(iron, CuB, propionate), the same water bridges, and the same hydrogen-bond-link sets at one, two, and
three links. The check is `manuscript/controls/cni_validation.py` in the manuscript repository.

## Limits

Geometry is judged from heavy atoms; no hydrogens are used except in `modules.protonation`, so `polar`
means hydrogen-bond-capable. Waters are crystallographic and their orientation is unknown. Group
labels, sequence identity, and the residue correspondence are the caller's; the vote is only as good as
the correspondence.
