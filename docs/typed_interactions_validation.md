# Typed interaction validation notes

This validation uses the official PDBe Arpeggio 1.4.4 executable (commit
`258855b8ba13447f2776b232ca32884d637c6a9c`; Python 3.10.14, Gemmi 0.7.5,
Open Babel 3.1.0). Native Arpeggio JSON is the reference. The adapter emits the
same supported target-type flags and retains its raw JSON, while splitting
protein, cofactor/ligand, and water endpoint records. It intentionally drops
non-target labels such as clashes, proximal, VDW, generic polar, and carbonyl
flags. There was no disagreement in target type calls in these runs; this is
an adapter round-trip check, not an independent reimplementation of CREDO.

## Curated structures

| Structure | Scope | Typed atom rows | Layer counts |
| --- | --- | ---: | --- |
| 1AG6 | Full compact plastocyanin mmCIF | 396 | 303 direct protein, 6 cofactor/ligand, 87 water endpoint |
| 1AR1 | Arpeggio neighborhood around all endpoints in the manuscript's 18-edge core | 400 | 230 direct protein, 150 cofactor/ligand, 20 water endpoint |

The 1AG6 adapter output contains 285 hydrophobic, 88 hydrogen-bond, 16
aromatic atom-contact, 3 covalent, 3 metal-coordination, and 1 ionic atom-pair
rows. Direct non-water hydrogen-bond angles were measured from Arpeggio's
prepared H coordinates when endpoint mapping was possible. Water hydrogen
bonds retain Arpeggio's water-specific rule, which does not require an explicit
donor-H angle.

For the 18 mapped 1AR1 residue pairs, an independent heavy-atom distance check
found all 18 at ≤4 Å (2.776–3.751 Å). Thus the baseline topology remains intact.
At this representative structure, typed direct-protein evidence appeared on
4/18 pairs: hydrophobic evidence on 272–326, 389–393, and 408–412; aromatic
atom-contact evidence on 272–326 and 276–325; 272–326 has both. The other
14 pairs are `proximity_only_unclassified`. No direct hydrogen bond, ionic,
metal, disulfide, cation–π, donor–π, sulfur–π, or ring-stacking label was
assigned to these 18 edges in this structure. These are structure-level calls,
not sequence-balanced support estimates.

## Unit and edge-case checks

`tests/test_interaction_typing.py` covers positive and near-miss geometry for
each emitted type, including distance limits, hydrogen-bond angles, aromatic
plane geometry, and covalent/disulfide chemistry. It also checks multiple
labels on one atom pair, altloc ambiguity reporting, absent typed evidence,
water/cofactor separation, atom order independence, and type-specific eligible
denominators with the manuscript's half-of-replicates sequence vote.

Arpeggio's JSON omits model identifiers. The adapter rejects multi-model
structures instead of silently combining model contacts; split models before
analysis. Arpeggio examines alternate conformers but omits the altloc ID from
contact JSON, so the adapter records the available altlocs as ambiguous when a
key has multiple conformers. Missing atoms are not built. If no H/D exists,
Open Babel adds hydrogens at the manifest pH; if even partial H/D exists,
addition is skipped for the whole structure. Histidine state therefore follows
the supplied H atoms or Open Babel's pH 7.4 behavior, rather than an explicit
tautomer-selection workflow. Mapping uses auth chain, residue number,
insertion code, residue name, and atom name, so it is independent of atom file
order; the 1AR1 edge check used chain A and its auth residue numbers.

The 1AR1 example confirms that the distance-defined 18-edge graph is much
denser than any one typed layer in a representative structure. A manuscript
update could show a per-edge interaction-type matrix and separate sequence-
balanced support plots by type and group. The latter requires the complete
structure-specific endpoint ledger and chemistry-specific eligibility flags;
do not infer those group-level results from this one representative or reuse
the 70% contact-network cutoff for typed interactions. No functional or proton-
transfer interpretation follows from these labels alone.
