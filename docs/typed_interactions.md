# Chemically typed interactions

The typed layer is opt-in and independent of CNI shell analysis. It wraps
PDBe Arpeggio 1.4.4 (source commit
`258855b8ba13447f2776b232ca32884d637c6a9c`), which implements CREDO geometry
and atom typing. Install Open Babel and the pinned optional Python dependency,
then run:

```text
conda install -c conda-forge openbabel=3 gemmi biopython
pip install -e '.[interactions]'
type-coordination-interactions structure.cif --output typed_output/
type-coordination-interactions structure_1.cif structure_2.cif --mapped-edges mapped_edges.csv --output typed_batch/
```

The legacy SSCNA command, distance contacts, cofactor exclusions, CSV schemas,
and existing fixture outputs are not changed. Arpeggio's native JSON is retained
alongside `atom_interaction_evidence.csv`, `direct_protein_protein_interactions.csv`,
`protein_cofactor_interactions.csv`, `water_endpoint_interactions.csv`,
`residue_pair_interactions.csv`, `per_structure_interaction_summary.csv`, and
`interaction_manifest.json`. The evidence table emits one row per supported
atom-pair type. A residue pair can have several rows/types. An external baseline
edge can be added to `summarize_residue_pairs(..., baseline_pairs=...)`; it is
reported as `proximity_only_unclassified` only when no typed evidence was found.

## Rules and parameters

The pinned Arpeggio defaults used here are: pH 7.4 for added hydrogens,
van-der-Waals compensation 0.1 Å, interacting-neighborhood cutoff 5.0 Å,
sequence-adjacent interactions excluded, and ambiguity definitions disabled.
The neighborhood cutoff is a candidate search parameter, not an interaction
label criterion. The pinned `arpeggio/core/config.py` criteria include:

| Type | Chemistry and geometry rule |
| --- | --- |
| Hydrogen bond | typed donor/acceptor; donor–acceptor ≤3.9 Å; donor–H–acceptor ≥90°; Arpeggio's polar fallback is retained in its JSON but not promoted to hydrogen bond |
| Weak hydrogen bond | typed weak donor/acceptor; ≤3.6 Å; H-angle ≥130°; donor substituent angle 30–150° |
| Ionic | oppositely ionisable atom types and ≤4.0 Å |
| Hydrophobic | both atoms typed as hydrophobes and ≤4.5 Å |
| Aromatic atom contact | both atoms typed aromatic and ≤4.0 Å |
| Aromatic stacking | aromatic ring planes; centroid distance ≤6.0 Å and Arpeggio plane geometry category |
| Cation–π / donor–π | atom-ring centroid distance ≤4.5 Å and ring-normal theta ≤30° |
| Methionine sulfur–π | methionine sulfur to aromatic centroid ≤6.0 Å |
| Metal coordination | typed acceptor and metal atom, ≤2.8 Å |
| Covalent / disulfide | Open Babel perceived bond; Cys SG–SG is separately labeled disulfide |

The manifest records parameters, software, output row counts, hydrogenation
assumptions, and altloc policy. Arpeggio adds hydrogens with Open Babel at pH
7.4 only if the input contains no H/D anywhere; a partial-H input suppresses
addition for the entire structure. Existing H coordinates are therefore a
user-controlled protonation assumption. Hydrogen-bond angles are measured from
the hydrogenated coordinate file when endpoint mapping is unambiguous;
unavailable angles are explicitly marked. The Arpeggio JSON export does not
include measured aromatic plane angles, so those rows retain the classified
plane category and state that the angle was unavailable. Arpeggio evaluates
observed alternate conformers without selecting one. It does not reconstruct
missing atoms or resolve occupancy alternatives.

For a mapped edge ledger, supply `--mapped-edges edges_by_structure.csv`. It
accepts one row per structure × edge with `structure_id`, `edge_id`, `group`,
`sequence_hash`, `src_resname`, `src_chain`, `src_resnum`, optional
`src_insertion_code`, corresponding `dst_*` fields, `baseline_contact`,
`eligible`, and optional JSON `eligible_by_type` (for example,
`{"hydrogen_bond": true, "disulfide": false}`). Run several structures in
one command to create `sequence_balanced_edge_type_summary.csv`; each structure
also receives its mapped edge observation table. Residue numbers and chains
are auth identifiers for that structure, so alignment positions must first be
mapped to observed residue identifiers.

For a compact reference analysis, Arpeggio can be restricted to mapped residue
endpoints with repeatable `--selection /A/272/` options. This limits its search
to their local neighborhoods; it does not change the geometric criteria.

## Separate graph layers and sequence-balanced summaries

Rows are partitioned as direct protein–protein, protein-cofactor/ligand, or
water endpoint evidence. Water-mediated paths are not folded into direct
protein–protein links. To summarize a mapped graph, pass one row per
structure × edge with `edge_id`, `group`, `sequence_hash`, `eligible`,
`eligible_by_type`, and `interaction_types` to
`sequence_balanced_type_summary`. If `eligible_by_type` is omitted, the runner
checks the mapped residues and their observed atom names, then records a
type-specific capability screen (for example, complete observed aromatic ring
atoms or charged-group atoms). Supplied `eligible_by_type` values override
those defaults. This denominator screen does not create interaction labels;
Arpeggio still decides positive calls using its atom typing and geometry. Check
the emitted eligibility table when the structure contains unusual residues or
protonation states. For each unique sequence and group, the positive vote follows the
manuscript rule: at least half of its eligible replicate structures are
positive. Output denominators count only eligible structures/sequences for
that interaction type. The standalone summary function consumes the eligibility
flags already present on its observation rows.

Typed interaction labels are structural evidence only. They do not establish
function, proton transfer, or reaction mechanism. Do not apply the manuscript's
18-edge support threshold to each chemical type without a separate analysis.
