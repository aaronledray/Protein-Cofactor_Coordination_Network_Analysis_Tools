# Scope check 1: known answer on plastocyanin (1AG6)

Question: on a small, well-understood copper protein, does the tool recover the known site, and does it
report anything the structure or the literature contradicts? Structure: spinach plastocyanin, 1.6 A X-ray,
one Cu, 99 residues, 64 waters, no alternate locations (`reference_structures/0_Plastocyanin/1ag6.cif`).
Regression test: `tests/test_known_answer_plastocyanin.py` (expected values computed in the test from the
coordinates, not from tool output).

## Figures

`python -m modules.site_figures reference_structures/0_Plastocyanin/1ag6.cif CU` writes four figures and a checklist to
`site_check_output/1ag6/` (git-ignored): `1_distance_ladder`, `2_site_geometry`, `3_sulfur_environment`, `4_contact_map`. Each carries the
question a person should answer. The same set is produced for any structure with a metal cofactor, so it is the standard first look
at a site (module `modules/site_figures.py`, test `tests/test_site_figures.py`).

## Independent ground truth (Biopython and numpy only)

| Cu to | Distance (A) |
|---|---|
| His37 ND1 | 1.96 |
| His87 ND1 | 2.01 |
| Cys84 SG | 2.15 |
| Met92 SD | 2.88 |

Every other heavy atom is 2.86 A or more away, and all of those are carbons of the ligands themselves (His CE1 and
CG, Cys CB, Met CE) or Pro36 O at 3.47 A. Cys84 SG has two backbone amide nitrogens within 3.6 A: Asn38 N (3.45 A) and Ser85 N (3.58 A), and no water within 3.8 A.
With the amide hydrogens placed from the backbone, Asn38 N-H points at the sulfur (H...S 2.49 A, 158 degrees) and
Ser85 N-H does not (H...S 3.63 A, 79 degrees): Ser85 N is the next residue's amide, bonded through the peptide link,
and is a packing contact (see docs/scope_checks/03_sulfur_acceptor_survey.md: the next residue's amide is H-bond-like
for 1 of 269 metal-bound Cys sulfurs). Asn38 OD1 is 2.75 A from
Ser85 N.

## What the tool reports

| Layer | Result | Agrees |
|---|---|---|
| `analyze_structure` (defaults) | Primary shell His37, Cys84, His87, Met92 and Pro36 O; secondary shell Asn38 N and Ser85 N; distances equal the independent ones to 3 decimals | yes |
| `analyze_structure(direct_coordination=True)` | Direct coordination: His37, Cys84, His87. Met92 is a shell contact (2.88 A > 2.6 A cutoff) | three of four |
| `cofactor_contacts` (`coord` <= 2.8 A) | His37, Cys84, His87 | three of four |
| `cofactor_contacts` (`coordination_distance=3.0`) | adds Met92 | yes |
| `residue_contacts` | 263 residue pairs, identical to a brute-force recomputation (same pairs, distance difference 0.0) | yes |
| `residue_contacts` (38, 85) | polar side-chain contact, 2.75 A (the Asn38 to Ser85 clasp) | yes |

## Literature

- Four ligands His37, Cys84, His87, Met92 with a weak, variable Cu-S(Met92) bond. Reported Cu(II)-S(Met92) is 2.69 A in
  *Ulva pertusa* ([1IUZ](https://rcsb.org/structure/1iuz)); older poplar values are near 2.8 to 2.9 A. 2.88 A here is in that range.
- The Asn38 backbone NH donates a hydrogen bond to the Cys84 thiolate, and the Asn38 side chain forms a clasp with Ser85
  ([Asn38 mutants, Princeton](https://collaborate.princeton.edu/en/publications/h-bonding-maintains-the-active-site-of-type-1-copper-proteins-sit/)).
  The 1AG6 geometry shows the first (N...S 3.45 A) and one bond of the clasp (OD1...N 2.75 A).
- Not found in the sources searched: Ser85 N...S as a hydrogen bond, and the placed-hydrogen geometry above agrees that it is not one.

## Findings

1. The tool recovers the site with no tuning. Nothing it reports contradicts the structure or the literature.
2. "Coordination" is threshold-defined and the two modules disagree on the default: 2.8 A in `cofactor_contacts`, 2.6 A in
   `analyze_structure(direct_coordination=True)`. Both leave the axial Met92 (2.88 A) out, so a four-ligand site reads as
   three ligands plus a shell contact. Defaults unchanged (pinned); a log should print the cutoff next to every `coord` and
   keep a "weak axial" tier (for example 2.8 to 3.2 A) rather than a yes/no.
3. The 3.6 A primary-shell cutoff pulls in Pro36 O (3.47 A), a backbone carbonyl. It is labelled a shell contact, not a
   ligand, which is correct, but "primary shell" is not "ligand".
4. Sulfur is not typed as a hydrogen-bond acceptor. `polar` needs an N/O pair, so the Asn38 N-H...S(Cys84) bond, the
   textbook second-shell feature of a cysteine ligand, is not flagged in `residue_contacts` (the pair is typed bb_bb
   on its closest atoms, O...CA 3.34 A). The shell analysis does find Asn38 N as secondary shell. The same gap applies
   to Cys and Met sulfur in heme and Fe-S proteins. Ser85 N is a sequence neighbour of Cys84 and is excluded from
   `residue_contacts` by the separation rule, so only the cofactor-centred view sees it; geometrically it is a packing
   contact, not an N-H...S donor, so the exclusion does not hide an NH...S bond here (it can still hide other neighbour contacts, TODO D).
   Fixed for the Asn38 bond by the opt-in `sulfur_acceptor` flags (item C of the follow-ups).
5. 15 of 263 contacts lie within 0.05 A of the 4.0 A cutoff, so contact counts on one structure move with coordinate
   precision. Conservation across many structures averages this out; a single structure does not.

## Multi-species follow-up (done)

One structure cannot test the one-vote-per-sequence rule or conservation, so the same machinery was run on every X-ray
plastocyanin with a copper ion at 2.0 A or better: 30 structures, 12 sequences (8 eukaryotic, 4 cyanobacterial), with
each chain aligned to 1AG6 so that a position means the same residue everywhere
(`scripts/plastocyanin_vote_check.py`; full report in `docs/scope_checks/04_plastocyanin_vote.md`). Results:

- Both audits flag nothing: Cys84 lies within 8 A of the copper in all 30 structures, and none is backbone-only.
- The vote recovers the copper site without being told the site (the cofactor is the only input). This is a consistency check, not an independent test of the tiers: the tier cutoffs come from a survey of type 1 copper sites that includes these proteins, so recovering their ligands is expected. His37, Cys84 and His87 are conserved in the `bond` tier in all 12
  sequences and in both lineages.
- Met92 shows what a convention does to a vote. It is a ligand in every sequence (`bond` or `weak` in 12 of 12) but the
  2.7 A boundary splits it: `bond` in 17% of sequences, `weak` in 83%, and the single 2.8 A cutoff catches it in only 58%,
  so neither class alone passes. A vote on "bond or weak" (`ligated`) passes at 100%. For a log of conserved ligands, vote
  on `ligated`; the tier labels are descriptive.
- The Asn38 N-H...S(Cys84) bond is conserved as `polar_s`, `polar_s_bb` and `polar_s_geom` in all 12 sequences, and the
  existing `polar` flag never sees it (0% of sequences). This is the finding that item C of the follow-ups fixes.
- The Cys84-His87 backbone N...S contact is not conserved (`polar_s_bb` 25%, `polar_s_geom` 33%): a real contact in some
  sequences and not others, which is the kind of thing a single structure cannot reveal.
- Only 16 of 99 positions are identical in all 12 sequences, yet 168 of 450 candidate residue pairs are conserved contacts.
  In 97 of those 168 neither residue is invariant: packing that survives substitution, which a sequence-only analysis misses.
  In 18 both residues are invariant.
- Voted separately, eukaryotic and cyanobacterial plastocyanins keep 114 of the 168 contacts (`min_sequences` lowered to 4
  so the 4-sequence lineage can vote).

Not covered: this is one family, so it tests the machinery and the copper site, not whether the conserved contacts are the
ones that matter for design (the TEM-1 result stands: contact-derived directives are descriptive).
