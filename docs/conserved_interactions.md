# Conserved interactions across structure sets

`modules.conserved_interactions` (experimental) asks which interactions recur in nearly every
sequence of every group of a structure set, and what chemistry they have. It generalizes the analysis
of the heme-copper oxidase comparison and uses the same rules, so results can be checked against it.

## What the output means

The output is a description: which contacts, with what chemistry, recur in nearly every sequence of every group of the set you gave it. It does not
say that a conserved contact is important, that it stabilizes or enables anything, or that changing it would matter. Conservation can reflect shared
ancestry as well as constraint; a contact can be absent from the vote because the sequences differ, the alignment is off, or one small group is
strict (see "How stable the vote is"). The project's own TEM-1 test (reported with the design work, outside this repository, and not rerun here) was
negative for using contact-derived directives predictively, so any directive built from this output should be worded as a description of the
structures, with the group sizes and margins beside it, and should not claim a functional or design effect without separate evidence. The flags are
heavy-atom geometry (`polar` is a hydrogen-bond-capable pair, not a proven bond); where a hydrogen is placed from the backbone or ring the placement is a
model, not an observation.

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

## Coordination tiers by metal and donor (optional)

One global cutoff cannot suit every metal-donor pair (type 1 copper: Cu-N(His) about 2.0 A, Cu-S(Cys) about 2.2 A, and the axial Cu-S(Met) anywhere from 2.5 to 3.5 A across families). `modules/coordination_table.py` keys a `bond` and a `weak` cutoff on (metal element, donor class) and records the source of each row; the type 1 copper rows come from the PDB survey in `docs/scope_checks/02_cu_type1_survey.md`. The Cu-S(Met) tiers (bond up to 2.7 A, weak up to 3.3 A) are conventions: the survey has no natural break.

`cofactor_contacts(..., coordination_model="by_pair")` adds `coord_tier` (`bond`, `weak`, `default`, `none`), `coord_donor_class` and `coord_rule` (the cutoffs that decided it) to each residue-feature row. For `analyze_structure`, `coordination_table.annotate_contacts(result["contacts"])` adds the same columns to the contacts table. A pair that is not in the table follows the global cutoff and reads `default`, so unlisted metals behave as before and say so. Backbone amide N is not a donor class: it stays a contact and a shell member but is never tagged `bond` or `weak`. The `coord` and `direct_coordination` columns are not changed, and with the option off the output and the manuscript reproduction are identical. The table grows one metal at a time, when a scope check needs it.

## Sulfur as a hydrogen-bond acceptor (optional)

`polar` needs an N/O pair, so a backbone N-H pointing at a cysteine or methionine sulfur (the standard second-shell feature of a thiolate ligand) is never typed. `residue_contacts(..., sulfur_acceptor=True)` adds three columns. `polar_s`: a backbone N other than Pro within 3.7 A of a Cys SG or 3.6 A of a Met SD in another residue, or a side-chain Ser OG, Thr OG1, Tyr OH, Asn ND2, Gln NE2, Arg N, Lys NZ, Trp NE1 or either His ring N (the tautomer is not known from heavy atoms) within 3.6 A. `polar_s_bb`: backbone N donors only. `polar_s_geom`: the hydrogen placed from the heavy atoms (backbone amide, Trp NE1, Arg, Asn/Gln amide, and each His ring N as if protonated) must point at the sulfur, H...S <= 2.9 A and N-H...S >= 120 degrees; hydroxyl and Lys NZ hydrogens rotate and are not tested, and a backbone amide needs the preceding residue in the table to place its H. A disulfide SG (another SG within 2.3 A) is never an acceptor. A donor atom that a metal coordinates is deprotonated, so it is not a donor: pass `metals=load_metals(path)` and every N/O inside the `bond` tier of `coordination_table` of a metal (the global 2.8 A for pairs not in the table) is excluded. Without `metals` nothing is excluded, so pass it for any structure with metal ions; the metal-bound His ligand N atoms would otherwise read as donors.

The windows and the donor choice come from `docs/scope_checks/03_sulfur_acceptor_survey.md`: in 2,048 atomic-resolution structures with the amide H placed from the backbone, 90% of backbone-N contacts within 3.7 A of a metal-bound Cys sulfur point the N-H at it (53% for free Cys, 48% within 3.6 A of a free Met sulfur, 19% for a disulfide Cys). Side-chain donors are weaker: with a hydrogen placed on the planar groups, the share of contacts within the window (non-disulfide Cys) that point an N-H at the sulfur is 58% for Gln NE2, 40% for Asn ND2, 32% for Trp NE1 and for Arg N, and 22% for His (as if protonated; the tautomer that works is usually only one of the two), against 73% for backbone N; hydroxyl and Lys NZ hydrogens rotate and cannot be tested, and deposited Asn, Gln and His orientations can be flipped. So `polar_s` is a capable-pair flag and `polar_s_bb` is the flag the evidence supports. Like `polar`, a flag marks a capable heavy-atom pair, not a proven hydrogen bond. `polar`, `hydrophobic` and all other columns are unchanged, and with the option off the output and the manuscript reproduction are identical.

Sequence neighbors stay excluded by `min_separation`. This does not hide N-H...S bonds: with the amide hydrogens placed, the next residue's amide points its N-H at the sulfur in 1 of 269 metal-bound Cys, 3 of 925 free Cys and 18 of 271 free Met (`docs/scope_checks/03_sulfur_acceptor_survey.md`), because it is bonded through the peptide link. Other neighbor contacts can still be hidden (TODO item D).

## How stable the vote is

Two controls, both on the heme-copper oxidase set (54 sequences in groups A, B1 and C of 31, 18 and 5) and on 12 plastocyanins
(`docs/scope_checks/07_shuffle_control.md`, `08_leave_one_out.md`):

- **Shuffle:** permuting the position labels within each structure (geometry unchanged) leaves 0 conserved contacts in every
  replicate for both sets (10 and 3 replicates), and the plastocyanin copper ligands are conserved in 0 of 10. Permuting only inside blocks
  of 3 positions leaves a handful (at most 7 of 168, at most 4 of 490); blocks of 5 or more leave none. The vote needs the position labels
  right to within a couple of residues, and it does not manufacture conservation from structures with no correspondence.
- **One sequence at a time:** a contact survives the 90% rule with at most `floor(0.1 n)` negative sequences, so 3 of 31 in group A,
  1 of 18 in B1 and none of 5 in C. Group C has exactly the minimum of 5 voting sequences, so all 490 conserved contacts sit at the minimum
  and leaving out any one C sequence loses all 490 (none are lost if the minimum is relaxed for that group); equally, one scrambled or
  atypical C sequence removes the whole set (a shuffle replicate that scrambled one sequence left 6 of 490). The three C sequences that block the most
  contacts alone hold back 7 to 11 contacts each. In group A, 16 of the 490 have no spare negative and in B1 36. In the plastocyanins 41 of 168 have none.

So a conserved contact is as reliable as the smallest group is deep. `conserved_contacts(..., coverage=True)` (opt-in) adds that depth to every
contact: `spare_negatives_<group>` and `spare_voters_<group>` (negative when the group fails), `conserved_in` (the groups it passes, for example `A+B1`),
`failing_groups`, `held_back_by` (when exactly one group fails and leaving out one of its negative sequences would make it pass, those sequences) and
`one_negative_away` (the same ignoring the minimum number of voting sequences, which is what names the group C sequence a 4-of-5 contact lacks). The
conserved set and the default columns do not change. `vote_coverage` computes the same from the table of `sequence_votes`.

**Group C (decision, 2026-11-02).** On the oxidase set the headline stays the manuscript's rule: all three groups at 90% with at least five voting
sequences (490 contacts). Treatments of group C change the claim, not just the count: requiring only A and B1 gives 721 (231 of them not conserved in C),
letting C have one negative gives 523, and 642 contacts are conserved in C but not in A and B1 (`docs/scope_checks/11_group_options.md`). A log therefore carries the
490 as the headline and, beside each contact, the groups it is conserved in and the depth columns above, and says that group C conservation means conserved in each of
five sequences. Changing the headline would be a change to the manuscript's rule and has to be made there too.

## Rules travel with the results

Every table returned by `residue_contacts`, `cofactor_contacts`, `water_bridges`, `relay_steps` and `conserved_contacts` carries a legend in `frame.attrs["rules"]`: one row per reported column with the rule in words, the distance that decides it and whether it is a default or an opt-in (`rules_of(frame)` returns it as a table, `format_rules(...)` as Markdown for a log). The legend is built from the parameters of the call that made the table, so a non-default cutoff shows up in it, and a flag that was not computed (the sulfur flags with the option off) says so. `class_rules(table, **parameters)` builds the same legend without running an analysis. Tests put a pair just inside and just outside each documented cutoff and check that the flag changes there. The legend is per table, not per row, except where the rule varies by row: `coord_rule` names the table row behind each coordination tier. `attrs` do not always survive merges or file writes, so write the legend next to the table (`format_rules`, or the `checklist.md` of `site_figures`).

## Wires across structure sets

`modules.conserved_wires` (experimental) asks which *wires* recur, not only which residues are reached.
`relay_graph` builds a directed hydrogen-bond-capable graph for one structure (source atoms, modeled waters, and
residues by position; a residue passes the wire on only if it was reached through a side-chain N/O atom, a
residue reached through its backbone only is a terminal). `shortest_signatures` lists the shortest paths to a
target as tuples such as `("w", 404)` (source, water, residue 404); `disjoint_path_count` counts vertex-disjoint
paths (each residue and water used once); `conserved_wires` runs these over a structure set and votes
reach (within k links), each signature, and reports independent-path counts per target and group. Signatures
treat waters as interchangeable (`w`), since modeled waters have no shared identity across structures. Several
start atoms (for example four propionate oxygens) raise the independent-path count compared with a
single-atom source. In the oxidase comparison it shows that Asp404 is reached by two conserved wires (a direct
link through His403 and a one-water link), His326 by one conserved one-water wire, and His325 by alternative
routes none of which is conserved alone.

## The alignment audit

A residue correspondence that is wrong in one region makes every contact there look absent in the
affected sequences, which caps conservation below any support threshold. Run `anchor_audit` with a
residue that must lie next to a known site (the metal-binding histidine and the metal) before voting,
and decide whether to drop the flagged sequences. In the oxidase comparison eight of the 63 sequences of the
archived sample failed this check.

`backbone_only_audit` flags structures that are backbone-only models: more than half of the mapped residues
other than Gly and Ala have at most five distinct heavy atoms (the default). Such a structure cannot show a
side-chain contact, so each side-chain contact reads as absent for its sequence. In the oxidase comparison four
structures (5GPN, 5J4Z, 5J7Y, 8GZU) score 1.0 and every other structure 0.034 or less. Neither audit runs
inside `conserved_contacts`; call them before voting and drop the flagged structures or sequences.

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

On the 249 structures (54 sequences: A 31, B1 18, C 5) of the oxidase comparison, with positions from the
alignment, the functions reproduce the manuscript analysis exactly: the same 490 conserved contacts (no extra,
none missing), every chemistry class identical for every contact, the same conserved cofactor contacts
(iron, CuB, propionate), the same water bridges, and the same hydrogen-bond-link sets at one, two, and
three links. The check is `manuscript/controls/cni_validation.py` in the manuscript repository; the link-count sets of `conserved_wires` agree with those of `relay_steps` (one to three links and the six-link Trp272 set).
The earlier figures (253 structures, 55 sequences, 484 contacts) included the four backbone-only models.
The check builds its atom table from the manuscript's own residue map, not with `load_atoms`, so it does not test
the loader; the rules in the next section are covered by `tests/test_conserved_interactions.py`.

## Loading structures

`load_atoms(path, chain, position_map)` reads one chain of the first model (`model_index`) with Biopython. For
mmCIF files the chain and residue numbers are the author's (`auth_asym_id`, `auth_seq_id`, with
`pdbx_PDB_ins_code`), so `position_map` keys are `(author number, insertion code)`. Residues without a map
entry are skipped. Water, ligands, and any other hetero residue (including modified amino acids such as MSE)
are not loaded; cofactor features and waters come in through the `features` and `waters` tables. Hydrogens
are dropped. Of an atom with alternate locations the first conformer is kept (`A` when there is one, otherwise the lowest letter). `seq_index` follows the
residue numbering, so unobserved residues and skipped hetero residues count as sequence steps. Choosing
the chain, removing duplicate chains, and building the position map are the caller's.

## Limits

Geometry is judged from heavy atoms; hydrogens are placed only for the opt-in flags that say so (`polar_s_geom`, the neighbor `polar` class) and in `modules.protonation`, so `polar`
means hydrogen-bond-capable. Waters are crystallographic and their orientation is unknown. Group
labels, sequence identity, and the residue correspondence are the caller's; the vote is only as good as
the correspondence.

Residue pairs closer than `min_separation` (2) in sequence are not reported. This hides few hydrogen bonds: in 2,048
atomic-resolution structures a side-chain acceptor takes its hydrogen bond from the backbone N-H of the neighboring
residue (offset +-1) in 0.2% to 3.1% of residues by type (never at -1), and those bonds are 9.1% of all side-chain-acceptor /
backbone-N-H bonds, the rest being at offsets 2 to 4 that are listed (`docs/scope_checks/05_sequence_neighbors_survey.md`).
This is a statement about frequency, not importance: a side chain that hydrogen-bonds a neighbor, or two adjacent side
chains that pack, can stabilize a site, so `residue_contacts(..., neighbor_contacts=True)` (opt-in) adds the contacts
between residues one apart, judged per atom pair. A pair of atoms counts when it is at least 4 bonds apart through the
peptide bond (bonded 1-2, 1-3 and 1-4 atoms never count; the count is a guard, since the polar and ionic pairs are already
at least 4 bonds apart) and belongs to a class in `neighbor_chemistry` (default `polar` and `ionic`;
`NEIGHBOR_CHEMISTRY_ALL` adds `hydrophobic` and `aromatic`):

- `polar`: a hydrogen-bond candidate with at least one side-chain atom. A donor (backbone N except Pro, Ser OG, Thr OG1, Tyr OH, Asn ND2,
  Gln NE2, Arg N, Lys NZ, Trp NE1, a His ring N) to an acceptor (any O, a His ring N) within 3.5 A, and where the donor's
  hydrogen can be placed it must point at the acceptor (H...A <= 2.5 A, N-H...A >= 120 degrees). A first version counted any
  N/O pair and called 9.6% of neighbor pairs polar (60% of Asp followed by Pro, whose N has no hydrogen); the survey of real bonds
  gives 1.7%.
- `ionic`: Lys/Arg to Asp/Glu within 4.0 A.
- `hydrophobic`: a C/S pair within 4.0 A with both atoms side-chain atoms beyond C-beta.
- `aromatic`: a ring atom (Phe, Tyr, Trp, His) within 4.0 A of a side-chain atom beyond C-beta. A ring beside the next residue's
  backbone or C-beta is not counted: it holds for 43% of all neighbor pairs that include an aromatic residue, so it would
  mark every aromatic residue, while the class fires in about 9% of them.

A residue pair with a qualifying atom pair gets a row with `neighbor` True (the column exists only with the option),
`min_distance` and `closest` of its closest qualifying atom pair, the selected flags typed as usual, and the other flags
and the sulfur flags False. The vote treats a neighbor pair like any other (negative when both positions are present and
the pair has no such atom pair). Rates in 2,048 atomic-resolution structures are in `docs/scope_checks/06_neighbor_chemistry_survey.md`
and `.csv` (polar 1.7%, ionic 0.14%, hydrophobic 2.6%, aromatic 2.0% of neighbor pairs): hydrophobic and aromatic contacts
depend strongly on the residue-type pair (Ile or Val followed by Pro, 42 to 47%), so judge a conserved one against that
background (the plastocyanin report does). Of the 5,751 polar neighbor pairs, Thr, Ser and Asp as the first residue account for 45% and Thr or Ser as the second residue for 34%
(the highest pair rates are Thr then His, 15%, and Ser then His, 13%).

On the 12 plastocyanins, 17 neighbor pairs occur in some structure with all four classes selected and none is conserved (none
reaches half of the sequences; `docs/scope_checks/04_plastocyanin_vote.md`). An earlier, looser version of the polar class reported a conserved Asp42 OD1 to Glu43 N-H
bond; the placed hydrogen points away (H...O 3.01 A, 79 degrees), so it was only a close pair, and the rule was tightened. Pro36 O is 3.30 A from His37 ND1, which is bonded to the copper; it is not a
hydrogen bond. Use `sulfur_acceptor` and `metals` in `conserved_contacts` (opt-in) to vote the N-H...S flags,
and vote the copper tiers as `bond` or `weak` together: a boundary between them (Cu-S(Met) 2.7 A) is a convention and splits a
conserved ligand across sequences (`docs/scope_checks/04_plastocyanin_vote.md`).
