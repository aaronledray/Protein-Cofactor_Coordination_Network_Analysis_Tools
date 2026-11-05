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

## State-stratified sequence vote (optional)

`sscna conserve --stratify-state` adds `state_contacts` while retaining the pooled `conserved_contacts` table. The manifest must have at least two
non-empty `state` labels (letters, digits and underscores) with the same meaning across all sequences and groups. For each state, the command applies the existing sequence vote
independently: within a sequence, at least half of eligible structures must show the contact; across sequences, each group still needs at least five
voting sequences and 90% support by default. `state_structure_votes.csv` records the eligible and positive structure counts behind each sequence-state vote,
including zero-coverage rows with a missing vote. `state_contrasts.csv` compares each unordered pair of states as passing, failing support, or lacking
coverage; these labels summarize the vote only and are not a test for a state effect. Existing opt-in sulfur-acceptor and sequence-neighbor rules also apply when requested. A structure
missing either position is not a negative. The output includes the contact candidates seen
in any state, so a state with no eligible sequences has zero voters and fails the minimum explicitly.

This is a descriptive within-state comparison, not a test that state causes a contact change. Structures are still equally weighted inside each
sequence-state and may be correlated. The first version covers residue-residue contacts only; it does not stratify cofactor contacts, waters, or
relay links. State-specific directives are not generated; `--directives` cannot be combined with `--stratify-state`. See
[`state_stratification_design.md`](state_stratification_design.md) for the method contract.

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

## Metal-bound waters and ligands

`load_atoms` keeps standard residues only, so a water or small molecule bound to a metal is not in the atom table or the vote. `load_metal_ligands(path)` reports them
per structure, descriptively: one row per (metal atom, non-protein atom within `max_distance`, default 2.8 A) with the metal and ligand residue, atom, element, distance and a
water flag; the metal's own residue (a porphyrin) is excluded, an atom with alternate locations keeps its first conformer, and a modified amino acid counts as a ligand.
It is never voted: in the 277 oxidase binuclear centres, 37% of heme irons and 42% of CuB have such an atom, in 16 different species (peroxide, CO, water, azide, peroxo,
cyanide, hydroxide, NO, ...; `docs/scope_checks/12_oxidase_metal_ligands.md`), and the species reflects what each experiment modelled, so a vote on it would largely record
how the structures were prepared, not conservation. For a catalytic metal in a single structure (the zinc-bound water of carbonic anhydrase) it is part of the site and belongs in the log beside the
residues; across a family it is state information, which the unbuilt state stratification can use.

## The log

`design_log.write_log(directory, tables, structures=..., title=..., parameters=...)` (experimental) writes the result of an analysis as a folder: one tidy CSV per table you give
it (for example `conserved_contacts(coverage=True)`, the cofactor vote with its coordination tiers, `load_metal_ligands`, the audits), `structures.csv` (structure, group,
sequence), `manifest.json` and `SUMMARY.md`. The manifest records the log format and version, the date, the repository commit, the dataset (structures, sequences, and per group the
sequences, structures and negative sequences a contact can take), the parameters you pass, and for each table its file, row count, columns (with types and the rule behind each
column from the legend that travels with the table) and the legends themselves. `SUMMARY.md` is the same for a person: the fixed statement of what the log means, the dataset and group
depth (a group with no tolerated negative sequence is called out), the tables, the headline counts, metal ligands by species, the parameters and the rules in force. The writer
refuses to replace an existing log unless asked and never touches other files in the folder. The tests check that nothing the log writes claims importance or prediction outside
its fixed statement. Composing the analysis functions into a pipeline (a command, single-protein mode) is not part of the writer.

**Directives (decided 2026-11-02; kinds renamed 2026-11-03).** Optional: `contact_directives(contacts, structures)` and `cofactor_directives(vote, structures, flags)` turn the vote tables into
`conserved_contact` and `conserved_cofactor_contact` directives, which `write_log(..., directives=...)` writes as `directives.csv` with the policy in the manifest and summary. A directive
only restates one row of a table with its evidence: the groups in which the contact (or cofactor class) is conserved, and for every group the voting sequences, the positives and the
spare negatives. Its `claim_level` is `descriptive`: it makes no claim about the effect of keeping or changing the position. Tiers, by fixed rules: `conserved_in_all` (every group passes the
support rule: the headline), `family_specific` (exactly one group passes, and every other group has voting sequences and shows the contact in under half of them) and
`conserved_in_subset` (everything else conserved somewhere). A group with no voting sequences for a pair is never counted as showing it absent, and the statement names it. The statement
is generated from a template; nothing is written by hand, a directive with another tier or claim level is refused, and a test scans the statements and the policy for importance wording. On
the oxidase set the 2,787 pairs conserved in some group give 490 conserved_in_all, 208 family_specific and 2,089 conserved_in_subset; on the 12 plastocyanins all 182 directives (168 contacts,
14 cofactor classes) are conserved_in_all.

Thresholds kept at 90% support and 0.5 for family-specific (2026-11-03, "for now"): the support rule matches the manuscript. On the oxidase set the all-groups count is 552 / 490 / 406 / 360 at 80 / 90 / 95 / 100% support, and the family-specific count is 181 / 208 / 277 at 0.25 / 0.5 / 0.75 (90% support); both are parameters of `contact_directives` and `cofactor_directives` and are recorded in the manifest.

## Command line: `sscna conserve` (experimental)

`sscna conserve --manifest structures.csv --output-dir log/` runs the vote on a set of structures and writes the log folder (`modules/conserve_cli.py`, a thin
layer over `conserved_contacts`, `cofactor_contacts`, `sequence_vote`, `load_metal_ligands` and `write_log`; it adds no analysis of its own and its tests check that
the log equals the API's tables).

- **Manifest CSV**, one row per structure: `structure_id`, `path` (relative paths are taken from the manifest's folder) and `sequence_id` (structures of one sequence
  share one vote); optional `group` (default `all`) and `chain` (default `A`; the first model is read).
- **Position map CSV** (optional), long format: `structure_id`, `residue_number`, optional `insertion_code`, `position`. Residues without a row are left out. Without a
  map the residue number is the position, which is right only when the numbering is already common to the set. Building the map (an alignment) stays the caller's.
- **Defaults** are the manuscript's rules: support 0.90, five voting sequences per group, 4.0 A. Written without options: `conserved_contacts`
  (with the depth of the vote, `coverage=True`), `audit_backbone_only` and `structures`.
- **Opt-in**: `--metals CU,FE` (cofactor features, one per metal element, which adds `cofactor_contacts`; two metals of one element in a structure are one feature, so use the
  API with a `features` table for sites such as CuA and CuB), `--coordination-tiers` (bond/weak/ligated flags), `--anchor POSITION` (the anchor audit against the first
  metal), `--metal-ligands`, `--sulfur-acceptor`, `--neighbor-contacts`, `--directives`. `--overwrite` replaces an existing log; without it an existing log is an error.
- A group with fewer sequences than `--min-sequences` is named in the manifest notes; if directives are requested and nothing is conserved, `directives.csv` is not
  written and the notes say so. Errors in the inputs are reported on stderr with exit code 2.
- It runs in one process and holds every structure's atoms in memory; the 249-structure oxidase set peaked near 2.5 GB.
- Not included: state stratification across a set of sequences, the water-bridge and relay tables, and wires across structure sets (use the API).

### Single-protein mode (`--single-protein`, decided 2026-11-03)

With one sequence the vote has nothing to count (one vote; the minimum of five voting sequences per group can never be met), so the structures of the protein
(entries, chains, states) are the evidence instead (`modules/single_protein.py`). They are not independent (the same crystal form, laboratory and refinement),
so **nothing in this mode is evidence of conservation**, and the log says so in its manifest and summary.

- **Input:** the same manifest, with exactly one `sequence_id` (several are refused: use the vote), an optional `state` column (for example `ox`, `red`; labels are
  letters, digits and underscores) and `group` ignored. One structure alone is allowed; the log then notes that persistence cannot be assessed.
  `--support`, `--min-sequences`, `--sulfur-acceptor`, `--neighbor-contacts` and `--anchor` are refused (they belong to the vote or the alignment).
- **`contact_persistence`:** one row per residue pair in contact (the `residue_contacts` rules) in at least one structure: `structures_present` (structures that map both
  positions; one that lacks a position is not a negative), `structures_in_contact`, `fraction_in_contact`, the closest heavy-atom distance in each structure as
  `closest_min` / `closest_median` / `closest_max` (a distance beyond cutoff + 2 A is recorded as that limit), `chemistry` (classes true in every structure where the pair
  is in contact) and `near_cutoff` (the median within `--near-margin`, default 0.3 A, of the cutoff). With two or more states, `present_<state>` and `in_contact_<state>`.
- **`cofactor_persistence`** (with `--metals`, optionally `--coordination-tiers`): per (feature, position), the structures in which each class holds out of the structures
  in which the residue is recorded, and the distance to the feature. There is no `near_cutoff` here (the classes use different cutoffs).
- **Why distance is reported:** within one sequence most of the variation between structures sits at the cutoff. In the plastocyanin set 12 to 26% of the candidate pairs of a
  sequence were in contact in some of its structures but not all (P00299, 8 structures: 72 of 307), the closest approach of those pairs had a median of 3.75 to 3.9 A in the three sequences with several structures (the cutoff is 4.0), and about a third of
  the pairs present in every structure lay within 0.3 A of the cutoff. A yes/no list mostly records where the cutoff falls.
- **Directives (opt-in, claim level `observed`):** `observed_contact` and `observed_cofactor_contact`, restating a row with its counts and distances. Tiers by fixed rules:
  `observed_in_all_structures` (in contact in every structure that shows both positions, and at least two), `observed_in_one_structure` (in contact in exactly one structure,
  whatever the number of structures that show both positions) and `observed_in_some_structures` (the rest). The wording scan of the tests covers this mode too.
- **Log:** the same folder with `"mode": "single_protein"` in the manifest and the states recorded in `structures.csv`; there is no `conserved_contacts` table.

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
