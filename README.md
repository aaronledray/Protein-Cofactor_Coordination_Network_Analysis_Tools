# Coordination Network Identifier

[![Tests](https://github.com/aaronledray/Protein-Cofactor_Coordination_Network_Analysis_Tools/actions/workflows/tests.yml/badge.svg)](https://github.com/aaronledray/Protein-Cofactor_Coordination_Network_Analysis_Tools/actions/workflows/tests.yml)

Coordination Network Identifier finds protein residues and atoms around a
cofactor, labels primary and secondary coordination shells, and compares
coordination networks across structures. It is designed for reproducible,
headless use from Python or the command line.

## What it provides

- Single-structure coordination analysis for PDB and mmCIF files.
- Legacy CSV/PNG/HTML reports for interactive inspection.
- Side-effect-free pandas tables for downstream workflows.
- Batch processing with multiprocessing and per-structure error isolation.
- Optional per-site analysis, arbitrary shell depth, direct metal-ligand
  annotations, motif-aware atom-pair contacts, carbon-seed inclusion, and
  class-specific distance cutoffs.
- Chain deduplication through `modules/deduplicate_chains.py`.

By default, PCS and SCS seed selection excludes carbon atoms, chooses at most
one atom per chemical moiety except for the configured multi-atom moieties,
and uses a 3.6 Å cutoff. Existing legacy runs retain their current defaults
and output filenames.

## Installation

Python 3.9 or newer is supported.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

## Single-structure analysis

The original SSCNA command remains available:

```bash
python 1_Single_Structure_Cofactor_Network_Analysis_SSCNA_v0.0.2.py \
  --template reference_structures/0_Plastocyanin/1ag6.cif \
  --cofactor CU \
  --distance 3.6 \
  --exclude-moieties alanine_sidechain \
  --mode Coord_Network
```

The stable installed command is equivalent and uses the same legacy defaults,
configuration files, output names, and report schema:

```bash
sscna analyze \
  --input reference_structures/0_Plastocyanin/1ag6.cif \
  --cofactor CU \
  --distance 3.6 \
  --exclude-moieties alanine_sidechain \
  --mode Coord_Network
```

All legacy analyzer options can follow `sscna analyze`; `--template` is also
accepted as a compatibility alias for `--input`. The comparison workflows are
available as `sscna compare pairwise ...` and `sscna compare profile ...`;
`sscna wire` and `sscna substrate` are described below.

The legacy plotting path opens matplotlib windows and a browser tab only when
run from a terminal; piped, CI, cron, and subprocess runs write the same files
without opening or blocking on anything. `--show` and `--no-show` override
that choice.

Outputs are written under `SSCNA_output/`:

- `<structure>_Coord_Breakdown.csv`: residue summary plus atom coordinates.
- `<structure>_Coord_Breakdown_atoms.csv`: atom-level category report.
- `<structure>_Coord_Links.csv`: nearest cofactor→PCS and PCS→SCS links.
- PNG plots and `<structure>_1_template_coordination_network.html`, the
  original interactive viewer, unless `--no-plots` is supplied.

With `--cohesive-viewer` (opt-in), the HTML is instead written as
`<structure>_coordination_network.html`, together with
`<structure>_Coord_Contacts.csv` and `<structure>_Coordination_Atoms.csv`
(shell membership with separate `coordination_role` values for primary
coordinators, network context, and geometric active-site components). The
tidy modes (`--shells N > 2`, `--per-site`, `--direct-coordination`, class
cutoffs) always use the cohesive viewer and these tables.

The cohesive viewer is the all-in-one interactive viewer with full-protein context, focused shell views, motif coloring, and exact
primary/secondary/tertiary contact overlays. Its layer controls independently
show network atoms, network-scoped motifs/backbone, active-site context, and
broader focused-site motifs/backbone. View presets provide quick Network,
Network + motifs, Structural context, and Full protein context layouts.
Cofactor bonds remain visible in black across bond modes; protein/context
bonds use the selected structural representation. Contact overlays use
dotted grey lines with progressively thinner lines for deeper shells. The
viewer also includes axis units in Å, a compact contact legend, click-to-view
atom details, and an image-export toolbar action. All viewer controls are
arranged in a vertical toolbar on the left side of the HTML window.
Recognized motif-level interactions, including heme-propionate O contacts
with nearby polar N/O/S atoms, are included as primary contacts with
explanatory hover text. These remain distinct from direct Fe coordination.

Useful additive options include:

```text
--no-plots                         Skip PNG/HTML rendering.
--cohesive-viewer                  Use the all-in-one viewer and write contact/atom tables.
--show / --no-show                Force or suppress figure windows and browser tabs.
--compact-html                     Use a CDN-backed Plotly bundle for smaller HTML output.
--first-model                      Analyze only the first model.
--shells N                         Add TCS and deeper shells for N > 2.
--per-site                         Separate cofactor sites in tidy output.
--site-model-mode per-model        Keep otherwise identical sites separate by model.
--include-carbon-seeds             Allow carbon atoms to seed shells.
--direct-coordination              Annotate direct metal-ligand links.
--direct-coordination-cutoff 2.6   Set the direct-link distance in Å.
--cofactor-class-cutoff metal=2.8  Override a class cutoff; repeatable.
--cofactor-class-config rules.yaml Load named cofactor-family cutoff rules.
--verbose                          Show progress logs.
```

Named cofactor-class rules are additive and do not change the legacy 3.6 Å
default unless explicitly configured. Built-in families include `metal_ion`,
`heme`, `iron_sulfur_cluster`, `metallo_cluster`, and `organic_cofactor`.
The aliases `metal`, `iron-sulfur`, `cluster`, and `organic` are accepted for
backward-compatible command lines. A simple override remains valid:

```bash
sscna analyze --input structure.pdb --cofactor HEM \
  --cofactor-class-cutoff heme=3.3
```

For reusable or custom families, pass a YAML or JSON file:

```yaml
cofactor_classes:
  heme:
    residues: [HEM]
    cutoff_A: 3.3
  custom_cluster:
    residues: [ABC, ABD]
    cutoff_A: 2.9
```

The same `--cofactor-class-config` and `--cofactor-class-cutoff` options are
available on `sscna compare`. If a mixed analysis contains both configured
and unconfigured families, the legacy fallback remains included so the
analysis does not silently discard network atoms.

## Cofactor-to-target protein wires

Protein-wire analysis is an opt-in graph mode, separate from PCS/SCS/TCS shell
expansion. It builds one shared chemical relay network from cofactor atoms and
chemically plausible relay atoms, then ranks paths to explicit target residues
or atoms under optional electron-transfer, proton-transfer, or PCET
interpretations. Each node and hop retains its residue, motif, atom identity,
distance, interaction class, and capability/support annotations.

For example, to inspect KatG relay candidates from heme to two tryptophans:

```bash
sscna wire \
  --input katg_1sj2.pdb \
  --cofactor HEM \
  --target TRP:A:107:NE1 \
  --target TRP:A:91:NE1 \
  --distance 3.6
```

Targets use `RESNAME:CHAIN:RESNUM[:ATOM]`. Results are written to
`protein_wire_output/` as `wire_nodes.csv`, `wire_edges.csv`,
`wire_paths.csv`, `wire_targets.csv`, `wire_summary.json`, and an interactive
`wire_network.html`. Water relays are included by default; use `--no-water`
to exclude them. The default `--mode relay` exposes the shared network;
`--mode proton` favors polar/water-mediated hops, `--mode electron` favors
aromatic/redox and sulfur hops, and `--mode pcet` favors hops supported by
both interpretations; `--mode redox` enables the residue-level redox-relay
model described below. Polar contacts are
marked as inferred when hydrogens are absent; explicit hydrogen geometry is
validated when deposited hydrogens support it. Aromatic edges retain
face-to-face or edge-to-face ring-orientation evidence when ring atoms are
available. The importable equivalent is
`modules.protein_wires.analyze_protein_wires()`.

For a residue-level electron-hole relay hypothesis, use the opt-in redox
mode. It restricts relay nodes to cofactor atoms and redox-capable residues,
allows aromatic contacts out to 6.0 Å by default, prevents a path from
revisiting the same residue, and records unresolved target labels explicitly:

```bash
sscna wire \
  --input katg_1sj2.pdb \
  --cofactor HEM \
  --target TRP:A:91:NE1 \
  --target TRP:A:107:NE1 \
  --target TRP:A:321:NE1 \
  --mode redox
```

The redox output remains a structural candidate ranking, not a transfer-rate
or quantum-coupling calculation. `wire_targets.csv` and the summary JSON
report residue-name mismatches such as requesting `TRP:A:93:NE1` when the
coordinate file contains a different residue at that number.

The legacy command currently scans all structure models unless
`--first-model` is supplied. Alternate locations are not assigned a custom
policy; Biopython's normal disordered-atom selection behavior is used.

For `--shells N` with `N > 2`, the command writes tidy
`Coordination_Residues.csv`, `Coordination_Atoms.csv`, and `Coord_Links.csv`
tables plus `Coord_Contacts.csv`; the cohesive viewer is also available when
plots are enabled.

## Importable API

```python
from modules.coordination_api import analyze_structure

tables = analyze_structure(
    "reference_structures/0_Plastocyanin/1ag6.cif",
    "CU",
    exclude_moieties=["alanine_sidechain"],
)

residues = tables["residues"]
links = tables["links"]
```

`analyze_structure()` creates no files and returns `residues`, `atoms`,
`links`, and `contacts` pandas DataFrames. The atom table includes a chemical
motif label. The `contacts` table preserves every qualifying atom pair between
adjacent shells, including source/target motifs, exact atom names, distance,
model identity, and direct-coordination status. The residue table has one row per residue per site
and shell (pass `include_model_id=True` for an opt-in `model_id` column and one row per model;
the default schema is unchanged), including structure ID, cofactor identity, residue identity,
atoms involved, and minimum distance to the previous shell. Insertion codes
and hetero flags are retained in the tidy tables.

Motif labels come from the shared residue-specific registry in
`modules/motif_registry.py`, rather than from viewer-only annotations. It
provides explicit atom membership for histidine imidazole rings, ASP/GLU
carboxylates (`COO`), cysteine thiols, heme iron/porphyrin and propionate
oxygen atoms, Fe-S clusters, the FeMo cofactor, and the oxygen-evolving
cluster. The registry is also directly inspectable:

```python
from modules.motif_registry import motif_atoms, motif_for_atom

motif_for_atom("HIE", "NE2")             # "imidazole"
motif_for_atom("GLU", "OE1")             # "COO"
motif_atoms("HEM", "heme_propionate")    # O1A/O2A/O1D/O2D
```

The resolver preserves exact atom identity in every contact, so chemically
equivalent alternatives such as OE1/OE2 remain distinguishable in the tidy
tables and comparison signatures.

Use `site_mode="per-site"` to separate cofactor copies. Site identity includes
cofactor residue name, residue number, chain, insertion code, and hetero flag.
By default, identical identities across models are pooled to preserve the
legacy all-model behavior. Set `site_model_mode="per-model"` (or use
`--site-model-mode per-model`) when each model must define separate sites;
the resulting IDs include the model number. In pooled mode a residue that
repeats across models is represented by its single closest atom per moiety (to
the pooled cofactor), so later models' copies are not re-seeded; this is
intentional legacy behavior, and `site_model_mode="per-model"` is the way to
keep every model's residues. Alternate locations follow
Biopython's selected conformer behavior and are not treated as separate
cofactor sites. With `combinatorial=True`, selected cofactor residues within
the combinatorial cutoff are clustered into one site, but proximity is only
compared within the same model. This is useful for multi-cofactor systems
such as nitrogenase without allowing cross-model distances to create a false
cluster.

Legacy CSV reports can be adapted into the same four-table contract used by
the current API without rewriting or modifying the original files:

```python
from modules.legacy_adapter import legacy_csvs_to_analysis, legacy_csvs_to_signature

tables = legacy_csvs_to_analysis(
    "legacy/1ag6_Coord_Breakdown.csv",
    "legacy/1ag6_Coord_Links.csv",
    structure_id="1ag6",
)
signature = legacy_csvs_to_signature(
    "legacy/1ag6_Coord_Breakdown.csv",
    "legacy/1ag6_Coord_Links.csv",
    structure_id="1ag6",
)
```

The adapter understands the mixed summary/atom breakdown format and the
separate atom-level report. Legacy `cofactor->pcs` links become canonical
primary-shell edges with `legacy_primary_contact` provenance; they are not
silently relabeled as validated direct coordination. Legacy reports lack
model IDs, insertion codes, hetero flags, and element fields, so those values
remain empty rather than being guessed.

## ML-facing residue labels

`modules/ml_export.py` reshapes `analyze_structure()` tables into one
deterministic label row per structure, site, model, and residue. It is an
opt-in adapter and does not alter any legacy output:

```python
from modules.coordination_api import analyze_structure
from modules.ml_export import export_residue_labels, labeling_contract

params = dict(shells=3, direct_coordination=True, site_mode="per-site")
tables = analyze_structure("structure.pdb", "CU", **params)
labels = export_residue_labels(tables, model_id=0, site_id="site_1")
contract = labeling_contract(**params)  # store beside the labels
```

Each residue carries its shallowest shell, all shells it occupies, primary
and direct-coordination flags, motifs, coordinating atoms, minimum distance,
contact count, and the labeling-contract version (`LABELING_CONTRACT_VERSION`).
`contacts` is the full atom-pair representation; `links` is the nearest-atom
edge list retained for legacy compatibility. Pooled multi-model sites are
rejected unless `model_id` is given, so shell labels are never attributed to
models ambiguously. Rows sort stably regardless of input order.

### Labeling contract 1.0

`analyze_labeling_example()` runs the analysis and export under a frozen
contract and returns `(labels, contract)`:

```python
from modules.ml_export import analyze_labeling_example

labels, contract = analyze_labeling_example("1a6m.pdb", "HEM")
```

Defaults: 3.6 Å cutoff, 3 shells, no carbon seeds, direct coordination at
2.6 Å, no residue expansion, crystallographic waters kept as ordinary shell members (motif `water`; there is no exclusion option), and
`per-site` boundaries with `per-model` separation, so every example is one
explicit site and model. Labels include water rows by default; pass
`exclude_solvent=True` to drop HOH/WAT/DOD/H2O rows from the export (waters
still take part in shell propagation, so water-mediated residues keep their
shell labels, and the contract records the flag). `site_mode`, `first_model_only`, and the
combinatorial cutoff cannot be overridden; other overrides are allowed and
recorded in the returned contract with the effective cutoff and cofactor
family. One cofactor family per analysis is required: mixed families would
share a single cutoff (the largest configured), so they raise
`MixedCofactorError` unless `allow_mixed_cofactors=True`, which is recorded as
`mixed_cofactor_classes`. Combinatorial clusters of one family, such as
nitrogenase ICS+CLF+HCA, are a single family. Golden fixtures in
`tests/fixtures/labels/` pin the 1.0 output; a deliberate change requires
bumping `LABELING_CONTRACT_VERSION`. Distances are rounded to 4 decimals.

### Assemblies, symmetry copies, and split grouping

Labels come from the deposited asymmetric unit only. Biological-assembly
transforms (`BIOMT`, `_pdbx_struct_assembly`) and crystal-symmetry mates are
never applied, so contacts across symmetry-related chains are not seen. The
contract records this as `symmetry_context` (`assembly_policy`, space group,
assembly count, and `has_unapplied_transforms` when non-identity operators
exist). Copies physically present in the file are separate sites and models.

`modules.assembly_policy.equivalent_site_groups(labels)` assigns an
`equivalence_group` so related sites can stay on one side of a train/test
split. The default `coordination` level groups sites with the same primary
coordinator residue names; it is intentionally coarse and also merges distinct
sites with the same chemistry (the six SF4 sites of 2zvs form one group).
`level="full"` signs every shell and is only suitable for near-exact
duplicates, because real non-crystallographic copies differ by contact noise.
Grouping never proves symmetry relatedness; use sequence or chain identity
for that.

## Hypothetical substrate-point seeds

`modules/substrate_seeds.py` is an opt-in, side-effect-free API for tracing
chains from a modelled substrate position (for example a hypothetical O2 or
H2O2 site that is not in the deposited structure) through the coordination
network to the cofactor:

```python
from modules.substrate_seeds import analyze_substrate_seed

result = analyze_substrate_seed(
    "1a6m.pdb", "HEM",
    {"cofactor_atom": "FE", "offset": (1.775, -0.841, 0.378)},  # or (x, y, z)
    contact_cutoff=4.0,
)
result["substrate_contacts"]  # network atoms within the cutoff of the point
result["chains"]              # ranked chains: entry residue, hops, length, shells
result["steps"]               # one row per node: point -> network atoms -> cofactor
```

The point is either explicit coordinates or an offset from a named cofactor
atom, resolved separately within each cofactor site (per-site analysis is the
default here). Contacts are limited to atoms already in the shell network.
For each contacted residue the best chain (fewest hops, then shortest summed
distance) follows adjacent-shell contacts to a cofactor atom, and chains are
ranked per site. On oxy-myoglobin 1a6m, a point 2 Å from Fe on the distal side
yields the chain point → distal His64 NE2 (SCS) → bound O2 (PCS) → Fe. Chains
are structural hypotheses about contact connectivity, not binding, reactivity,
or energetic predictions.

The same analysis is available from the command line, with an interactive 3D
view (network atoms by shell, the point, dotted point contacts, and
legend-toggled chains; the top three are visible initially):

```bash
sscna substrate --input 1a6m.pdb --cofactor HEM \
  --from-atom FE --offset 1.775 -0.841 0.378      # or --point X Y Z
```

It writes `substrate_seed.csv`, `substrate_contacts.csv`,
`substrate_chains.csv`, `substrate_steps.csv`, `substrate_summary.json`, and
`substrate_chains.html` under `substrate_seed_output/` (`--no-html` and
`--compact-html` are supported). Wire-graph integration is not available yet.

## Batch analysis

```bash
python batch_coordination_network.py \
  --input reference_structures/0_Plastocyanin \
  --cofactor CU \
  --workers 4 \
  --output-dir coordination_batch_output
```

The batch command writes combined `coordination_residues.csv`,
`coordination_atoms.csv`, `coordination_links.csv`,
`coordination_contacts.csv`, and
`coordination_errors.csv`. A malformed or missing structure is recorded in
the error table without stopping other structures.

## Network comparison and family profiles

The analysis tables can be converted into a canonical, JSON-friendly network
signature for comparison or profile building:

```python
from modules.coordination_api import analyze_structure
from modules.network_comparison import (
    build_network_profile,
    build_network_signature,
    compare_network_signatures,
    score_signature_against_profile,
)

reference_tables = analyze_structure("reference.pdb", "CU", shells=3)
query_tables = analyze_structure("query.pdb", "CU", shells=3)
reference = build_network_signature(reference_tables, structure_id="reference")
query = build_network_signature(query_tables, structure_id="query")

pairwise = compare_network_signatures(reference, query)
profile = build_network_profile([reference])
profile_score = score_signature_against_profile(query, profile)
```

Signatures retain primary, secondary, and tertiary atom-level contacts,
motif and atom identity, residue features, distances, cofactor identity, and
explainable unmatched features. Raw residue numbering is intentionally not a
comparison key by default, so homologous structures can be compared without
pre-alignment; pass a `residue_map` to `build_network_signature()` when an
alignment supplies shared positions. Family profiles retain feature support
and repeated-contact multiplicity rather than collapsing equivalent ligands
to a single presence/absence flag.

This comparison layer is independent of Plotly and file output. That keeps a
future public interface straightforward: an upload or PDB-ID adapter can
produce analysis tables, the comparison layer can return JSON, and the
existing cohesive viewer can render the selected structure separately. Remote
structure retrieval and web request handling are deliberately kept outside
the core analysis functions.

For file-oriented workflows, `compare_coordination_networks.py` provides the
same operations without requiring Python code:

```bash
python compare_coordination_networks.py pairwise \
  --reference reference.pdb \
  --query query_a.pdb query_b.pdb \
  --cofactor CU \
  --output-dir comparison_output

python compare_coordination_networks.py profile \
  --reference known_family/ \
  --query candidates/ \
  --template known_family/template.pdb \
  --cofactor OEX \
  --output-dir profile_output
```

The runner writes canonical signatures, JSON results, ranked CSV summaries,
an `errors.csv` manifest, and both static PNG plus interactive HTML similarity
heatmaps for multi-structure runs. The earlier
`2_Single_Template_Single_Query_Network_Comparison...` and
`3_Single_Template_Multiple_Query_Network_Comparison...` scripts remain
available for legacy CSV/Kabsch alignment workflows; the new runner uses the
current motif-aware analysis tables as its source of truth. Alignment is on by
default in the runner: it tries cofactor anchors first, falls back to shared
network atoms, and records the alignment method, anchor RMSD, matched atom
count, and residue-position mappings. Use `--no-alignment` for a strictly
chemistry/topology-only comparison or `--alignment-cutoff` to adjust the
post-alignment mapping tolerance.

Use `--template` to reproduce the legacy template-numbering convention: the
template supplies the alignment coordinate system and the residue numbers
shown across the 2D conservation map. When omitted, the first reference
structure serves as the template.

Profile runs also write `profile_conservation_viewer.html`. This is a cohesive
single-structure viewer built from the first reference structure, with the
family profile overlaid as a selectable `Color by -> Family conservation`
mode. Focused atoms retain the normal network/layer/bond controls, while their
hover text reports the percentage of reference structures supporting the
corresponding residue. Use `--compact-html` when a smaller self-contained
viewer is preferred.

Aligned profile runs also write `residue_conservation_map.csv`,
`residue_conservation_map.png`, and `residue_conservation_map.html`. These are
two-dimensional frequency maps in the legacy orientation: residue types run
across the x-axis, template residue numbers run down the y-axis, and cell color
is the fraction of reference structures carrying that residue type at that
template position. Interactive hover includes the contributing structures and
full mapped position identity. If alignment is disabled or cannot establish
template positions, the map is recorded as unavailable in
`residue_conservation_map_details.json` rather than being inferred from raw
residue numbering.

## Stability

Modules are grouped into a stable core, frozen legacy compatibility, and
experimental features; see [`docs/stability.md`](docs/stability.md). Release
changes are listed in [`CHANGELOG.md`](CHANGELOG.md).

## Reference structures and validation

The repository includes small reference cases for regression and chemistry
checks:

- Plastocyanin `1ag6`: His37, Cys84, His87, and Met92 are PCS contacts.
- Myoglobin `1a6m`: proximal His93 is PCS and distal His64 is SCS.
- Carbonic anhydrase `1ca2`: the single Zn site includes His94, His96, and
  His119 as PCS contacts.
- Ferredoxin `2zvs`: six SF4 cofactors are parsed and cysteine ligands are
  detected in the PCS output.
- OEC `4ub6` and nitrogenase `3u7q_monomer` provide multi-cofactor regression
  coverage.

Run the regression and validation suite with:

```bash
python -m unittest discover -s tests -v
```

Measured chemistry checks and headless timings are summarized in
[`docs/validation.md`](docs/validation.md).

## Repository layout

```text
modules/                         Core analysis, comparison, wire, and labeling code
  coordination_api.py            Headless analysis tables and batch API
  motif_registry.py, cofactor_classes.py   Motif vocabulary and cofactor-class cutoffs
  network_comparison.py, network_alignment.py, comparison_runner.py,
  comparison_plots.py, conservation_viewer.py   Comparison and family profiles
  protein_wires.py, wire_viewer.py   Cofactor-to-target relay networks
  ml_export.py, assembly_policy.py, sensitivity.py   Labeling contract and splits
  substrate_seeds.py             Hypothetical substrate-point chains
  legacy_adapter.py              Legacy CSV -> current tables
benchmarks/                      Smoke, dataset-scale, and sensitivity benchmarks
tests/                           Regression and chemistry validation tests
reference_structures/            Public example structures
1_Single_..._SSCNA_v0.0.2.py    Legacy single-structure CLI
sscna_cli.py                     Stable `sscna` command (analyze, compare, wire)
batch_coordination_network.py   Headless batch CLI
compare_coordination_networks.py Pairwise and family-profile comparison CLI
```

Generated analysis outputs belong in ignored output directories. Do not add
credentials, private structures, or local development notes to the public
repository.

## License

[BSD Zero Clause (0BSD)](LICENSE): use, copy, modify, and distribute this
software for any purpose, with or without fee, with no attribution
requirement. Earlier releases were published under the GPLv3.
