# Coordination Network Identifier

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
  annotations, carbon-seed inclusion, and class-specific distance cutoffs.
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
python -m pip install -r requirements.txt
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

Outputs are written under `SSCNA_output/`:

- `<structure>_Coord_Breakdown.csv`: residue summary plus atom coordinates.
- `<structure>_Coord_Breakdown_atoms.csv`: atom-level category report.
- `<structure>_Coord_Links.csv`: nearest cofactor→PCS and PCS→SCS links.
- PNG and HTML plots, unless `--no-plots` is supplied.

Useful additive options include:

```text
--no-plots                         Skip PNG/HTML rendering.
--first-model                      Analyze only the first model.
--shells N                         Add TCS and deeper shells for N > 2.
--per-site                         Separate cofactor sites in tidy output.
--include-carbon-seeds             Allow carbon atoms to seed shells.
--direct-coordination              Annotate direct metal-ligand links.
--direct-coordination-cutoff 2.6   Set the direct-link distance in Å.
--cofactor-class-cutoff metal=2.8  Override a class cutoff; repeatable.
--verbose                          Show progress logs.
```

The legacy command currently scans all structure models unless
`--first-model` is supplied. Alternate locations are not assigned a custom
policy; Biopython's normal disordered-atom selection behavior is used.

For `--shells N` with `N > 2`, the command writes tidy
`Coordination_Residues.csv`, `Coordination_Atoms.csv`, and `Coord_Links.csv`
tables and does not attempt the legacy two-shell plots.

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

`analyze_structure()` creates no files and returns `residues`, `atoms`, and
`links` pandas DataFrames. The residue table has one row per residue per site
and shell, including structure ID, cofactor identity, residue identity,
atoms involved, and minimum distance to the previous shell. Insertion codes
and hetero flags are retained in the tidy tables.

Use `site_mode="per-site"` to separate cofactor copies. With
`combinatorial=True`, cofactor residues within the combinatorial cutoff are
clustered into one site; this is useful for multi-cofactor systems such as
nitrogenase.

## Batch analysis

```bash
python batch_coordination_network.py \
  --input reference_structures/0_Plastocyanin \
  --cofactor CU \
  --workers 4 \
  --output-dir coordination_batch_output
```

The batch command writes combined `coordination_residues.csv`,
`coordination_atoms.csv`, `coordination_links.csv`, and
`coordination_errors.csv`. A malformed or missing structure is recorded in
the error table without stopping other structures.

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
modules/                         Core analysis and reporting code
tests/                           Regression and chemistry validation tests
reference_structures/            Public example structures
1_Single_..._SSCNA_v0.0.2.py    Legacy single-structure CLI
batch_coordination_network.py   Headless batch CLI
```

Generated analysis outputs belong in ignored output directories. Do not add
credentials, private structures, or local development notes to the public
repository.

## License

See [LICENSE](LICENSE).
