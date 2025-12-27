# Coordination Network Identifier

Tools for extracting coordination networks around protein cofactors and comparing them across structures. The codebase is centered on Biopython utilities with plotting helpers for quick inspection.



---



## What even is this?
- Single-structure coordination analysis (primary and secondary spheres) with optional residue expansion, combinatorial multi-cofactor filtering, and moiety exclusions.
- CSV, HTML, and PNG outputs that capture coordination atoms, moiety labels, and simple network plots.
- Template-to-query alignment that takes a template coordination CSV plus a candidate structure, aligns common atoms, and records RMSD.
- Chain deduplication to clean PDB/mmCIF inputs by sequence or by atomic coordinates, with ligand-based masking.




## Motivation:

- I got tired of manually finding interactions between proteins and their cofactors.
- Needed a robust definition of coordination networks; literature's "seconady coordination sphere" was too vague for protein design workflows (rather for my applications, thereof)!



---







## Requirements:
- Python 3.9 or newer.
- Packages: biopython, numpy, pandas, matplotlib, plotly.








## Repository layout
- `modules/`: core logic for coordination finding, IO, plotting, deduplication.
- `output/`: sample outputs from prior runs.
- `matches/`: alignment job folders created by the comparison script.
- `data/`: example structures.
- Legacy single-file scripts are kept in the project root and `Archived/` for reference.








## Set up a virtual environment and install:

```bash
python -m venv .venv
source .venv/bin/activate
pip install biopython numpy pandas matplotlib plotly
```



Activate:

``` zsh
source .venv/bin/activate
```

then when ready to revert to SYSTEM python, use:

``` zsh
deactivate
```









---

# Running:

## Single-structure analysis (SSCNA)
Script: `1_Single_Structure_Cofactor_Network_Analysis_SSCNA_v0.0.2.py`

Example (non-interactive):
```bash
python 1_Single_Structure_Cofactor_Network_Analysis_SSCNA_v0.0.2.py \
  --template ./data/1_OEX/0_Aligned_Reduced/4ub6.pdb \
  --cofactor OEX \
  --distance 3.6 \
  --combinatorial \
  --combinatorial-cutoff 20.0 \
  --exclude-moieties alanine_sidechain \
  --mode Coord_Network
```


Key flags:
- `--template` path to the structure (PDB/mmCIF, gz accepted).
- `--cofactor` (and `--cofactor2` when using combinatorial mode) comma-separated residue names.
- `--distance` cutoff for moiety interactions (Å).
- `--expand-residues` to include full residues for PCS/SCS atoms.
- `--combinatorial` plus `--combinatorial-cutoff` to keep cofactors near one another.
- `--exclude-moieties` comma-separated labels to ignore (e.g., alanine_sidechain).
- `--mode` `Coord_Network` (default) or `Residues_of_Interest`.
- `--interactive` prompts for missing parameters.

Outputs (written to the working directory):
- `Coord_Breakdown.csv` and `Coord_Links.csv` style tables describing cofactor, PCS, and SCS atoms.
- Static plots: `1_static_*` PNG files.
- Interactive Plotly: `1_template_coordination_network.html` (and ROI variant when using `Residues_of_Interest` mode).

## Template vs. query alignment
Script: `2_Single_Template_Single_Query_Network_Comparison_v0.0.2.py`

This consumes a template coordination CSV and a query structure, aligns common atoms with Kabsch, and logs RMSD.

Example:
```bash
python 2_Single_Template_Single_Query_Network_Comparison_v0.0.2.py \
  --template-csv output/4ub6.pdb_Coord_Breakdown.csv \
  --query-structure data/3_OEX_Designed_Protein_Candidates/.../model.cif \
  --job-name OEX_testmatch1
```

Behavior:
- Writes `matches/<job-name>/accessory_names.txt` with allowed residue/atom synonyms (edit if needed).
- Appends alignment results to `matches/<job-name>/alignment_report.csv` with atoms used and RMSD.
- Requires at least three shared atom names between template and query.

## Chain deduplication helper
Script: `modules/deduplicate_chains.py` (run directly).

Use it to collapse duplicate chains before coordination analysis.
```bash
python modules/deduplicate_chains.py 5jqr.cif
python modules/deduplicate_chains.py ./structures --by coords --mask-res HEM,CU1 --yes
```

Options:
- `--by` `sequence` (default) or `coords` for exact coordinate signatures.
- `--mask-res` comma-separated residue names; when duplicates exist, keeps the chain containing those residues (`--mask-mode` all|any).
- `--yes` for non-interactive runs; `--dry-run` to report without writing.
- Writes `<input>_dedup.pdb` (or similar) next to the source or to `--outdir`.
