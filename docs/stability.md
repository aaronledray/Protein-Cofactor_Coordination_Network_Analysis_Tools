# Stability tiers

Modules fall into three tiers. A pipeline that needs reproducible labels
(for example model training) should depend only on the stable core and pin a
release tag.

## Stable core

Covered by the labeling contract and golden fixtures; breaking changes bump
`LABELING_CONTRACT_VERSION` or the package's minor version and are listed in
`CHANGELOG.md`. These modules import no plotting code and work without plotly
or matplotlib installed (enforced by `tests/test_module_boundaries.py`).

| Module | Role |
|---|---|
| `coordination_api` | `analyze_structure()`, `batch_analyze()`, the four tables |
| `ml_export` | `analyze_labeling_example()`, `export_residue_labels()`, contract 1.0 |
| `assembly_policy` | symmetry context and equivalent-site grouping |
| `sensitivity` | label stability under alternative parameters |
| `cofactor_classes`, `motif_registry` | cofactor families, cutoffs, motif vocabulary |
| `structure_processing`, `chemistry`, `moieties`, `structure_utils`, `io_utils` | shell identification and chemistry tables |

## Legacy compatibility

Frozen unless a change is explicitly approved: the versioned single-structure
script, `sscna analyze`, `batch_coordination_network.py`, and their default
output filenames and CSV schemas (pinned by `tests/fixtures/baseline/` and
`tests/test_stable_cli.py`). New behavior here is opt-in only
(`--cohesive-viewer`, `--shells`, `--per-site`, ...). The one default change in
0.2.0 is that figure windows are not opened for non-interactive runs.

## Experimental

Useful and tested, but their interfaces and outputs may change between minor
releases and they are not part of the labeling contract:

- comparison and profiles: `network_comparison`, `network_alignment`,
  `comparison_runner`, `comparison_plots`, `conservation_viewer`,
  `legacy_adapter`, `compare_coordination_networks.py`
- protein wires: `protein_wires`, `wire_viewer`, `sscna wire`
- substrate seeds: `substrate_seeds`, `substrate_viewer`, `sscna substrate`
- plotting and display: `plotting` (including the cohesive viewer), `analysis`,
  `display`

Wire and substrate results are structural hypotheses, not rate, binding, or
reactivity calculations.

## Adding a module

Place it in a tier in `tests/test_module_boundaries.py` and in this file. A
stable-core module may import only other stable-core modules.
