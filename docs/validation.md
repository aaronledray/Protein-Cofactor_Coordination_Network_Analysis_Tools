# Validation notes

These checks are intentionally small and are also encoded in
`tests/test_chemistry_validation.py`.

## Known sites

- 1ag6 plastocyanin: His37, Cys84, His87, and Met92 are PCS contacts.
- 1a6m myoglobin: proximal His93 is PCS and distal His64 is SCS.
- 1ca2 carbonic anhydrase: the single Zn site identifies His94, His96, and
  His119 as PCS contacts.
- 2zvs ferredoxin: six SF4 residues are parsed, and at least twelve Cys
  residues are identified in the PCS table across the three chains.

The Zn and Fe–S reference entries are public RCSB PDB structures:
[1CA2](https://www.rcsb.org/structure/1CA2) and
[2ZVS](https://www.rcsb.org/structure/2ZVS).

## Headless timing

Measured on the development machine with `modules.coordination_api` and
default two-shell settings, excluding plot generation and file writes:

| Structure | PCS | SCS | Links | Runtime |
|---|---:|---:|---:|---:|
| 1ag6 | 5 | 2 | 7 | 0.032 s |
| 1a6m | 17 | 15 | 32 | 0.029 s |
| 4ub6 | 21 | 16 | 37 | 0.038 s |
| 3u7q | 38 | 62 | 100 | 0.164 s |

The combined elapsed time for the four sequential API calls was 0.265 s.
These are indicative local timings, not performance guarantees across
hardware, Python versions, or file systems.

## Interpretation notes

The 2zvs case contains six SF4 cofactors, so union mode reports a combined
network. Use `site_mode="per-site"` or `--per-site` when labels must remain
attached to individual cofactor sites. The tool completed both the Zn and
Fe–S cases without parser or classification exceptions.

Site-boundary regression tests also verify that insertion codes remain
distinct, pooled model identity remains backward-compatible, per-model site
IDs can be requested explicitly, and combinatorial cofactors are never merged
solely because they are close in different structure models.

## Template-anchored profile validation

The comparison runner was exercised with explicit template structures and
self-query profiles across the available representative families. All cases
completed with zero recorded input or analysis errors and used cofactor-first
alignment anchors.

| Family / case | Template | Reference structures | Query structures | Anchor RMSD | Conservation positions |
|---|---|---:|---:|---:|---:|
| OEX family | 4ub6 | 13 | 13 | 0.000–0.400 Å | 49 |
| Heme | 1a6m | 1 | 1 | 0.000 Å | 46 |
| Zn | 1ca2 | 1 | 1 | 0.000 Å | 13 |
| SF4 / ferredoxin | 2zvs | 1 | 1 | 0.000 Å | 109 |
| Nitrogenase, combinatorial ICS + CLF + HCA | 3u7q monomer | 1 | 1 | 0.000 Å | 136 |

The OEX family self-query profile scores ranged from 0.583 to 0.807 across
the 13 structures; single-reference cases scored 1.000 against themselves.
The generated artifacts are written under the ignored
`network_profile_output/validation/` directory during local validation.

## Batch smoke benchmark

`python benchmarks/batch_smoke_benchmark.py` runs `analyze_structure` (3 shells)
over Cu, heme, Zn, per-site Fe-S, combinatorial nitrogenase, and OEC cases,
plus a deliberately missing input. It writes CSV/JSON to the ignored
`benchmark_output/` directory and checks that `batch_analyze` with 2 workers
equals the serial result and records the missing file in the error table.

Measured on one Apple-silicon laptop (single run, so treat as indicative):

| Case | Sites | Residue rows | Contacts | Wall (s) | Peak Python MB |
|---|---:|---:|---:|---:|---:|
| Cu plastocyanin 1ag6 | 1 | 9 | 8 | 0.08 | 2.4 |
| Heme myoglobin 1a6m | 1 | 46 | 76 | 0.11 | 5.8 |
| Zn carbonic anhydrase 1ca2 | 1 | 13 | 12 | 0.08 | 5.4 |
| SF4 ferredoxin 2zvs (per-site) | 6 | 110 | 157 | 0.27 | 5.9 |
| Nitrogenase 3u7q monomer (combinatorial) | 1 | 136 | 224 | 0.67 | 41.9 |
| OEC 4ub6 (per-site, first model) | 1 | 49 | 109 | 0.16 | 7.9 |

These small cases say nothing about dataset-scale throughput; the larger
benchmark on the intended PDB subset is still outstanding.

## Label sensitivity

`modules/sensitivity.py` compares exported residue labels under alternative
parameters against a baseline (3.6 Å, no carbon seeds, 3 shells, per-site), and
`python benchmarks/sensitivity_report.py` writes `benchmark_output/label_sensitivity.csv`
for the benchmark cases. Primary-shell Jaccard overlap versus baseline:

| Case | 3.2 Å | 3.4 Å | 3.8 Å | 4.0 Å | Carbon seeds |
|---|---:|---:|---:|---:|---:|
| Cu 1ag6 | 0.80 | 0.80 | 1.00 | 1.00 | 1.00 |
| Heme 1a6m | 0.71 | 0.77 | 0.90 | 0.74 | 0.85 |
| Zn 1ca2 | 1.00 | 1.00 | 1.00 | 0.67 | 1.00 |
| SF4 2zvs | 0.84 | 0.91 | 0.87 | 0.68 | 0.70 |
| Nitrogenase 3u7q | 0.62 | 0.79 | 0.94 | 0.85 | 0.81 |
| OEC 4ub6 | 1.00 | 1.00 | 0.94 | 0.88 | 1.00 |

Findings:
- Primary labels are most stable for single-metal sites at nearby cutoffs; the
  heme, Fe-S, and nitrogenase cases move noticeably (Jaccard 0.62–0.91) for
  ±0.2 Å, so cutoff and carbon-seed policy must be recorded with exported labels.
- Tightening the cutoff only removes residues and loosening it only adds them in
  the Cu case; shell reassignment of kept residues is small but nonzero.
- `expand_residues` and `exclude_moieties=["alanine_sidechain"]` did not change
  residue-level labels in any benchmark case (atom-level tables may still differ).
- Water is not a variant: there is no water option. Crystallographic waters are ordinary shell members (motif `water`; 1a6m has 18 in PCS-TCS), and only the moiety-bond expansion skips solvent. Exported labels therefore include HOH rows.

## Dataset-scale benchmark

`python benchmarks/dataset_benchmark.py --input DIR --cofactor HEM --limit N --workers W`
runs `analyze_labeling_example` (contract 1.0: 3 shells, per-site, per-model)
over a directory, with per-structure failure isolation, and writes
`benchmark_output/dataset_benchmark_*.{csv,json}`.

Measured on the 2,755 heme-protein PDB entries in the local
`Labtools_apl/.../4_Xlinks_Heme/heme` collection (8 workers, 10-core Apple
silicon laptop). This is the designated benchmark set; if the training
subset later differs, repeat the run on it.

| Pass | Structures | OK | No HEM | Failed | Wall | Throughput | Median / p95 / max per structure | Max worker RSS |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| Serial subset | 50 | 50 | 0 | 0 | 12 s | 4.2/s | 0.09 / 1.0 / 3.8 s | 585 MB |
| 8 workers, sample | 200 | 200 | 0 | 0 | 11.5 s | 17.4/s | 0.14 / 1.3 / 8.5 s | 500 MB |
| 8 workers, full | 2,755 | 2,753 | 2 | 0 | 103 s | 26.7/s | 0.13 / 1.3 / 10.9 s | 600 MB |
| 8 workers, full (rerun, final code) | 2,755 | 2,753 | 2 | 0 | 125 s | 22.0/s | 0.14 / 1.4 / 10.6 s | 480 MB |
| Serial subset (rerun) | 50 | 50 | 0 | 0 | 12 s | 4.1/s | 0.11 / 1.1 / 2.8 s | 339 MB |

The full run was repeated on the final code and produced the same 214,637 residue-label rows, which confirms the labels are deterministic; wall time varied from 103 s to 125 s between runs. The full run produced 214,637 residue-label rows across up to 20 sites per
structure with no exceptions; the two structures without HEM were reported as
`no_cofactor`, not failures. Cost is heavy-tailed: the slowest entries are
large or many-model structures (for example 6NIN, 2GB8, 2N18), and memory
scales with them. Parallel speedup over serial on the sample is about 4x on 10
cores, so a dataset 10x larger would take roughly 20 minutes at this rate.
