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
