#!/usr/bin/env python3
"""Headless batch entry point for coordination-network analysis."""

import argparse
from pathlib import Path

from modules.coordination_api import batch_analyze


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze PDB/mmCIF structures and write combined coordination tables."
    )
    parser.add_argument("--input", nargs="+", required=True, help="Structure files or directories")
    parser.add_argument("--cofactor", required=True, help="Cofactor residue names, comma-separated")
    parser.add_argument("--cofactor2", help="Second cofactor residue names, comma-separated")
    parser.add_argument("--distance", type=float, default=3.6)
    parser.add_argument("--shells", type=int, default=2)
    parser.add_argument("--per-site", action="store_true", help="Analyze each cofactor site separately")
    parser.add_argument("--expand-residues", action="store_true")
    parser.add_argument("--first-model", action="store_true", help="Analyze only the first structure model")
    parser.add_argument("--combinatorial", action="store_true")
    parser.add_argument("--combinatorial-cutoff", type=float, default=20.0)
    parser.add_argument("--exclude-moieties", default="alanine_sidechain")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output-dir", default="coordination_batch_output")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result = batch_analyze(
        args.input,
        args.cofactor,
        workers=args.workers,
        distance_cutoff=args.distance,
        shells=args.shells,
        site_mode="per-site" if args.per_site else "union",
        expand_residues=args.expand_residues,
        first_model_only=args.first_model,
        combinatorial=args.combinatorial,
        combinatorial_cofactor_cutoff=args.combinatorial_cutoff,
        cofactor_resname2=args.cofactor2,
        exclude_moieties=[part.strip() for part in args.exclude_moieties.split(",") if part.strip()],
    )
    result["residues"].to_csv(output_dir / "coordination_residues.csv", index=False)
    result["atoms"].to_csv(output_dir / "coordination_atoms.csv", index=False)
    result["links"].to_csv(output_dir / "coordination_links.csv", index=False)
    result["errors"].to_csv(output_dir / "coordination_errors.csv", index=False)
    print(
        f"Processed {result['residues']['structure_id'].nunique()} structure(s); "
        f"errors={len(result['errors'])}; output={output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
