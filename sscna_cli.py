#!/usr/bin/env python3
"""Stable installed command-line interface for SSCNA.

``sscna analyze`` builds the tidy coordination tables and the cohesive
interactive viewer for one structure (``modules/analyze_cli.py``).
``sscna compare`` exposes the newer pairwise and
family-profile workflows; ``sscna wire`` exposes opt-in cofactor-to-target
protein-wire analysis; ``sscna substrate`` traces opt-in chains from a
hypothetical substrate point through the coordination network.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence


def _help_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sscna",
        description="Analyze and compare protein cofactor coordination networks.",
    )
    parser.add_argument(
        "command",
        choices=("analyze", "compare", "wire", "substrate"),
        help="analyze one structure, compare structures/profiles, find cofactor-to-target wires, or trace substrate-point chains",
    )
    return parser


def _analyze(argv: Sequence[str]) -> int:
    from modules.analyze_cli import main as analyze_main

    return analyze_main(list(argv))


def _parse_wire_target(value: str) -> Dict[str, Any]:
    """Parse ``RESNAME:CHAIN:RESNUM[:ATOM]`` into a target selector."""
    parts = [part.strip() for part in str(value).split(":")]
    if len(parts) not in {3, 4} or not all(parts):
        raise argparse.ArgumentTypeError(
            "target must use RESNAME:CHAIN:RESNUM[:ATOM], e.g. TRP:A:107:NE1"
        )
    selector: Dict[str, Any] = {
        "label": value,
        "residue": parts[0],
        "chain": parts[1],
        "residue_number": parts[2],
    }
    if len(parts) == 4:
        selector["atom"] = parts[3]
    return selector


def _wire(argv: Sequence[str]) -> int:
    from modules.protein_wires import analyze_protein_wires

    parser = argparse.ArgumentParser(
        prog="sscna wire",
        description="Find opt-in cofactor-to-target relay or redox-relay paths.",
    )
    parser.add_argument("--input", required=True, help="Input PDB/mmCIF structure.")
    parser.add_argument("--cofactor", required=True, help="Cofactor residue name(s), comma-separated.")
    parser.add_argument(
        "--target",
        action="append",
        required=True,
        type=_parse_wire_target,
        metavar="RES:CHAIN:NUM[:ATOM]",
        help="Target selector; repeat for multiple targets, e.g. TRP:A:107:NE1.",
    )
    parser.add_argument("--distance", type=float, default=3.6, help="Maximum wire-hop distance in Å (default: 3.6).")
    parser.add_argument("--max-hops", type=int, default=8, help="Maximum hops per path (default: 8).")
    parser.add_argument(
        "--mode",
        choices=("relay", "generic", "proton", "electron", "pcet", "redox"),
        default="relay",
        help="Shared relay graph or interpretation: relay, proton, electron, pcet, or residue-level redox (default: relay).",
    )
    parser.add_argument(
        "--redox-distance",
        type=float,
        default=6.0,
        help="Maximum aromatic/redox hop distance in Å for --mode redox (default: 6.0).",
    )
    parser.add_argument("--no-water", action="store_true", help="Exclude water molecules as relay nodes.")
    parser.add_argument("--include-backbone", action="store_true", help="Allow backbone atoms as relay nodes.")
    parser.add_argument("--first-model", action="store_true", help="Analyze only the first model.")
    parser.add_argument("--no-html", action="store_true", help="Skip the interactive HTML viewer.")
    parser.add_argument(
        "--output-dir",
        default="protein_wire_output",
        help="Directory for wire_nodes.csv, wire_edges.csv, and wire_paths.csv.",
    )
    args = parser.parse_args(list(argv))

    result = analyze_protein_wires(
        args.input,
        args.cofactor,
        args.target,
        max_hop_distance=args.distance,
        max_hops=args.max_hops,
        include_water=not args.no_water,
        include_backbone=args.include_backbone,
        first_model_only=args.first_model,
        wire_mode=args.mode,
        redox_hop_distance=args.redox_distance,
    )
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    result["nodes"].to_csv(output_dir / "wire_nodes.csv", index=False)
    result["edges"].to_csv(output_dir / "wire_edges.csv", index=False)
    result["paths"].to_csv(output_dir / "wire_paths.csv", index=False)
    result["target_status"].to_csv(output_dir / "wire_targets.csv", index=False)
    if not args.no_html:
        from modules.wire_viewer import plot_interactive_protein_wire_network

        plot_interactive_protein_wire_network(
            result["nodes"],
            result["edges"],
            result["paths"],
            output_filename=str(output_dir / "wire_network.html"),
            pdb_name=Path(args.input).name,
            cofactor_resname=args.cofactor,
            wire_mode=args.mode,
        )
    summary = {
        "input": str(Path(args.input).expanduser().resolve()),
        "cofactor": args.cofactor,
        "targets": args.target,
        "max_hop_distance_A": args.distance,
        "redox_hop_distance_A": args.redox_distance if args.mode == "redox" else None,
        "max_hops": args.max_hops,
        "wire_mode": args.mode,
        "include_water": not args.no_water,
        "include_backbone": args.include_backbone,
        "html_written": not args.no_html,
        "node_count": len(result["nodes"]),
        "edge_count": len(result["edges"]),
        "path_count": len(result["paths"]),
        "target_status": result["target_status"].to_dict(orient="records"),
        "unresolved_targets": result["target_status"].loc[
            ~result["target_status"]["matched"].astype(bool), "target_label"
        ].tolist(),
    }
    (output_dir / "wire_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


def _substrate(argv: Sequence[str]) -> int:
    from modules.substrate_seeds import DEFAULT_CONTACT_CUTOFF_A, analyze_substrate_seed

    parser = argparse.ArgumentParser(
        prog="sscna substrate",
        description=(
            "Trace chains from a hypothetical substrate point through the "
            "coordination network to the cofactor (structural hypotheses only)."
        ),
    )
    parser.add_argument("--input", required=True, help="Input PDB/mmCIF structure.")
    parser.add_argument("--cofactor", required=True, help="Cofactor residue name(s), comma-separated.")
    where = parser.add_mutually_exclusive_group(required=True)
    where.add_argument("--point", nargs=3, type=float, metavar=("X", "Y", "Z"),
                       help="Explicit substrate-point coordinates in Å.")
    where.add_argument("--from-atom", metavar="ATOM",
                       help="Cofactor atom to offset from (resolved within each site); needs --offset.")
    parser.add_argument("--offset", nargs=3, type=float, metavar=("DX", "DY", "DZ"),
                        help="Offset in Å from --from-atom.")
    parser.add_argument("--contact-cutoff", type=float, default=DEFAULT_CONTACT_CUTOFF_A,
                        help=f"Point-to-network-atom contact distance in Å (default: {DEFAULT_CONTACT_CUTOFF_A}).")
    parser.add_argument("--max-chains", type=int, default=10, help="Maximum ranked chains per site (default: 10).")
    parser.add_argument("--distance", type=float, default=3.6, help="Shell distance cutoff in Å (default: 3.6).")
    parser.add_argument("--shells", type=int, default=3, help="Number of shells to build (default: 3).")
    parser.add_argument("--first-model", action="store_true", help="Analyze only the first model.")
    parser.add_argument("--no-html", action="store_true", help="Skip the interactive HTML viewer.")
    parser.add_argument("--compact-html", action="store_true", help="Use a CDN-backed Plotly bundle.")
    parser.add_argument("--output-dir", default="substrate_seed_output",
                        help="Directory for substrate_*.csv, substrate_summary.json, and the viewer.")
    args = parser.parse_args(list(argv))
    if args.from_atom and not args.offset:
        parser.error("--from-atom requires --offset DX DY DZ")
    if args.offset and not args.from_atom:
        parser.error("--offset requires --from-atom")

    point: Any = (
        {"cofactor_atom": args.from_atom, "offset": tuple(args.offset)}
        if args.from_atom
        else tuple(args.point)
    )
    result = analyze_substrate_seed(
        args.input,
        [name.strip() for name in args.cofactor.split(",") if name.strip()],
        point,
        contact_cutoff=args.contact_cutoff,
        max_chains=args.max_chains,
        shells=args.shells,
        distance_cutoff=args.distance,
        first_model_only=args.first_model,
    )
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    result["seed"].to_csv(output_dir / "substrate_seed.csv", index=False)
    result["substrate_contacts"].to_csv(output_dir / "substrate_contacts.csv", index=False)
    result["chains"].to_csv(output_dir / "substrate_chains.csv", index=False)
    result["steps"].to_csv(output_dir / "substrate_steps.csv", index=False)
    html_written = False
    if not args.no_html and not result["seed"].empty:
        from modules.substrate_viewer import plot_interactive_substrate_chains

        plot_interactive_substrate_chains(
            result,
            str(output_dir / "substrate_chains.html"),
            pdb_name=Path(args.input).name,
            compact_html=args.compact_html,
        )
        html_written = True
    summary = {
        "input": str(Path(args.input).expanduser().resolve()),
        "cofactor": args.cofactor,
        "point": point,
        "contact_cutoff_A": args.contact_cutoff,
        "shells": args.shells,
        "site_count": len(result["seed"]),
        "contact_count": len(result["substrate_contacts"]),
        "chain_count": len(result["chains"]),
        "html_written": html_written,
        "note": "Chains are structural hypotheses about contact connectivity, not binding or reactivity predictions.",
    }
    (output_dir / "substrate_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in {"-h", "--help"}:
        _help_parser().print_help()
        return 0

    command, command_args = args[0], args[1:]
    if command == "analyze":
        return _analyze(command_args)
    if command == "compare":
        from compare_coordination_networks import main as compare_main

        return compare_main(command_args)
    if command == "wire":
        return _wire(command_args)
    if command == "substrate":
        return _substrate(command_args)
    _help_parser().error(f"invalid command: {command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
