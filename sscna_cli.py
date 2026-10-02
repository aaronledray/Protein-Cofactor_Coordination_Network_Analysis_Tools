#!/usr/bin/env python3
"""Stable installed command-line interface for SSCNA.

``sscna analyze`` intentionally delegates to the legacy single-structure
entry point.  This keeps its defaults, sidecar configuration behavior, and
output schema stable while providing a command name that does not depend on
the versioned filename.  ``sscna compare`` exposes the newer pairwise and
family-profile workflows; ``sscna wire`` exposes opt-in cofactor-to-target
protein-wire analysis; ``sscna substrate`` traces opt-in chains from a
hypothetical substrate point through the coordination network.
"""

from __future__ import annotations

import argparse
import json
from importlib.metadata import distribution, PackageNotFoundError
import runpy
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence


LEGACY_FILENAME = "1_Single_Structure_Cofactor_Network_Analysis_SSCNA_v0.0.2.py"


def _legacy_script_path() -> Path:
    """Find the legacy source both in a checkout and in an installed wheel."""
    checkout_path = Path(__file__).with_name(LEGACY_FILENAME)
    if checkout_path.is_file():
        return checkout_path

    try:
        installed_distribution = distribution("coordination-network-identifier")
    except PackageNotFoundError as exc:
        raise RuntimeError("The coordination-network-identifier distribution is not installed") from exc
    installed_files = installed_distribution.files or ()
    for installed_file in installed_files:
        if installed_file.name == LEGACY_FILENAME:
            candidate = Path(installed_distribution.locate_file(installed_file))
            if candidate.is_file():
                return candidate
    raise RuntimeError(f"Installed legacy analyzer is missing: {LEGACY_FILENAME}")


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


def _run_legacy_analyzer(args: Sequence[str]) -> int:
    legacy_script = _legacy_script_path()

    previous_argv = sys.argv
    try:
        sys.argv = [str(legacy_script), *args]
        try:
            runpy.run_path(str(legacy_script), run_name="__main__")
        except SystemExit as exc:
            if exc.code in (None, 0):
                return 0
            if isinstance(exc.code, int):
                return exc.code
            print(exc.code, file=sys.stderr)
            return 1
    finally:
        sys.argv = previous_argv
    return 0


def _analyze(argv: Sequence[str]) -> int:
    if any(argument in {"-h", "--help"} for argument in argv):
        return _run_legacy_analyzer(["--help"])

    parser = argparse.ArgumentParser(
        prog="sscna analyze",
        description=(
            "Run the legacy single-structure analyzer with a stable command name. "
            "All legacy options are forwarded unchanged."
        ),
        add_help=False,
    )
    parser.add_argument("--input", help="Input PDB/mmCIF structure.")
    parser.add_argument(
        "--template",
        help="Legacy alias for --input; retained for compatibility.",
    )
    parsed, forwarded = parser.parse_known_args(list(argv))

    if parsed.input and parsed.template:
        parser.error("use only one of --input or --template")
    structure = parsed.input or parsed.template
    if not structure:
        parser.error("the following argument is required: --input")

    return _run_legacy_analyzer(["--template", structure, *forwarded])


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
