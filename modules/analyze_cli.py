"""Single-structure coordination analysis command (``sscna analyze``).

Builds the tidy coordination tables with :func:`modules.coordination_api.analyze_structure` and, unless plots are
skipped, renders the cohesive interactive viewer. This replaces the versioned legacy script and its
``Coord_Breakdown`` report schema.

Outputs, written to ``--output-dir`` (default ``SSCNA_output`` in the working directory) with the input file name as
prefix: ``Coordination_Residues.csv``, ``Coordination_Atoms.csv``, ``Coord_Links.csv``, ``Coord_Contacts.csv``, and
``coordination_network.html`` (not written with ``--no-plots``).
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from typing import Any, Dict, List, Optional, Sequence

from .cofactor_classes import load_cofactor_class_config, merge_cofactor_class_configs
from .io_utils import unpack_pdb_file

logger = logging.getLogger(__name__)

DEFAULT_SHELLS = 3
DEFAULT_DISTANCE = 3.6
DEFAULT_COMBINATORIAL_CUTOFF = 20.0
DEFAULT_DIRECT_COORDINATION_CUTOFF = 2.6
DEFAULT_EXCLUDED_MOIETIES = ["alanine_sidechain"]


def resolve_structure_path(user_path: Optional[str]) -> Optional[str]:
    """Return a path to an existing structure file, trying common extensions and ``.gz`` variants."""
    if not user_path:
        return None
    path = os.path.abspath(os.path.expanduser(user_path))
    if os.path.isfile(path):
        return path

    base, ext = os.path.splitext(path)
    if ext == ".gz":
        base = os.path.splitext(base)[0]

    extensions = (".pdb", ".cif", ".mmcif")
    candidates = []
    for extension in extensions:
        candidates.extend([base + extension, base + extension + ".gz"])
    # A bare PDB identifier such as "1ag6" is looked up in the working directory.
    if os.path.sep not in user_path and len(os.path.basename(base)) in (4, 5):
        for extension in extensions:
            name = os.path.basename(base) + extension
            candidates.extend([os.path.join(os.getcwd(), name), os.path.join(os.getcwd(), name + ".gz")])
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    return None


def split_names(value: Optional[str]) -> List[str]:
    """Split a comma-separated command-line value into stripped, non-empty names."""
    return [item.strip() for item in (value or "").split(",") if item.strip()]


def parse_class_cutoffs(values: Optional[Sequence[str]]) -> Dict[str, float]:
    """Parse repeated ``CLASS=ANGSTROMS`` options."""
    cutoffs: Dict[str, float] = {}
    for value in values or []:
        if "=" not in value:
            raise ValueError("--cofactor-class-cutoff must use CLASS=ANGSTROMS")
        name, cutoff = value.split("=", 1)
        cutoffs[name.strip().lower()] = float(cutoff)
    return cutoffs


def atom_groups_by_shell(tables: Dict[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    """Group the tidy atom rows by shell label for the viewer."""
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for row in tables["atoms"].to_dict("records"):
        groups.setdefault(str(row.get("shell", "Unknown")), []).append(row)
    return groups


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sscna analyze",
        description="Identify the coordination network around a cofactor in one structure.",
    )
    parser.add_argument("--input", "--template", dest="input", required=True, help="Input PDB/mmCIF structure (.gz accepted).")
    parser.add_argument("--cofactor", required=True, help="Cofactor residue name(s), comma-separated.")
    parser.add_argument("--cofactor2", help="Second cofactor residue name(s), comma-separated.")
    parser.add_argument("--distance", type=float, default=DEFAULT_DISTANCE, help=f"Moiety interaction cutoff in Å (default {DEFAULT_DISTANCE}).")
    parser.add_argument("--shells", type=int, default=DEFAULT_SHELLS, help=f"Number of coordination shells (default {DEFAULT_SHELLS}).")
    parser.add_argument("--per-site", action="store_true", help="Analyze each cofactor site separately.")
    parser.add_argument("--site-model-mode", choices=["pooled", "per-model"], default="pooled",
                        help="Boundary policy for per-site analyses across structure models.")
    expand = parser.add_mutually_exclusive_group()
    expand.add_argument("--expand-residues", dest="expand_residues", action="store_true", default=False,
                        help="Expand shell members to entire residues.")
    expand.add_argument("--no-expand-residues", dest="expand_residues", action="store_false")
    combinatorial = parser.add_mutually_exclusive_group()
    combinatorial.add_argument("--combinatorial", dest="combinatorial", action="store_true", default=False,
                               help="Enable combinatorial cofactor filtering.")
    combinatorial.add_argument("--no-combinatorial", dest="combinatorial", action="store_false")
    parser.add_argument("--combinatorial-cutoff", dest="combinatorial_cofactor_cutoff", type=float,
                        default=DEFAULT_COMBINATORIAL_CUTOFF, help="Cutoff for combinatorial cofactor filtering in Å.")
    parser.add_argument("--exclude-moieties", default=",".join(DEFAULT_EXCLUDED_MOIETIES),
                        help="Comma-separated moieties to exclude (default alanine_sidechain; empty string for none).")
    parser.add_argument("--first-model", action="store_true", help="Analyze only the first structure model.")
    parser.add_argument("--include-carbon-seeds", action="store_true", help="Allow carbon atoms to seed shells.")
    parser.add_argument("--direct-coordination", action="store_true", help="Annotate direct metal-ligand links.")
    parser.add_argument("--direct-coordination-cutoff", type=float, default=DEFAULT_DIRECT_COORDINATION_CUTOFF,
                        help="Direct metal-ligand distance cutoff in Å.")
    parser.add_argument("--cofactor-class-cutoff", action="append", help="Class-specific cutoff, e.g. metal=2.8 (repeatable).")
    parser.add_argument("--cofactor-class-config", help="YAML/JSON file defining named cofactor families and cutoff rules.")
    parser.add_argument("--output-dir", default="SSCNA_output", help="Output directory (default SSCNA_output).")
    parser.add_argument("--no-plots", action="store_true", help="Write the tables only; skip the interactive viewer.")
    parser.add_argument("--compact-html", action="store_true", help="Use a CDN-backed Plotly bundle for a smaller HTML file.")
    parser.add_argument("--external-residue-network", dest="external_residue_network_csv",
                        help="Overlay an external residue-pair CSV (source_chain, source_residue_number, target_chain, "
                             "target_residue_number) in the viewer.")
    show = parser.add_mutually_exclusive_group()
    show.add_argument("--show", action="store_true", help="Open the viewer in a browser tab (default: only from a terminal).")
    show.add_argument("--no-show", action="store_true", help="Never open a browser tab; files are still written.")
    parser.add_argument("--verbose", action="store_true", help="Show informational progress logs.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING, format="[%(levelname)s] %(message)s")

    path = resolve_structure_path(args.input)
    if not path:
        print(f"[ERROR] Structure file not found (tried PDB/CIF variants): {args.input}", file=sys.stderr)
        return 1
    cofactor, cofactor2 = split_names(args.cofactor), split_names(args.cofactor2) or None
    if not cofactor:
        print("[ERROR] At least one cofactor residue name is required.", file=sys.stderr)
        return 1

    class_sources = []
    if args.cofactor_class_config:
        class_sources.append(load_cofactor_class_config(args.cofactor_class_config))
    if args.cofactor_class_cutoff:
        class_sources.append(parse_class_cutoffs(args.cofactor_class_cutoff))
    class_cutoffs = merge_cofactor_class_configs(*class_sources) if class_sources else {}

    if args.show or args.no_show:
        from .display import set_show_figures

        set_show_figures(bool(args.show))

    from .coordination_api import analyze_structure

    tables = analyze_structure(
        path,
        cofactor,
        distance_cutoff=args.distance,
        expand_residues=args.expand_residues,
        combinatorial=args.combinatorial,
        combinatorial_cofactor_cutoff=args.combinatorial_cofactor_cutoff,
        cofactor_resname2=cofactor2,
        exclude_moieties=split_names(args.exclude_moieties),
        first_model_only=args.first_model,
        shells=args.shells,
        site_mode="per-site" if args.per_site else "union",
        site_model_mode=args.site_model_mode,
        include_carbon_seeds=args.include_carbon_seeds,
        direct_coordination=args.direct_coordination,
        direct_coordination_cutoff=args.direct_coordination_cutoff,
        cofactor_class_cutoffs=class_cutoffs,
    )

    os.makedirs(args.output_dir, exist_ok=True)
    name = os.path.basename(path)
    prefix = os.path.join(args.output_dir, f"{name}_")
    for table, suffix in (("residues", "Coordination_Residues"), ("atoms", "Coordination_Atoms"),
                          ("links", "Coord_Links"), ("contacts", "Coord_Contacts")):
        tables[table].to_csv(f"{prefix}{suffix}.csv", index=False)

    if not args.no_plots:
        from .moieties import atom_type_colors, bond_lookup
        from .plotting import plot_interactive_cohesive_network

        logger.info("Rendering cohesive interactive viewer...")
        structure, _ = unpack_pdb_file(path)
        groups = atom_groups_by_shell(tables)
        plot_interactive_cohesive_network(
            structure=structure,
            cofactor_atoms=groups.get("Cofactor", []),
            pcs_atoms=groups.get("PCS", []),
            scs_atoms=groups.get("SCS", []),
            focused_atoms_by_shell=groups,
            contacts=tables["contacts"],
            bond_lookup_table=bond_lookup,
            pdb_name=name,
            cofactor_resname=",".join(cofactor),
            atom_type_colors=atom_type_colors,
            output_filename=f"{prefix}coordination_network.html",
            first_model_only=args.first_model,
            compact_html=args.compact_html,
            external_residue_network_csv=args.external_residue_network_csv,
        )
        from .display import figures_should_show

        if figures_should_show():
            import webbrowser
            from pathlib import Path

            webbrowser.open(Path(f"{prefix}coordination_network.html").resolve().as_uri())
    logger.info("Wrote %d-shell coordination tables to %s.", args.shells, args.output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
