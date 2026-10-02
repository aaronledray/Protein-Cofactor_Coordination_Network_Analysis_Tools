#!/usr/bin/env python3
"""Compare one or many canonical coordination-network signatures.

Examples
--------
Pairwise template/query comparison::

    python compare_coordination_networks.py pairwise \
        --reference reference.pdb --query query_a.pdb query_b.pdb \
        --cofactor CU --output-dir comparison_output

Family profile and query scoring::

    python compare_coordination_networks.py profile \
        --reference family/*.pdb --query candidates/*.pdb \
        --cofactor OEX --output-dir profile_output
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from modules.comparison_runner import compare_reference_to_queries, profile_reference_set
from modules.conservation_viewer import write_profile_conservation_viewer
from modules.coordination_api import discover_structure_paths
from modules.comparison_plots import (
    write_residue_conservation_heatmap,
    write_similarity_heatmap,
)
from modules.cofactor_classes import load_cofactor_class_config, merge_cofactor_class_configs


def _add_analysis_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cofactor", required=True, help="Cofactor residue name(s), comma-separated.")
    parser.add_argument("--cofactor2", help="Optional second cofactor residue name(s), comma-separated.")
    parser.add_argument("--distance", type=float, default=3.6, help="Shell/contact cutoff in Å.")
    parser.add_argument("--shells", type=int, default=3, help="Number of coordination shells to retain.")
    parser.add_argument(
        "--site-model-mode",
        choices=["pooled", "per-model"],
        default="pooled",
        help="Boundary policy for per-site analyses across structure models.",
    )
    parser.add_argument(
        "--alignment-cutoff",
        type=float,
        default=2.0,
        help="Maximum post-alignment atom distance used for residue mapping in Å.",
    )
    parser.add_argument(
        "--no-alignment",
        dest="align",
        action="store_false",
        default=True,
        help="Disable Kabsch/nearest-atom alignment and use chemistry/topology only.",
    )
    parser.add_argument("--include-carbon-seeds", action="store_true")
    parser.add_argument("--direct-coordination", action="store_true")
    parser.add_argument("--direct-coordination-cutoff", type=float, default=2.6)
    parser.add_argument(
        "--cofactor-class-cutoff",
        action="append",
        default=[],
        help="Class-specific cutoff, e.g. metal=2.8; repeatable.",
    )
    parser.add_argument(
        "--cofactor-class-config",
        help="YAML/JSON file defining named cofactor families and cutoff rules.",
    )
    parser.add_argument("--expand-residues", action="store_true")
    parser.add_argument("--first-model", action="store_true")
    parser.add_argument("--combinatorial", action="store_true")
    parser.add_argument("--combinatorial-cutoff", type=float, default=20.0)
    parser.add_argument(
        "--exclude-moieties",
        default="alanine_sidechain",
        help="Comma-separated motif exclusions; use an empty string to disable.",
    )
    parser.add_argument(
        "--compact-html",
        action="store_true",
        help="Use a CDN-backed Plotly bundle for the conservation viewer.",
    )


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and compare motif-aware coordination-network signatures."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    pairwise = subparsers.add_parser(
        "pairwise",
        help="Compare one reference structure against one or many queries.",
    )
    pairwise.add_argument("--reference", required=True, help="One reference PDB/mmCIF file.")
    pairwise.add_argument("--query", required=True, nargs="+", help="Query files or directories.")
    pairwise.add_argument("--output-dir", default="network_comparison_output")
    pairwise.add_argument("--distance-tolerance", type=float, default=0.5)
    _add_analysis_arguments(pairwise)

    profile = subparsers.add_parser(
        "profile",
        help="Build a family profile and optionally score query structures.",
    )
    profile.add_argument("--reference", required=True, nargs="+", help="Known family files or directories.")
    profile.add_argument(
        "--template",
        help="Optional template structure whose residue numbering anchors alignment and conservation-map columns.",
    )
    profile.add_argument("--query", nargs="*", help="Optional query files or directories.")
    profile.add_argument("--output-dir", default="network_profile_output")
    profile.add_argument("--required-support", type=float, default=0.8)
    profile.add_argument("--distance-tolerance", type=float, default=0.5)
    _add_analysis_arguments(profile)
    return parser.parse_args(argv)


def _analysis_options(args: argparse.Namespace) -> Dict[str, Any]:
    class_cutoffs: Dict[str, Any] = {}
    for value in args.cofactor_class_cutoff:
        if "=" not in value:
            raise ValueError("--cofactor-class-cutoff must use CLASS=ANGSTROMS")
        name, cutoff = value.split("=", 1)
        class_cutoffs[name.strip().lower()] = float(cutoff)
    config_sources = []
    if args.cofactor_class_config:
        config_sources.append(load_cofactor_class_config(args.cofactor_class_config))
    if class_cutoffs:
        config_sources.append(class_cutoffs)
    return {
        "distance_cutoff": args.distance,
        "shells": args.shells,
        "site_model_mode": args.site_model_mode,
        "include_carbon_seeds": args.include_carbon_seeds,
        "direct_coordination": args.direct_coordination,
        "direct_coordination_cutoff": args.direct_coordination_cutoff,
        "cofactor_class_cutoffs": merge_cofactor_class_configs(*config_sources),
        "expand_residues": args.expand_residues,
        "first_model_only": args.first_model,
        "combinatorial": args.combinatorial,
        "combinatorial_cofactor_cutoff": args.combinatorial_cutoff,
        "cofactor_resname2": args.cofactor2,
        "exclude_moieties": [
            value.strip() for value in args.exclude_moieties.split(",") if value.strip()
        ],
    }


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _safe_filename(value: str) -> str:
    cleaned = "".join(character if character.isalnum() or character in "-_." else "_" for character in value)
    return cleaned or "structure"


def _write_signatures(output_dir: Path, records: Iterable[Mapping[str, Any]], prefix: str = "signature") -> None:
    signature_dir = output_dir / "signatures"
    signature_dir.mkdir(parents=True, exist_ok=True)
    for index, record in enumerate(records, start=1):
        signature = record["signature"]
        structure_id = str(signature.get("structure_id") or f"structure_{index}")
        path = signature_dir / f"{index:03d}_{_safe_filename(structure_id)}_{prefix}.json"
        _write_json(path, signature)


def _write_errors(path: Path, errors: Iterable[Mapping[str, Any]]) -> None:
    rows = list(errors)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["structure_id", "structure_path", "error"])
        writer.writeheader()
        writer.writerows(rows)


def _write_pairwise_artifacts(output_dir: Path, result: Mapping[str, Any]) -> None:
    _write_json(output_dir / "pairwise_comparison.json", result)
    _write_signatures(output_dir, [result["reference"]], prefix="reference")
    _write_signatures(output_dir, result["queries"], prefix="query")
    rows: List[Dict[str, Any]] = []
    for comparison in result["comparisons"]:
        scores = comparison["scores"]
        alignment = comparison.get("alignment", {})
        rows.append({
            "reference_id": comparison["reference_id"],
            "query_id": comparison["query_id"],
            "reference_path": comparison["reference_path"],
            "query_path": comparison["query_path"],
            "cofactor_compatible": comparison["cofactor_compatible"],
            "alignment_status": alignment.get("status", "not_requested"),
            "alignment_method": alignment.get("method", "none"),
            "anchor_count": alignment.get("anchor_count", 0),
            "anchor_rmsd_A": alignment.get("anchor_rmsd_A"),
            "matched_atom_count": alignment.get("matched_atom_count", 0),
            "residue_mapping_count": alignment.get("residue_mapping_count", 0),
            **scores,
        })
    with (output_dir / "pairwise_comparisons.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "reference_id", "query_id", "reference_path", "query_path",
            "cofactor_compatible", "overall", "edge", "residue", "distance",
            "primary", "secondary", "tertiary", "alignment_status",
            "alignment_method", "anchor_count", "anchor_rmsd_A",
            "matched_atom_count", "residue_mapping_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    _write_errors(output_dir / "errors.csv", result["errors"])
    heatmap = write_similarity_heatmap(
        [result["reference"], *result["queries"]],
        output_dir,
        stem="pairwise_similarity_heatmap",
        title="Pairwise coordination-network similarity",
    )
    _write_json(output_dir / "pairwise_similarity_heatmap_details.json", heatmap)


def _write_profile_artifacts(output_dir: Path, result: Mapping[str, Any]) -> Dict[str, Any]:
    _write_json(output_dir / "network_profile.json", result["profile"])
    _write_json(output_dir / "profile_comparison.json", result)
    _write_signatures(output_dir, result["references"], prefix="reference")
    _write_signatures(output_dir, result["queries"], prefix="query")
    rows: List[Dict[str, Any]] = []
    for score in result["scores"]:
        alignment = score.get("alignment", {})
        rows.append({
            "query_id": score["query_id"],
            "query_path": score["query_path"],
            "cofactor_compatible": score["cofactor_compatible"],
            "alignment_status": alignment.get("status", "not_requested"),
            "alignment_method": alignment.get("method", "none"),
            "anchor_count": alignment.get("anchor_count", 0),
            "anchor_rmsd_A": alignment.get("anchor_rmsd_A"),
            "matched_atom_count": alignment.get("matched_atom_count", 0),
            "residue_mapping_count": alignment.get("residue_mapping_count", 0),
            **score["scores"],
            "primary_coverage": score["edges"]["primary"]["coverage"],
            "secondary_coverage": score["edges"]["secondary"]["coverage"],
            "tertiary_coverage": score["edges"]["tertiary"]["coverage"],
            "residue_coverage": score["residues"]["coverage"],
        })
    with (output_dir / "profile_scores.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "query_id", "query_path", "cofactor_compatible", "overall", "edge",
            "residue", "distance", "primary_coverage", "secondary_coverage",
            "tertiary_coverage", "residue_coverage",
            "alignment_status", "alignment_method", "anchor_count",
            "anchor_rmsd_A", "matched_atom_count", "residue_mapping_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    _write_errors(output_dir / "errors.csv", result["errors"])
    heatmap = write_similarity_heatmap(
        result["references"],
        output_dir,
        stem="reference_similarity_heatmap",
        title="Reference coordination-network similarity",
    )
    _write_json(output_dir / "reference_similarity_heatmap_details.json", heatmap)
    try:
        conservation_map = write_residue_conservation_heatmap(
            result,
            output_dir,
            stem="residue_conservation_map",
            title="Family residue conservation",
        )
    except ValueError as exc:
        conservation_map = {
            "status": "unavailable",
            "reason": str(exc),
        }
    _write_json(output_dir / "residue_conservation_map_details.json", conservation_map)
    return conservation_map


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    options = _analysis_options(args)
    if args.command == "pairwise":
        reference_paths = list(discover_structure_paths([args.reference]))
        if len(reference_paths) != 1:
            raise SystemExit("--reference must resolve to exactly one structure file")
        result = compare_reference_to_queries(
            reference_paths[0],
            args.query,
            args.cofactor,
            analysis_options=options,
            distance_tolerance_A=args.distance_tolerance,
            align=args.align,
            alignment_cutoff_A=args.alignment_cutoff,
        )
        _write_pairwise_artifacts(output_dir, result)
        print(f"Compared {len(result['comparisons'])} query structure(s); errors={len(result['errors'])}")
        print(f"Output: {output_dir}")
        return 0

    template_paths = discover_structure_paths([args.template]) if args.template else []
    if len(template_paths) > 1:
        raise SystemExit("--template must resolve to exactly one structure file")
    template_path = template_paths[0] if template_paths else None
    if template_path is not None and not template_path.is_file():
        raise SystemExit(f"Template structure not found: {template_path}")
    result = profile_reference_set(
        args.reference,
        args.cofactor,
        template_path=template_path,
        query_inputs=args.query or None,
        analysis_options=options,
        required_support=args.required_support,
        distance_tolerance_A=args.distance_tolerance,
        align=args.align,
        alignment_cutoff_A=args.alignment_cutoff,
    )
    conservation_map = _write_profile_artifacts(output_dir, result)
    conservation_viewer = write_profile_conservation_viewer(
        Path(result["references"][0]["structure_path"]),
        args.cofactor,
        result,
        output_dir / "profile_conservation_viewer.html",
        analysis_options=options,
        compact_html=args.compact_html,
    )
    print(
        f"Built profile from {len(result['references'])} reference structure(s); "
        f"scored {len(result['scores'])} query structure(s); errors={len(result['errors'])}"
    )
    print(f"Output: {output_dir}")
    print(f"Conservation viewer: {conservation_viewer}")
    if conservation_map.get("html"):
        print(f"2D conservation map: {conservation_map['html']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
