"""File-oriented workflows for coordination-network comparisons.

The lower-level comparison functions are intentionally pure and operate on
already-built signatures.  This module provides the small amount of orchestration
needed by a CLI or future web endpoint: discover structure files, analyze each
one, isolate per-file errors, and attach source paths to JSON/CSV-ready results.
"""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple

from .coordination_api import analyze_structure, discover_structure_paths
from .network_comparison import (
    build_network_profile,
    build_network_signature,
    compare_network_signatures,
    score_signature_against_profile,
)
from .network_alignment import align_analysis_tables, reference_position_map, reference_position_order


def _structure_id(path: Path) -> str:
    """Return a readable ID without common PDB/mmCIF suffixes."""
    name = path.name
    for _ in range(2):
        lowered = name.lower()
        for suffix in (".pdb.gz", ".mmcif.gz", ".cif.gz", ".pdb", ".mmcif", ".cif"):
            if lowered.endswith(suffix):
                name = name[: -len(suffix)]
                break
        else:
            break
    return name or path.stem


def analysis_from_path(
    structure_path: Path,
    cofactor_resname: str,
    *,
    analysis_options: Optional[Mapping[str, Any]] = None,
):
    """Analyze one structure without emitting progress output."""
    options = dict(analysis_options or {})
    # The analysis API has useful progress output for interactive use.  A
    # comparison runner should keep stdout reserved for its own summary.
    with redirect_stdout(StringIO()):
        return analyze_structure(structure_path, cofactor_resname, **options)


def signature_from_path(
    structure_path: Path,
    cofactor_resname: str,
    *,
    analysis_options: Optional[Mapping[str, Any]] = None,
    residue_map: Optional[Mapping[Any, Any]] = None,
) -> Dict[str, Any]:
    """Analyze one structure and return its canonical signature."""
    tables = analysis_from_path(
        structure_path,
        cofactor_resname,
        analysis_options=analysis_options,
    )
    return build_network_signature(
        tables,
        structure_id=_structure_id(structure_path),
        residue_map=residue_map,
    )


def _collect_analyses(
    inputs: Iterable[Path],
    cofactor_resname: str,
    *,
    analysis_options: Optional[Mapping[str, Any]] = None,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, str]]]:
    records: List[Dict[str, Any]] = []
    errors: List[Dict[str, str]] = []
    for path in discover_structure_paths(list(inputs)):
        try:
            analysis = analysis_from_path(
                path,
                cofactor_resname,
                analysis_options=analysis_options,
            )
        except Exception as exc:
            errors.append({
                "structure_id": _structure_id(path),
                "structure_path": str(path),
                "error": f"{type(exc).__name__}: {exc}",
            })
            continue
        records.append({"structure_path": str(path), "path": path, "analysis": analysis})
    return records, errors


def collect_signatures(
    inputs: Iterable[Path],
    cofactor_resname: str,
    *,
    analysis_options: Optional[Mapping[str, Any]] = None,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, str]]]:
    """Build signatures for files while isolating per-structure failures."""
    analysis_records, errors = _collect_analyses(
        inputs,
        cofactor_resname,
        analysis_options=analysis_options,
    )
    signatures = [
        {
            "structure_path": record["structure_path"],
            "signature": build_network_signature(
                record["analysis"],
                structure_id=_structure_id(record["path"]),
            ),
        }
        for record in analysis_records
    ]
    return signatures, errors


def compare_reference_to_queries(
    reference_path: Path,
    query_inputs: Iterable[Path],
    cofactor_resname: str,
    *,
    analysis_options: Optional[Mapping[str, Any]] = None,
    distance_tolerance_A: float = 0.5,
    align: bool = True,
    alignment_cutoff_A: float = 2.0,
) -> Dict[str, Any]:
    """Compare one reference structure against one or many query structures."""
    reference_analysis = analysis_from_path(
        reference_path,
        cofactor_resname,
        analysis_options=analysis_options,
    )
    query_analysis_records, errors = _collect_analyses(
        query_inputs,
        cofactor_resname,
        analysis_options=analysis_options,
    )
    alignments: Dict[str, Dict[str, Any]] = {}
    reference_map: Optional[Mapping[Any, Any]] = None
    if align:
        for record in query_analysis_records:
            alignment = align_analysis_tables(
                reference_analysis,
                record["analysis"],
                match_cutoff_A=alignment_cutoff_A,
            )
            alignments[record["structure_path"]] = alignment
            if alignment["metrics"]["status"] == "aligned" and reference_map is None:
                reference_map = alignment["reference_residue_map"]
    reference = build_network_signature(
        reference_analysis,
        structure_id=_structure_id(reference_path),
        residue_map=reference_map,
    )
    query_records: List[Dict[str, Any]] = []
    comparisons: List[Dict[str, Any]] = []
    for record in query_analysis_records:
        alignment = alignments.get(record["structure_path"])
        query_map = alignment["query_residue_map"] if alignment and alignment["metrics"]["status"] == "aligned" else None
        query_signature = build_network_signature(
            record["analysis"],
            structure_id=_structure_id(record["path"]),
            residue_map=query_map,
        )
        query_record = {
            "structure_path": record["structure_path"],
            "signature": query_signature,
        }
        if alignment:
            query_record["alignment"] = alignment["metrics"]
        query_records.append(query_record)
        result = compare_network_signatures(
            reference,
            query_signature,
            distance_tolerance_A=distance_tolerance_A,
        )
        result["reference_path"] = str(reference_path)
        result["query_path"] = record["structure_path"]
        if alignment:
            result["alignment"] = alignment["metrics"]
        comparisons.append(result)
    comparisons.sort(key=lambda row: row["scores"]["overall"], reverse=True)
    return {
        "mode": "pairwise",
        "reference": {"structure_path": str(reference_path), "signature": reference},
        "queries": query_records,
        "comparisons": comparisons,
        "errors": errors,
    }


def profile_reference_set(
    reference_inputs: Iterable[Path],
    cofactor_resname: str,
    *,
    template_path: Optional[Path] = None,
    query_inputs: Optional[Iterable[Path]] = None,
    analysis_options: Optional[Mapping[str, Any]] = None,
    required_support: float = 0.8,
    distance_tolerance_A: float = 0.5,
    align: bool = True,
    alignment_cutoff_A: float = 2.0,
) -> Dict[str, Any]:
    """Build a family profile and optionally score one or many query structures.

    ``template_path`` supplies the numbering/alignment coordinate system for
    the family map.  It may be one of the reference structures or a separate
    template structure used only as the geometric anchor.
    """
    reference_analysis_records, reference_errors = _collect_analyses(
        reference_inputs,
        cofactor_resname,
        analysis_options=analysis_options,
    )
    if not reference_analysis_records:
        raise ValueError("No reference structures produced a valid signature")
    anchor_record = reference_analysis_records[0]
    anchor_analysis = anchor_record["analysis"]
    anchor_path = Path(anchor_record["path"]).resolve()
    if template_path is not None:
        anchor_path = Path(template_path).expanduser().resolve()
        anchor_analysis = analysis_from_path(
            anchor_path,
            cofactor_resname,
            analysis_options=analysis_options,
        )
    reference_map = reference_position_map(anchor_analysis) if align else None
    template_order = reference_position_order(anchor_analysis) if align else []
    candidate_alignments: Dict[str, Dict[str, Any]] = {}
    if align:
        for record in reference_analysis_records:
            record_path = Path(record["path"]).resolve()
            if template_path is None and record_path == anchor_path:
                continue
            alignment = align_analysis_tables(
                anchor_analysis,
                record["analysis"],
                match_cutoff_A=alignment_cutoff_A,
            )
            candidate_alignments[record["structure_path"]] = alignment

    query_analysis_records: List[Dict[str, Any]] = []
    query_errors: List[Dict[str, str]] = []
    if query_inputs is not None:
        query_analysis_records, query_errors = _collect_analyses(
            query_inputs,
            cofactor_resname,
            analysis_options=analysis_options,
        )
        for record in query_analysis_records:
            alignment = align_analysis_tables(
                anchor_analysis,
                record["analysis"],
                match_cutoff_A=alignment_cutoff_A,
            ) if align else None
            if alignment:
                candidate_alignments[record["structure_path"]] = alignment

    reference_records: List[Dict[str, Any]] = []
    for record in reference_analysis_records:
        alignment = candidate_alignments.get(record["structure_path"])
        residue_map = None
        if Path(record["path"]).resolve() == anchor_path and align:
            residue_map = reference_map
        elif alignment and alignment["metrics"]["status"] == "aligned":
            residue_map = alignment["query_residue_map"]
        signature = build_network_signature(
            record["analysis"],
            structure_id=_structure_id(record["path"]),
            residue_map=residue_map,
        )
        reference_record = {"structure_path": record["structure_path"], "signature": signature}
        if alignment:
            reference_record["alignment"] = alignment["metrics"]
        reference_records.append(reference_record)

    template_signature = build_network_signature(
        anchor_analysis,
        structure_id=_structure_id(anchor_path),
        residue_map=reference_map,
    )
    template_record = {
        "structure_path": str(anchor_path),
        "signature": template_signature,
        "position_order": template_order,
    }

    profile = build_network_profile(
        [record["signature"] for record in reference_records],
        required_support=required_support,
    )
    query_records: List[Dict[str, Any]] = []
    scores: List[Dict[str, Any]] = []
    for record in query_analysis_records:
        alignment = candidate_alignments.get(record["structure_path"])
        query_map = alignment["query_residue_map"] if alignment and alignment["metrics"]["status"] == "aligned" else None
        query_signature = build_network_signature(
            record["analysis"],
            structure_id=_structure_id(record["path"]),
            residue_map=query_map,
        )
        query_record = {"structure_path": record["structure_path"], "signature": query_signature}
        if alignment:
            query_record["alignment"] = alignment["metrics"]
        query_records.append(query_record)
        result = score_signature_against_profile(
            query_signature,
            profile,
            distance_tolerance_A=distance_tolerance_A,
        )
        result["query_path"] = record["structure_path"]
        if alignment:
            result["alignment"] = alignment["metrics"]
        scores.append(result)
    scores.sort(key=lambda row: row["scores"]["overall"], reverse=True)
    return {
        "mode": "profile",
        "template": template_record,
        "references": reference_records,
        "profile": profile,
        "queries": query_records,
        "scores": scores,
        "errors": reference_errors + query_errors,
    }


__all__ = [
    "analysis_from_path",
    "collect_signatures",
    "compare_reference_to_queries",
    "profile_reference_set",
    "signature_from_path",
]
