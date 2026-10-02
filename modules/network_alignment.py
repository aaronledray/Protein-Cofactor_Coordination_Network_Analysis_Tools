"""Geometry-aware alignment helpers for coordination-network comparisons.

The legacy template/query workflow used Kabsch alignment followed by nearest
atom matching.  This module keeps that useful idea, but applies it to the
side-effect-free analysis tables and returns explicit diagnostics plus residue
position maps for canonical signatures.
"""

from __future__ import annotations

from collections import defaultdict
from math import sqrt
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if value != value:
            return ""
    except Exception:
        pass
    return str(value).strip()


def _records(table: Any) -> List[Dict[str, Any]]:
    if table is None:
        return []
    if hasattr(table, "to_dict"):
        return [dict(row) for row in table.to_dict("records")]
    return [dict(row) for row in table]


def _coordinates(row: Mapping[str, Any]) -> np.ndarray:
    return np.asarray(
        [
            float(row.get("x", 0.0)),
            float(row.get("y", 0.0)),
            float(row.get("z", 0.0)),
        ],
        dtype=float,
    )


def _shell(row: Mapping[str, Any]) -> str:
    return _text(row.get("shell")).upper()


def _residue_key(row: Mapping[str, Any]) -> Tuple[str, str, str, str]:
    """Use a residue-name-aware key so hetero/protein numbering cannot collide."""
    return (
        _text(row.get("residue_name")).upper(),
        _text(row.get("chain")),
        _text(row.get("residue_number")),
        _text(row.get("insertion_code")),
    )


def reference_position_map(analysis: Mapping[str, Any]) -> Dict[Tuple[str, str, str, str], str]:
    """Assign stable position labels to every residue in a reference analysis."""
    rows = _records(analysis.get("atoms"))
    residue_shells: Dict[Tuple[str, str, str, str], set] = defaultdict(set)
    for row in rows:
        residue_shells[_residue_key(row)].add(_shell(row))
    mapping: Dict[Tuple[str, str, str, str], str] = {}
    for key in sorted(residue_shells):
        residue, chain, number, insertion = key
        shell = "+".join(sorted(residue_shells[key]))
        mapping[key] = f"ref:{shell}:{residue}:{chain}:{number}:{insertion}"
    return mapping


def reference_position_order(analysis: Mapping[str, Any]) -> List[str]:
    """Return template positions in the order encountered in its atom table."""
    mapping = reference_position_map(analysis)
    order: List[str] = []
    seen = set()
    for row in _records(analysis.get("atoms")):
        position = mapping.get(_residue_key(row))
        if position and position not in seen:
            seen.add(position)
            order.append(position)
    return order


def _anchor_key(row: Mapping[str, Any], *, include_shell: bool) -> Tuple[str, ...]:
    fields = (
        _text(row.get("residue_name")).upper(),
        _text(row.get("atom_name")).upper(),
        _text(row.get("element")).upper(),
        _text(row.get("motif")).lower(),
    )
    return ((_shell(row),) + fields) if include_shell else fields


def _grouped_rows(rows: Sequence[Mapping[str, Any]], *, include_shell: bool) -> Dict[Tuple[str, ...], List[Mapping[str, Any]]]:
    groups: Dict[Tuple[str, ...], List[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[_anchor_key(row, include_shell=include_shell)].append(row)
    return groups


def _anchor_pairs(
    reference_rows: Sequence[Mapping[str, Any]],
    query_rows: Sequence[Mapping[str, Any]],
    *,
    include_shell: bool,
) -> List[Tuple[Mapping[str, Any], Mapping[str, Any]]]:
    reference_groups = _grouped_rows(reference_rows, include_shell=include_shell)
    query_groups = _grouped_rows(query_rows, include_shell=include_shell)
    pairs: List[Tuple[Mapping[str, Any], Mapping[str, Any]]] = []
    for key in sorted(set(reference_groups) & set(query_groups)):
        reference_group = sorted(
            reference_groups[key],
            key=lambda row: tuple(_coordinates(row).tolist()),
        )
        query_group = sorted(
            query_groups[key],
            key=lambda row: tuple(_coordinates(row).tolist()),
        )
        pairs.extend(zip(reference_group, query_group))
    return pairs


def _kabsch(reference_coordinates: np.ndarray, query_coordinates: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    reference_center = reference_coordinates.mean(axis=0)
    query_center = query_coordinates.mean(axis=0)
    reference_centered = reference_coordinates - reference_center
    query_centered = query_coordinates - query_center
    covariance = reference_centered.T @ query_centered
    u_matrix, _, v_transpose = np.linalg.svd(covariance)
    rotation = v_transpose.T @ u_matrix.T
    if np.linalg.det(rotation) < 0:
        v_transpose[-1, :] *= -1
        rotation = v_transpose.T @ u_matrix.T
    translation = query_center - rotation @ reference_center
    aligned = (rotation @ reference_coordinates.T).T + translation
    rmsd = float(sqrt(np.mean(np.sum((aligned - query_coordinates) ** 2, axis=1))))
    return rotation, translation, rmsd


def _match_network_atoms(
    reference_rows: Sequence[Mapping[str, Any]],
    query_rows: Sequence[Mapping[str, Any]],
    rotation: np.ndarray,
    translation: np.ndarray,
    cutoff_A: float,
) -> List[Tuple[Mapping[str, Any], Mapping[str, Any], float]]:
    if not reference_rows or not query_rows:
        return []
    transformed_reference = [
        (rotation @ _coordinates(row)) + translation for row in reference_rows
    ]
    candidate_pairs: List[Tuple[float, int, int]] = []
    for reference_index, reference_row in enumerate(reference_rows):
        reference_element = _text(reference_row.get("element")).upper()
        reference_shell = _shell(reference_row)
        for query_index, query_row in enumerate(query_rows):
            if reference_shell != _shell(query_row):
                continue
            query_element = _text(query_row.get("element")).upper()
            if reference_element and query_element and reference_element != query_element:
                continue
            distance = float(np.linalg.norm(transformed_reference[reference_index] - _coordinates(query_row)))
            if distance <= cutoff_A:
                candidate_pairs.append((distance, reference_index, query_index))
    candidate_pairs.sort(key=lambda item: (item[0], item[1], item[2]))
    used_reference = set()
    used_query = set()
    matches: List[Tuple[Mapping[str, Any], Mapping[str, Any], float]] = []
    for distance, reference_index, query_index in candidate_pairs:
        if reference_index in used_reference or query_index in used_query:
            continue
        used_reference.add(reference_index)
        used_query.add(query_index)
        matches.append((reference_rows[reference_index], query_rows[query_index], distance))
    return matches


def align_analysis_tables(
    reference_analysis: Mapping[str, Any],
    query_analysis: Mapping[str, Any],
    *,
    match_cutoff_A: float = 2.0,
) -> Dict[str, Any]:
    """Align two analysis results and produce residue maps for signatures.

    Cofactor atoms are preferred as anchors. If fewer than three shared
    cofactor anchors exist, shared network atoms are attempted. With fewer than
    three usable anchors no rigid alignment is claimed and the returned maps
    are empty, leaving the chemistry/topology comparison available.
    """
    reference_atoms = _records(reference_analysis.get("atoms"))
    query_atoms = _records(query_analysis.get("atoms"))
    reference_cofactor = [row for row in reference_atoms if _shell(row) == "COFACTOR"]
    query_cofactor = [row for row in query_atoms if _shell(row) == "COFACTOR"]
    anchor_pairs = _anchor_pairs(reference_cofactor, query_cofactor, include_shell=False)
    method = "cofactor_kabsch"
    if len(anchor_pairs) < 3:
        anchor_pairs = _anchor_pairs(reference_atoms, query_atoms, include_shell=True)
        method = "network_kabsch"

    if len(anchor_pairs) < 3:
        return {
            "reference_residue_map": {},
            "query_residue_map": {},
            "metrics": {
                "status": "insufficient_anchors",
                "method": "none",
                "anchor_count": len(anchor_pairs),
                "anchor_rmsd_A": None,
                "matched_atom_count": 0,
                "residue_mapping_count": 0,
                "reference_atom_count": len(reference_atoms),
                "query_atom_count": len(query_atoms),
            },
        }

    reference_coordinates = np.vstack([_coordinates(reference) for reference, _ in anchor_pairs])
    query_coordinates = np.vstack([_coordinates(query) for _, query in anchor_pairs])
    rotation, translation, anchor_rmsd = _kabsch(reference_coordinates, query_coordinates)
    matches = _match_network_atoms(
        reference_atoms,
        query_atoms,
        rotation,
        translation,
        match_cutoff_A,
    )

    reference_map = reference_position_map(reference_analysis)
    query_map: Dict[Tuple[str, str, str, str], str] = {}
    residue_match_counts: Dict[str, int] = defaultdict(int)
    residue_distances: Dict[str, List[float]] = defaultdict(list)
    for reference, query, distance in matches:
        reference_key = _residue_key(reference)
        query_key = _residue_key(query)
        position = reference_map.get(reference_key)
        if position is None:
            continue
        existing = query_map.get(query_key)
        if existing is not None and existing != position:
            continue
        query_map[query_key] = position
        residue_match_counts[position] += 1
        residue_distances[position].append(distance)

    residue_matches = [
        {
            "position": position,
            "atom_match_count": residue_match_counts[position],
            "mean_distance_A": round(float(np.mean(residue_distances[position])), 6),
        }
        for position in sorted(residue_match_counts)
    ]
    return {
        "reference_residue_map": reference_map,
        "query_residue_map": query_map,
        "metrics": {
            "status": "aligned",
            "method": method,
            "anchor_count": len(anchor_pairs),
            "anchor_rmsd_A": round(anchor_rmsd, 6),
            "matched_atom_count": len(matches),
            "residue_mapping_count": len(query_map),
            "reference_atom_count": len(reference_atoms),
            "query_atom_count": len(query_atoms),
            "residue_matches": residue_matches,
        },
    }


__all__ = ["align_analysis_tables", "reference_position_map", "reference_position_order"]
