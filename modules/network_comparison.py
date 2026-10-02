"""Canonical coordination-network signatures and transparent comparisons.

This module deliberately sits below any web UI.  ``analyze_structure`` produces
structure-specific tables; the functions here turn those tables into stable,
JSON-friendly feature records that can be compared directly or aggregated into
a family profile.

The first scorer is intentionally explainable rather than a black-box model:
it compares rooted contact topology, chemical/motif identity, residue features,
and contact distances.  Sequence/structure alignment can be supplied later by
mapping residue identities to shared ``position`` labels.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from math import exp, sqrt
from statistics import mean, median, pstdev
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SIGNATURE_SCHEMA_VERSION = "1.0"
CONTACT_CATEGORIES = ("primary", "secondary", "tertiary")
DEFAULT_EDGE_WEIGHTS = {"primary": 0.55, "secondary": 0.25, "tertiary": 0.20}


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if value != value:  # NaN
            return ""
    except Exception:
        pass
    return str(value).strip()


def _bool(value: Any) -> bool:
    return value is True or _text(value).lower() in {"true", "1", "yes"}


def _float(value: Any) -> Optional[float]:
    try:
        if value is None:
            return None
        number = float(value)
        return number if number == number else None
    except (TypeError, ValueError):
        return None


def _json_value(value: Any) -> Any:
    """Convert common NumPy/Pandas scalar values to JSON-friendly values."""
    if value is None:
        return None
    if hasattr(value, "item"):
        try:
            return value.item()
        except ValueError:
            pass
    return value


def _records(table: Any) -> List[Dict[str, Any]]:
    if table is None:
        return []
    if hasattr(table, "to_dict"):
        return [dict(row) for row in table.to_dict("records")]
    return [dict(row) for row in table]


def _row_position(
    row: Mapping[str, Any],
    prefix: str,
    residue_map: Optional[Mapping[Any, Any]],
) -> Optional[str]:
    if not residue_map:
        return None
    chain = _text(row.get(f"{prefix}_chain"))
    number = _text(row.get(f"{prefix}_resnum"))
    insertion = _text(row.get(f"{prefix}_insertion_code"))
    residue = _text(row.get(f"{prefix}_resname")).upper()
    candidates = [
        (chain, number, insertion),
        (residue, chain, number, insertion),
        f"{chain}:{number}{insertion}",
        f"{residue}:{chain}:{number}{insertion}",
    ]
    for key in candidates:
        if key in residue_map:
            value = residue_map[key]
            return None if value is None else _text(value)
    return None


def _atom_token(
    row: Mapping[str, Any],
    prefix: str,
    residue_map: Optional[Mapping[Any, Any]],
) -> Dict[str, Any]:
    shell = _text(row.get(f"{prefix}_shell")).upper()
    residue = _text(row.get(f"{prefix}_resname")).upper()
    atom = _text(row.get(f"{prefix}_atom")).upper()
    element = _text(row.get(f"{prefix}_element")).upper()
    motif = _text(row.get(f"{prefix}_motif")).lower()
    position = _row_position(row, prefix, residue_map)
    token = {
        "shell": shell,
        "residue": residue,
        "atom": atom,
        "element": element,
        "motif": motif,
    }
    if position:
        token["position"] = position
    return token


def _token_key(token: Mapping[str, Any]) -> str:
    """Stable comparison key; numbering is omitted unless alignment supplied it."""
    fields = (
        _text(token.get("shell")).upper(),
        _text(token.get("residue")).upper(),
        _text(token.get("position")),
        _text(token.get("motif")).lower(),
        _text(token.get("atom")).upper(),
        _text(token.get("element")).upper(),
    )
    return "|".join(fields)


def _contact_category(row: Mapping[str, Any]) -> Optional[str]:
    role = _text(row.get("contact_role")).lower()
    link = _text(row.get("link_type")).lower()
    if _bool(row.get("direct_coordination")) or role in {
        "primary_motif_contact",
        "legacy_primary_contact",
    }:
        return "primary"
    if link == "pcs->scs" or _text(row.get("dst_shell")).upper() in {"SCS", "SHELL2"}:
        return "secondary"
    if link == "scs->tcs" or _text(row.get("dst_shell")).upper() in {"TCS", "SHELL3"}:
        return "tertiary"
    return None


def _edge_record(
    row: Mapping[str, Any],
    category: str,
    residue_map: Optional[Mapping[Any, Any]],
) -> Dict[str, Any]:
    source = _atom_token(row, "src", residue_map)
    destination = _atom_token(row, "dst", residue_map)
    distance = _float(row.get("distance_A"))
    record: Dict[str, Any] = {
        "category": category,
        "source": source,
        "destination": destination,
        "match_key": f"{category}:{_token_key(source)}->{_token_key(destination)}",
    }
    if distance is not None:
        record["distance_A"] = distance
    if _text(row.get("contact_role")):
        record["contact_role"] = _text(row.get("contact_role"))
    return record


def _atom_match_key(row: Mapping[str, Any], residue_map: Optional[Mapping[Any, Any]]) -> str:
    shell = _text(row.get("shell")).upper()
    residue = _text(row.get("residue_name")).upper()
    atom = _text(row.get("atom_name")).upper()
    element = _text(row.get("element")).upper()
    motif = _text(row.get("motif")).lower()
    chain = _text(row.get("chain"))
    number = _text(row.get("residue_number"))
    insertion = _text(row.get("insertion_code"))
    position = None
    if residue_map:
        for key in ((chain, number, insertion), (residue, chain, number, insertion)):
            if key in residue_map:
                position = _text(residue_map[key])
                break
    return "|".join((shell, residue, position or "", motif, atom, element))


def _atom_position(
    row: Mapping[str, Any],
    residue_map: Optional[Mapping[Any, Any]],
) -> Optional[str]:
    if not residue_map:
        return None
    residue = _text(row.get("residue_name")).upper()
    chain = _text(row.get("chain"))
    number = _text(row.get("residue_number"))
    insertion = _text(row.get("insertion_code"))
    for key in (
        (residue, chain, number, insertion),
        (chain, number, insertion),
    ):
        if key in residue_map:
            return _text(residue_map[key])
    return None


def build_network_signature(
    analysis: Mapping[str, Any],
    *,
    structure_id: Optional[str] = None,
    site_id: Optional[str] = None,
    residue_map: Optional[Mapping[Any, Any]] = None,
) -> Dict[str, Any]:
    """Build a canonical, JSON-friendly signature from ``analyze_structure`` output.

    ``residue_map`` is optional.  When supplied, use keys such as
    ``(chain, residue_number, insertion_code)`` and values such as an aligned
    reference position.  Without it, comparison intentionally ignores raw
    residue numbering so homologous structures can still be compared by
    chemistry and motif identity.
    """
    atoms = _records(analysis.get("atoms"))
    residues = _records(analysis.get("residues"))
    contacts = _records(analysis.get("contacts"))
    links = _records(analysis.get("links"))

    if structure_id is None:
        structure_id = _text(atoms[0].get("structure_id")) if atoms else ""
    if site_id is None:
        site_id = _text(atoms[0].get("site_id")) if atoms else "all"

    cofactor_names = sorted({
        _text(row.get("residue_name")).upper()
        for row in atoms
        if _text(row.get("shell")).upper() == "COFACTOR"
    })

    edge_groups: Dict[str, List[Dict[str, Any]]] = {category: [] for category in CONTACT_CATEGORIES}
    for row in contacts:
        category = _contact_category(row)
        if category:
            edge_groups[category].append(_edge_record(row, category, residue_map))
    for category in edge_groups:
        edge_groups[category].sort(key=lambda edge: edge["match_key"])

    atom_features: List[Dict[str, Any]] = []
    for row in atoms:
        feature = {
            "shell": _text(row.get("shell")).upper(),
            "residue": _text(row.get("residue_name")).upper(),
            "atom": _text(row.get("atom_name")).upper(),
            "element": _text(row.get("element")).upper(),
            "motif": _text(row.get("motif")).lower(),
            "coordination_role": _text(row.get("coordination_role")).lower(),
        }
        position = _atom_position(row, residue_map)
        if position:
            feature["position"] = position
        feature["match_key"] = _atom_match_key(row, residue_map)
        atom_features.append(feature)
    atom_features.sort(key=lambda feature: feature["match_key"])

    residue_groups: Dict[Tuple[str, str, str, str, str], Dict[str, Any]] = {}
    for row in atoms:
        key = (
            _text(row.get("shell")).upper(),
            _text(row.get("residue_name")).upper(),
            _text(row.get("chain")),
            _text(row.get("residue_number")),
            _text(row.get("insertion_code")),
        )
        group = residue_groups.setdefault(key, {
            "shell": key[0], "residue": key[1], "chain": key[2],
            "residue_number": key[3], "insertion_code": key[4],
            "atoms": set(), "motifs": set(), "roles": set(),
        })
        group["atoms"].add(_text(row.get("atom_name")).upper())
        group["motifs"].add(_text(row.get("motif")).lower())
        group["roles"].add(_text(row.get("coordination_role")).lower())

    residue_features: List[Dict[str, Any]] = []
    for group in residue_groups.values():
        position = None
        if residue_map:
            for key in (
                (group["chain"], group["residue_number"], group["insertion_code"]),
                (group["residue"], group["chain"], group["residue_number"], group["insertion_code"]),
            ):
                if key in residue_map:
                    position = _text(residue_map[key])
                    break
        feature = {
            "shell": group["shell"],
            "residue": group["residue"],
            "atoms": sorted(group["atoms"]),
            "motifs": sorted(group["motifs"]),
            "roles": sorted(group["roles"]),
        }
        if position:
            feature["position"] = position
        feature["match_key"] = "|".join((
            group["shell"], group["residue"], position or "",
            ",".join(feature["motifs"]), ",".join(feature["roles"]),
        ))
        residue_features.append(feature)
    residue_features.sort(key=lambda feature: feature["match_key"])

    return {
        "schema_version": SIGNATURE_SCHEMA_VERSION,
        "structure_id": _text(structure_id),
        "site_id": _text(site_id) or "all",
        "cofactor": {
            "residue_names": cofactor_names,
            "key": "+".join(cofactor_names),
        },
        "edges": edge_groups,
        "atoms": atom_features,
        "residues": residue_features,
        "summary": {
            "atom_count": len(atom_features),
            "residue_count": len(residue_features),
            "contact_count": sum(len(edges) for edges in edge_groups.values()),
            "link_count": len(links),
            "shell_counts": dict(Counter(feature["shell"] for feature in atom_features)),
        },
    }


def _multiset_jaccard(left: Sequence[str], right: Sequence[str]) -> float:
    left_counts, right_counts = Counter(left), Counter(right)
    intersection = sum((left_counts & right_counts).values())
    union = sum((left_counts | right_counts).values())
    return intersection / union if union else 1.0


def _edge_matches(reference_edges: Sequence[Mapping[str, Any]], query_edges: Sequence[Mapping[str, Any]]) -> Tuple[List[Tuple[Mapping[str, Any], Mapping[str, Any]]], List[Mapping[str, Any]], List[Mapping[str, Any]]]:
    reference_by_key: Dict[str, List[Mapping[str, Any]]] = defaultdict(list)
    query_by_key: Dict[str, List[Mapping[str, Any]]] = defaultdict(list)
    for edge in reference_edges:
        reference_by_key[_text(edge.get("match_key"))].append(edge)
    for edge in query_edges:
        query_by_key[_text(edge.get("match_key"))].append(edge)
    matches: List[Tuple[Mapping[str, Any], Mapping[str, Any]]] = []
    unmatched_reference: List[Mapping[str, Any]] = []
    unmatched_query: List[Mapping[str, Any]] = []
    for key in sorted(set(reference_by_key) | set(query_by_key)):
        ref_group, query_group = reference_by_key.get(key, []), query_by_key.get(key, [])
        matched_count = min(len(ref_group), len(query_group))
        for ref_edge, query_edge in zip(ref_group[:matched_count], query_group[:matched_count]):
            matches.append((ref_edge, query_edge))
        unmatched_reference.extend(ref_group[matched_count:])
        unmatched_query.extend(query_group[matched_count:])
    return matches, unmatched_reference, unmatched_query


def _distance_score(matches: Sequence[Tuple[Mapping[str, Any], Mapping[str, Any]]], tolerance_A: float) -> Optional[float]:
    scores = []
    for reference, query in matches:
        left, right = _float(reference.get("distance_A")), _float(query.get("distance_A"))
        if left is not None and right is not None:
            scores.append(exp(-abs(left - right) / max(tolerance_A, 1e-9)))
    return mean(scores) if scores else None


def compare_network_signatures(
    reference: Mapping[str, Any],
    query: Mapping[str, Any],
    *,
    distance_tolerance_A: float = 0.5,
    edge_weights: Optional[Mapping[str, float]] = None,
) -> Dict[str, Any]:
    """Compare two signatures and return scores plus explainable mismatches."""
    weights = dict(DEFAULT_EDGE_WEIGHTS)
    if edge_weights:
        weights.update({key: float(value) for key, value in edge_weights.items()})

    reference_cofactor = _text(reference.get("cofactor", {}).get("key"))
    query_cofactor = _text(query.get("cofactor", {}).get("key"))
    cofactor_compatible = reference_cofactor == query_cofactor

    category_scores: Dict[str, float] = {}
    category_distances: Dict[str, Optional[float]] = {}
    unmatched: Dict[str, Dict[str, List[Mapping[str, Any]]]] = {}
    all_matches: List[Tuple[Mapping[str, Any], Mapping[str, Any]]] = []
    active_weight = 0.0
    active_categories: List[str] = []
    for category in CONTACT_CATEGORIES:
        reference_edges = list(reference.get("edges", {}).get(category, []))
        query_edges = list(query.get("edges", {}).get(category, []))
        matches, missing, extra = _edge_matches(reference_edges, query_edges)
        category_scores[category] = _multiset_jaccard(
            [_text(edge.get("match_key")) for edge in reference_edges],
            [_text(edge.get("match_key")) for edge in query_edges],
        )
        category_distances[category] = _distance_score(matches, distance_tolerance_A)
        unmatched[category] = {"reference_only": missing, "query_only": extra}
        all_matches.extend(matches)
        if reference_edges or query_edges:
            active_weight += weights.get(category, 0.0)
            active_categories.append(category)

    edge_score = (
        sum(weights.get(category, 0.0) * category_scores[category] for category in active_categories)
        / active_weight
        if active_weight else 1.0
    )
    distances = [score for score in category_distances.values() if score is not None]
    distance_score = mean(distances) if distances else None
    residue_score = _multiset_jaccard(
        [_text(feature.get("match_key")) for feature in reference.get("residues", [])],
        [_text(feature.get("match_key")) for feature in query.get("residues", [])],
    )
    components = [(edge_score, 0.65), (residue_score, 0.15)]
    if distance_score is not None:
        components.append((distance_score, 0.20))
    denominator = sum(weight for _, weight in components)
    overall = sum(score * weight for score, weight in components) / denominator
    if not cofactor_compatible:
        overall *= 0.25

    return {
        "schema_version": SIGNATURE_SCHEMA_VERSION,
        "reference_id": _text(reference.get("structure_id")),
        "query_id": _text(query.get("structure_id")),
        "cofactor_compatible": cofactor_compatible,
        "scores": {
            "overall": round(float(overall), 6),
            "edge": round(float(edge_score), 6),
            "residue": round(float(residue_score), 6),
            "distance": None if distance_score is None else round(float(distance_score), 6),
            "primary": round(float(category_scores["primary"]), 6),
            "secondary": round(float(category_scores["secondary"]), 6),
            "tertiary": round(float(category_scores["tertiary"]), 6),
        },
        "matched_edge_count": len(all_matches),
        "unmatched": unmatched,
    }


def _distance_summary(distances: Sequence[float]) -> Dict[str, Optional[float]]:
    if not distances:
        return {"count": 0, "mean_A": None, "median_A": None, "std_A": None, "min_A": None, "max_A": None}
    return {
        "count": len(distances),
        "mean_A": round(mean(distances), 6),
        "median_A": round(median(distances), 6),
        "std_A": round(pstdev(distances), 6) if len(distances) > 1 else 0.0,
        "min_A": round(min(distances), 6),
        "max_A": round(max(distances), 6),
    }


def _count_summary(counts: Sequence[int]) -> Dict[str, Optional[float]]:
    """Summarize per-signature feature multiplicity for profile scoring."""
    if not counts:
        return {"count": 0, "mean": None, "median": None, "min": None, "max": None}
    return {
        "count": len(counts),
        "mean": round(mean(counts), 6),
        "median": round(median(counts), 6),
        "min": min(counts),
        "max": max(counts),
    }


def build_network_profile(
    signatures: Iterable[Mapping[str, Any]],
    *,
    required_support: float = 0.8,
) -> Dict[str, Any]:
    """Aggregate aligned or alignment-free signatures into a family profile."""
    signatures = list(signatures)
    if not signatures:
        raise ValueError("At least one signature is required to build a profile")
    if not 0.0 < required_support <= 1.0:
        raise ValueError("required_support must be in the interval (0, 1]")

    n = len(signatures)
    cofactor_counts = Counter(_text(signature.get("cofactor", {}).get("key")) for signature in signatures)
    profile_edges: Dict[str, Dict[str, Dict[str, Any]]] = {category: {} for category in CONTACT_CATEGORIES}
    for category in CONTACT_CATEGORIES:
        for signature in signatures:
            signature_edges = signature.get("edges", {}).get(category, [])
            signature_counts = Counter(_text(edge.get("match_key")) for edge in signature_edges)
            for key, occurrence_count in signature_counts.items():
                entry = profile_edges[category].setdefault(
                    key,
                    {
                        "feature_key": key,
                        "count": 0,
                        "signature_count": 0,
                        "multiplicities": [],
                        "distances_A": [],
                    },
                )
                entry["count"] += occurrence_count
                entry["signature_count"] += 1
                entry["multiplicities"].append(occurrence_count)
            for edge in signature_edges:
                key = _text(edge.get("match_key"))
                distance = _float(edge.get("distance_A"))
                if distance is not None:
                    profile_edges[category][key]["distances_A"].append(distance)
        for key, entry in profile_edges[category].items():
            entry["support_count"] = entry.pop("signature_count")
            entry["support"] = round(entry["support_count"] / n, 6)
            entry["required"] = entry["support"] >= required_support
            entry["multiplicity"] = _count_summary(entry.pop("multiplicities"))
            entry["distance"] = _distance_summary(entry.pop("distances_A"))
        profile_edges[category] = dict(sorted(profile_edges[category].items()))

    profile_residues: Dict[str, Dict[str, Any]] = {}
    for signature in signatures:
        signature_counts = Counter(
            _text(feature.get("match_key")) for feature in signature.get("residues", [])
        )
        for key, occurrence_count in signature_counts.items():
            entry = profile_residues.setdefault(
                key,
                {
                    "feature_key": key,
                    "count": 0,
                    "support_count": 0,
                    "multiplicities": [],
                },
            )
            entry["count"] += occurrence_count
            entry["support_count"] += 1
            entry["multiplicities"].append(occurrence_count)
    for entry in profile_residues.values():
        entry["support"] = round(entry["support_count"] / n, 6)
        entry["required"] = entry["support"] >= required_support
        entry["multiplicity"] = _count_summary(entry.pop("multiplicities"))

    return {
        "schema_version": SIGNATURE_SCHEMA_VERSION,
        "profile_type": "coordination_network",
        "n_signatures": n,
        "required_support": required_support,
        "cofactor_distribution": [
            {"cofactor_key": key, "count": count, "support": round(count / n, 6)}
            for key, count in cofactor_counts.most_common()
        ],
        "dominant_cofactor_key": cofactor_counts.most_common(1)[0][0],
        "edges": profile_edges,
        "residues": dict(sorted(profile_residues.items())),
    }


def score_signature_against_profile(
    signature: Mapping[str, Any],
    profile: Mapping[str, Any],
    *,
    distance_tolerance_A: float = 0.5,
) -> Dict[str, Any]:
    """Score a query signature against a family profile."""
    profile_cofactor = _text(profile.get("dominant_cofactor_key"))
    query_cofactor = _text(signature.get("cofactor", {}).get("key"))
    cofactor_compatible = query_cofactor == profile_cofactor
    edge_results: Dict[str, Any] = {}
    coverage_scores: List[float] = []
    precision_scores: List[float] = []
    distance_scores: List[float] = []
    for category in CONTACT_CATEGORIES:
        entries = profile.get("edges", {}).get(category, {})
        query_edges = signature.get("edges", {}).get(category, [])
        query_by_key: Dict[str, List[Mapping[str, Any]]] = defaultdict(list)
        for edge in query_edges:
            query_by_key[_text(edge.get("match_key"))].append(edge)
        required_keys = {key for key, value in entries.items() if value.get("required")}
        observed_keys = set(query_by_key)
        expected_multiplicity = {
            key: max(1, int(round(_float(entries[key].get("multiplicity", {}).get("median")) or 1)))
            for key in entries
        }
        matched_required = {
            key for key in required_keys
            if len(query_by_key.get(key, [])) >= expected_multiplicity[key]
        }
        known_observed = observed_keys & set(entries)
        coverage = len(matched_required) / len(required_keys) if required_keys else 1.0
        precision = len(known_observed) / len(observed_keys) if observed_keys else 1.0
        expected_known_count = sum(expected_multiplicity[key] for key in entries)
        matched_known_count = sum(
            min(len(query_by_key.get(key, [])), expected_multiplicity[key])
            for key in entries
        )
        multiplicity = matched_known_count / expected_known_count if expected_known_count else 1.0
        for key in known_observed:
            expected = entries[key].get("distance", {}).get("median_A")
            if expected is not None:
                for edge in query_by_key[key]:
                    actual = _float(edge.get("distance_A"))
                    if actual is not None:
                        distance_scores.append(exp(-abs(actual - expected) / max(distance_tolerance_A, 1e-9)))
        coverage_scores.append(coverage)
        precision_scores.append(precision)
        edge_results[category] = {
            "coverage": round(coverage, 6),
            "precision": round(precision, 6),
            "multiplicity": round(multiplicity, 6),
            "required_count": len(required_keys),
            "matched_required_count": len(matched_required),
            "expected_edge_count": expected_known_count,
            "matched_edge_count": matched_known_count,
            "unexpected_features": sorted(observed_keys - set(entries)),
            "missing_required_features": sorted(required_keys - observed_keys),
        }

    multiplicity_scores = [edge_results[category]["multiplicity"] for category in CONTACT_CATEGORIES]
    edge_score = mean(coverage_scores + precision_scores + multiplicity_scores) if coverage_scores else 1.0
    profile_residues = profile.get("residues", {})
    query_residue_counts = Counter(
        _text(feature.get("match_key")) for feature in signature.get("residues", [])
    )
    required_residue_keys = {
        key for key, value in profile_residues.items() if value.get("required")
    }
    expected_residue_multiplicity = {
        key: max(1, int(round(_float(value.get("multiplicity", {}).get("median")) or 1)))
        for key, value in profile_residues.items()
    }
    matched_required_residues = {
        key for key in required_residue_keys
        if query_residue_counts.get(key, 0) >= expected_residue_multiplicity[key]
    }
    observed_residue_keys = set(query_residue_counts)
    known_residue_keys = observed_residue_keys & set(profile_residues)
    residue_coverage = (
        len(matched_required_residues) / len(required_residue_keys)
        if required_residue_keys else 1.0
    )
    residue_precision = (
        len(known_residue_keys) / len(observed_residue_keys)
        if observed_residue_keys else 1.0
    )
    expected_residue_count = sum(expected_residue_multiplicity.values())
    matched_residue_count = sum(
        min(query_residue_counts.get(key, 0), expected_residue_multiplicity[key])
        for key in profile_residues
    )
    residue_multiplicity = (
        matched_residue_count / expected_residue_count
        if expected_residue_count else 1.0
    )
    residue_score = mean((residue_coverage, residue_precision, residue_multiplicity))
    distance_score = mean(distance_scores) if distance_scores else None
    components = [(edge_score, 0.75), (residue_score, 0.15)]
    if distance_score is not None:
        components.append((distance_score, 0.10))
    denominator = sum(weight for _, weight in components)
    overall = sum(score * weight for score, weight in components) / denominator
    if not cofactor_compatible:
        overall *= 0.25
    return {
        "schema_version": SIGNATURE_SCHEMA_VERSION,
        "query_id": _text(signature.get("structure_id")),
        "cofactor_compatible": cofactor_compatible,
        "profile_cofactor_key": profile_cofactor,
        "scores": {
            "overall": round(float(overall), 6),
            "edge": round(float(edge_score), 6),
            "residue": round(float(residue_score), 6),
            "distance": None if distance_score is None else round(float(distance_score), 6),
        },
        "edges": edge_results,
        "residues": {
            "coverage": round(residue_coverage, 6),
            "precision": round(residue_precision, 6),
            "multiplicity": round(residue_multiplicity, 6),
            "required_count": len(required_residue_keys),
            "matched_required_count": len(matched_required_residues),
            "expected_residue_count": expected_residue_count,
            "matched_residue_count": matched_residue_count,
            "unexpected_features": sorted(observed_residue_keys - set(profile_residues)),
            "missing_required_features": sorted(required_residue_keys - observed_residue_keys),
        },
    }


def build_residue_conservation_map(
    analysis: Mapping[str, Any],
    signature: Mapping[str, Any],
    profile: Mapping[str, Any],
    *,
    residue_map: Optional[Mapping[Any, Any]] = None,
) -> Dict[Tuple[str, str, str, str, str, str], float]:
    """Map analyzed residue identities to profile support values.

    The returned keys match the viewer's model-aware residue identity:
    ``(model_id, residue, residue_number, chain, insertion_code, hetero_flag)``.
    This keeps conservation coloring separate from atom coordinates while
    allowing the cohesive viewer to color every atom belonging to a residue.
    """
    profile_residues = profile.get("residues", {})
    exact_support: Dict[Tuple[str, str, str], float] = {}
    fallback_support: Dict[Tuple[str, str], List[float]] = defaultdict(list)
    for feature in signature.get("residues", []):
        profile_entry = profile_residues.get(_text(feature.get("match_key")), {})
        support = _float(profile_entry.get("support"))
        if support is None:
            support = 0.0
        identity = (
            _text(feature.get("shell")).upper(),
            _text(feature.get("residue")).upper(),
            _text(feature.get("position")),
        )
        exact_support[identity] = support
        fallback_support[(identity[0], identity[1])].append(support)

    mapped: Dict[Tuple[str, str, str, str, str, str], float] = {}
    for row in _records(analysis.get("atoms")):
        shell = _text(row.get("shell")).upper()
        residue = _text(row.get("residue_name")).upper()
        position = _atom_position(row, residue_map)
        support = exact_support.get((shell, residue, position))
        if support is None:
            candidates = fallback_support.get((shell, residue), [])
            support = max(candidates) if candidates else 0.0
        viewer_key = (
            _text(row.get("model_id")),
            _text(row.get("residue_name")),
            _text(row.get("residue_number")),
            _text(row.get("chain")),
            _text(row.get("insertion_code")),
            _text(row.get("hetero_flag")),
        )
        mapped[viewer_key] = round(float(support), 6)
    return mapped


__all__ = [
    "build_network_signature",
    "compare_network_signatures",
    "build_network_profile",
    "score_signature_against_profile",
    "build_residue_conservation_map",
]
