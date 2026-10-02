"""Cofactor-to-target protein-wire analysis.

Protein wires are deliberately separate from coordination shells.  A shell is
an expanding geometric neighborhood; a wire is a chemically filtered graph of
relay atoms used to find paths between an explicit cofactor source and one or
more explicit targets.

The first implementation is intentionally conservative:

* cofactor atoms are source anchors;
* non-hydrogen N/O/S atoms, aromatic relay atoms, and optional waters are
  candidate relay nodes;
* atoms from the same residue are not connected by through-space wire edges;
* every edge keeps its atom identities, motif labels, distance, and interaction
  class; and
* paths are ranked with a transparent distance/hop cost rather than a learned
  score.

The graph supports generic, proton, electron, PCET, and opt-in residue-level
redox-relay scoring policies without changing the legacy single-structure
analyzer.
"""

from __future__ import annotations

from collections import defaultdict, deque
import heapq
import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd

from .io_utils import unpack_pdb_file
from .motif_registry import motif_for_atom


PathLike = Union[str, Path]

WIRE_NODE_COLUMNS = [
    "node_id", "node_role", "target_labels", "model_id", "residue_name",
    "residue_number", "chain", "insertion_code", "hetero_flag", "atom_name",
    "element", "motif", "redox_capable", "protonatable", "pcet_candidate",
    "relay_classes", "x", "y", "z",
]

WIRE_EDGE_COLUMNS = [
    "edge_id", "src_node_id", "dst_node_id", "wire_mode",
    "interaction_type", "geometry_status", "geometry_score",
    "proton_support", "electron_support", "pcet_support",
    "distance_A", "edge_cost",
]

WIRE_PATH_COLUMNS = [
    "path_id", "target_label", "source_node_id", "target_node_id", "hops",
    "wire_mode", "path_cost", "path_score", "node_ids", "edge_types",
    "edge_geometry", "residue_path",
]

WIRE_TARGET_COLUMNS = [
    "target_label", "requested_residue", "requested_residue_number",
    "requested_chain", "requested_atom", "matched", "status",
    "observed_residues", "matched_node_ids",
]

WATER_RESIDUES = frozenset({"HOH", "WAT", "H2O", "DOD", "SOL"})
WIRE_ELEMENTS = frozenset({"N", "O", "S"})
WIRE_MODES = frozenset({"generic", "relay", "proton", "electron", "pcet", "redox"})

REDOX_RESIDUES = frozenset({"TRP", "TYR", "HIS", "CYS", "MET"})
PROTONATABLE_RESIDUES = frozenset({
    "ASP", "GLU", "HIS", "LYS", "ARG", "TYR", "SER", "THR", "CYS",
})
REDOX_ELEMENTS = frozenset({"FE", "CU", "MN", "CO", "NI", "MO", "V", "W"})

# Carbon atoms are admitted only when they belong to a recognizable aromatic
# relay motif.  This lets phenyl/tyrosine/tryptophan relays participate without
# turning every aliphatic carbon in the protein into a graph node.
AROMATIC_ATOMS: Mapping[str, frozenset] = {
    "PHE": frozenset({"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
    "TYR": frozenset({"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
    "TRP": frozenset({"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}),
    "HIS": frozenset({"CG", "ND1", "CD2", "CE1", "NE2"}),
    "HID": frozenset({"CG", "ND1", "CD2", "CE1", "NE2"}),
    "HIE": frozenset({"CG", "ND1", "CD2", "CE1", "NE2"}),
    "HIP": frozenset({"CG", "ND1", "CD2", "CE1", "NE2"}),
}


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if value != value:  # NaN
            return ""
    except Exception:
        pass
    return str(value).strip()


def _normalized_names(names: Union[str, Sequence[str]]) -> Set[str]:
    if isinstance(names, str):
        names = names.split(",")
    return {_text(name).upper() for name in names if _text(name)}


def _normalized_wire_mode(mode: Any) -> str:
    normalized = _text(mode).lower() or "relay"
    if normalized == "generic":
        return "relay"
    if normalized not in WIRE_MODES:
        raise ValueError(f"wire_mode must be one of {sorted(WIRE_MODES)}")
    return normalized


def _residue_key(atom: Mapping[str, Any]) -> Tuple[str, str, str, str, str, str]:
    return (
        _text(atom.get("model_id")),
        _text(atom.get("residue", atom.get("residue_name"))).upper(),
        _text(atom.get("residue_number")),
        _text(atom.get("chain")),
        _text(atom.get("insertion_code")),
        _text(atom.get("hetero_flag")),
    )


def _node_id(atom: Mapping[str, Any]) -> str:
    residue = _text(atom.get("residue", atom.get("residue_name"))).upper()
    return "|".join(
        (
            f"m{_text(atom.get('model_id'))}",
            residue,
            _text(atom.get("residue_number")),
            _text(atom.get("chain")),
            _text(atom.get("insertion_code")),
            _text(atom.get("hetero_flag")),
            _text(atom.get("name", atom.get("atom_name"))).upper(),
        )
    )


def _coordinates(atom: Mapping[str, Any]) -> np.ndarray:
    coordinates = atom.get("coordinates")
    if coordinates is None:
        coordinates = [atom.get("x"), atom.get("y"), atom.get("z")]
    values = np.asarray(coordinates, dtype=float)
    if values.shape != (3,) or not np.isfinite(values).all():
        raise ValueError(f"Atom {_node_id(atom)} does not have finite 3D coordinates")
    return values


def _normalized_atom(atom: Mapping[str, Any]) -> Dict[str, Any]:
    residue = _text(atom.get("residue", atom.get("residue_name"))).upper()
    name = _text(atom.get("name", atom.get("atom_name"))).upper()
    element = _text(atom.get("element")).upper()
    normalized = dict(atom)
    normalized.update(
        {
            "residue": residue,
            "residue_number": atom.get("residue_number", ""),
            "name": name,
            "element": element,
            "motif": motif_for_atom(residue, name),
            "coordinates": _coordinates(atom),
        }
    )
    normalized["node_id"] = _node_id(normalized)
    return normalized


def _capabilities(atom: Mapping[str, Any], cofactor_names: Set[str]) -> Dict[str, Any]:
    residue = _text(atom.get("residue")).upper()
    element = _text(atom.get("element")).upper()
    redox_capable = residue in cofactor_names or residue in REDOX_RESIDUES or element in REDOX_ELEMENTS
    protonatable = (
        _is_water(atom)
        or residue in PROTONATABLE_RESIDUES
        or (residue in cofactor_names and element in WIRE_ELEMENTS)
    )
    aromatic_relay = _is_aromatic(atom)
    polar_relay = element in WIRE_ELEMENTS
    relay_classes = []
    if redox_capable:
        relay_classes.append("redox_capable")
    if protonatable:
        relay_classes.append("protonatable")
    if aromatic_relay:
        relay_classes.append("aromatic_relay")
    if polar_relay:
        relay_classes.append("polar_relay")
    if _is_water(atom):
        relay_classes.append("water")
    return {
        "redox_capable": bool(redox_capable),
        "protonatable": bool(protonatable),
        "pcet_candidate": bool(redox_capable and protonatable),
        "relay_classes": ";".join(relay_classes),
    }


def _is_water(atom: Mapping[str, Any]) -> bool:
    return (
        _text(atom.get("residue")).upper() in WATER_RESIDUES
        or _text(atom.get("motif")).lower() == "water"
    )


def _is_backbone(atom: Mapping[str, Any]) -> bool:
    return _text(atom.get("motif")).lower() == "backbone"


def _is_aromatic(atom: Mapping[str, Any]) -> bool:
    residue = _text(atom.get("residue")).upper()
    name = _text(atom.get("name")).upper()
    return name in AROMATIC_ATOMS.get(residue, frozenset())


def _is_relay_candidate(
    atom: Mapping[str, Any],
    *,
    cofactor_names: Set[str],
    include_water: bool,
    include_backbone: bool,
    wire_mode: str = "relay",
) -> bool:
    residue = _text(atom.get("residue")).upper()
    element = _text(atom.get("element")).upper()
    if residue in cofactor_names:
        return True
    if element in {"H", "D"}:
        return False
    if _is_water(atom) and not include_water:
        return False
    if _is_backbone(atom) and not include_backbone:
        return False
    if wire_mode == "redox":
        # Redox-relay mode is intentionally residue-level in its chemistry:
        # keep only cofactor atoms and redox-capable side-chain atoms.  This
        # prevents an arbitrary Asp/Glu/water chain from being mistaken for a
        # sequence of redox hops.  Polar side-chain atoms are retained when
        # they belong to a redox-active residue (for example Tyr-OH or Cys-S).
        return residue in REDOX_RESIDUES and element in WIRE_ELEMENTS | {"S"}
    return element in WIRE_ELEMENTS or _is_aromatic(atom)


def _target_matches(atom: Mapping[str, Any], selector: Mapping[str, Any]) -> bool:
    residue = _text(selector.get("residue", selector.get("resname"))).upper()
    atom_name = _text(selector.get("atom", selector.get("atom_name"))).upper()
    chain = _text(selector.get("chain"))
    number = _text(selector.get("residue_number", selector.get("resnum")))
    insertion = _text(selector.get("insertion_code"))
    model_id = _text(selector.get("model_id"))

    if residue and _text(atom.get("residue")).upper() != residue:
        return False
    if atom_name and _text(atom.get("name")).upper() != atom_name:
        return False
    if chain and _text(atom.get("chain")) != chain:
        return False
    if number and _text(atom.get("residue_number")) != number:
        return False
    if insertion and _text(atom.get("insertion_code")) != insertion:
        return False
    if model_id and _text(atom.get("model_id")) != model_id:
        return False
    motif = _text(selector.get("motif")).lower()
    if motif and _text(atom.get("motif")).lower() != motif:
        return False
    return True


def _target_label(selector: Mapping[str, Any], index: int) -> str:
    explicit = _text(selector.get("label"))
    if explicit:
        return explicit
    residue = _text(selector.get("residue", selector.get("resname"))).upper() or "target"
    number = _text(selector.get("residue_number", selector.get("resnum")))
    chain = _text(selector.get("chain"))
    atom = _text(selector.get("atom", selector.get("atom_name"))).upper()
    suffix = f" {number}" if number else ""
    if chain:
        suffix += f" {chain}"
    if atom:
        suffix += f":{atom}"
    return f"{residue}{suffix}" if suffix else f"{residue}_{index}"


def _interaction_type(first: Mapping[str, Any], second: Mapping[str, Any]) -> str:
    if _is_water(first) or _is_water(second):
        return "water_mediated_candidate"
    first_element = _text(first.get("element")).upper()
    second_element = _text(second.get("element")).upper()
    if "S" in {first_element, second_element}:
        return "sulfur_contact"
    if _is_aromatic(first) or _is_aromatic(second):
        return "aromatic_redox_contact"
    if first_element in WIRE_ELEMENTS and second_element in WIRE_ELEMENTS:
        return "polar_contact"
    return "through_space_contact"


def _angle_degrees(first: np.ndarray, vertex: np.ndarray, second: np.ndarray) -> float:
    first_vector = first - vertex
    second_vector = second - vertex
    denominator = float(np.linalg.norm(first_vector) * np.linalg.norm(second_vector))
    if denominator == 0.0:
        return 0.0
    cosine = float(np.dot(first_vector, second_vector) / denominator)
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def _hydrogen_neighbors(
    donor: Mapping[str, Any],
    atoms_by_residue: Mapping[Tuple[str, str, str, str, str, str], Sequence[Mapping[str, Any]]],
) -> List[Mapping[str, Any]]:
    neighbors = []
    donor_coordinates = donor["coordinates"]
    for atom in atoms_by_residue.get(_residue_key(donor), ()):
        element = _text(atom.get("element")).upper()
        name = _text(atom.get("name")).upper()
        if element not in {"H", "D"} and not name.startswith(("H", "D")):
            continue
        if float(np.linalg.norm(donor_coordinates - atom["coordinates"])) <= 1.35:
            neighbors.append(atom)
    return neighbors


def _ring_geometry(
    first: Mapping[str, Any],
    second: Mapping[str, Any],
    aromatic_groups: Mapping[Tuple[str, str, str, str, str, str], Tuple[np.ndarray, np.ndarray]],
    *,
    center_distance_cutoff: float = 5.5,
) -> Tuple[str, float]:
    first_group = aromatic_groups.get(_residue_key(first))
    second_group = aromatic_groups.get(_residue_key(second))
    if first_group is None or second_group is None:
        return "aromatic_unvalidated", 0.25
    first_center, first_normal = first_group
    second_center, second_normal = second_group
    center_distance = float(np.linalg.norm(first_center - second_center))
    alignment = abs(float(np.dot(first_normal, second_normal)))
    if center_distance <= center_distance_cutoff and alignment >= 0.70:
        return "aromatic_face_to_face", 1.0
    if center_distance <= center_distance_cutoff and alignment <= 0.50:
        return "aromatic_edge_to_face", 0.75
    return "aromatic_unvalidated", 0.25


def _edge_geometry(
    first: Mapping[str, Any],
    second: Mapping[str, Any],
    *,
    all_atoms_by_residue: Mapping[Tuple[str, str, str, str, str, str], Sequence[Mapping[str, Any]]],
    aromatic_groups: Mapping[Tuple[str, str, str, str, str, str], Tuple[np.ndarray, np.ndarray]],
    distance: float,
    aromatic_center_distance_cutoff: float = 5.5,
) -> Tuple[str, float]:
    first_element = _text(first.get("element")).upper()
    second_element = _text(second.get("element")).upper()
    polar_pair = first_element in WIRE_ELEMENTS and second_element in WIRE_ELEMENTS

    if _is_aromatic(first) and _is_aromatic(second):
        return _ring_geometry(
            first,
            second,
            aromatic_groups,
            center_distance_cutoff=aromatic_center_distance_cutoff,
        )

    if polar_pair:
        # If deposited hydrogens are present, validate the donor-heavy-
        # atom/H/acceptor angle.  Without hydrogens, preserve the edge as an
        # explicitly inferred polar contact rather than pretending geometry
        # was observed.
        for donor, acceptor in ((first, second), (second, first)):
            for hydrogen in _hydrogen_neighbors(donor, all_atoms_by_residue):
                hydrogen_acceptor_distance = float(
                    np.linalg.norm(hydrogen["coordinates"] - acceptor["coordinates"])
                )
                angle = _angle_degrees(donor["coordinates"], hydrogen["coordinates"], acceptor["coordinates"])
                if hydrogen_acceptor_distance <= 2.70 and angle >= 120.0:
                    return "validated_hydrogen_bond", 1.0
        return "inferred_polar_contact", 0.55

    if _is_water(first) or _is_water(second):
        return "water_contact", 0.50
    if _is_aromatic(first) or _is_aromatic(second):
        return "aromatic_polar_contact", 0.45
    return "distance_only", 0.20


def _edge_support(
    first: Mapping[str, Any],
    second: Mapping[str, Any],
    interaction_type: str,
    geometry_status: str,
    cofactor_names: Set[str],
) -> Tuple[float, float, float]:
    """Return proton, electron, and PCET support for one graph edge.

    These are structural plausibility scores, not transfer probabilities.  A
    PCET score is intentionally conservative: both the proton and electron
    interpretations must have support.
    """
    first_capabilities = _capabilities(first, cofactor_names)
    second_capabilities = _capabilities(second, cofactor_names)
    redox_pair = first_capabilities["redox_capable"] or second_capabilities["redox_capable"]
    proton_pair = first_capabilities["protonatable"] or second_capabilities["protonatable"]

    if geometry_status == "validated_hydrogen_bond":
        proton_support = 1.0
    elif geometry_status == "inferred_polar_contact":
        proton_support = 0.65 if proton_pair else 0.25
    elif geometry_status == "water_contact":
        proton_support = 0.70
    elif interaction_type == "cofactor_contact" and proton_pair:
        proton_support = 0.45
    else:
        proton_support = 0.10

    if geometry_status in {"aromatic_face_to_face", "aromatic_edge_to_face"}:
        electron_support = 0.95 if redox_pair else 0.45
    elif interaction_type == "sulfur_contact":
        electron_support = 0.80 if redox_pair else 0.45
    elif interaction_type == "cofactor_contact":
        electron_support = 0.75 if redox_pair else 0.35
    elif geometry_status == "aromatic_polar_contact":
        electron_support = 0.55 if redox_pair else 0.25
    elif interaction_type == "polar_contact" and redox_pair:
        electron_support = 0.35
    else:
        electron_support = 0.10

    pcet_support = min(proton_support, electron_support)
    return proton_support, electron_support, pcet_support


def _edge_cost(
    distance: float,
    interaction_type: str,
    geometry_status: str,
    geometry_score: float,
    max_hop_distance: float,
    wire_mode: str,
    proton_support: float,
    electron_support: float,
    pcet_support: float,
) -> float:
    # Every hop has a cost, so paths do not grow arbitrarily long.  Distances
    # near the cutoff are penalized, while chemically typed contacts remain
    # inspectable rather than being hidden inside an opaque score.
    if wire_mode == "proton":
        interaction_penalty = {
            "polar_contact": 0.0,
            "sulfur_contact": 0.15,
            "aromatic_redox_contact": 0.75,
            "water_mediated_candidate": 0.05,
            "through_space_contact": 0.80,
            "cofactor_contact": 0.0,
        }.get(interaction_type, 0.80)
        geometry_penalty = 0.0 if geometry_status == "validated_hydrogen_bond" else 0.35
    elif wire_mode in {"electron", "redox"}:
        interaction_penalty = {
            "polar_contact": 0.35,
            "sulfur_contact": 0.05,
            "aromatic_redox_contact": 0.0,
            "water_mediated_candidate": 0.90,
            "through_space_contact": 0.70,
            "cofactor_contact": 0.0,
        }.get(interaction_type, 0.70)
        geometry_penalty = 0.0 if geometry_status in {
            "aromatic_face_to_face", "aromatic_edge_to_face",
        } else 0.20
    elif wire_mode == "pcet":
        interaction_penalty = 0.25 if interaction_type == "cofactor_contact" else 0.0
        geometry_penalty = 0.0 if pcet_support >= 0.60 else 0.75 * (1.0 - pcet_support)
    else:
        interaction_penalty = {
            "polar_contact": 0.0,
            "sulfur_contact": 0.10,
            "aromatic_redox_contact": 0.15,
            "water_mediated_candidate": 0.20,
            "through_space_contact": 0.60,
            "cofactor_contact": 0.0,
        }.get(interaction_type, 0.60)
        geometry_penalty = 0.15 * (1.0 - geometry_score)
    if wire_mode == "proton":
        interpretation_penalty = 1.0 - proton_support
    elif wire_mode in {"electron", "redox"}:
        interpretation_penalty = 1.0 - electron_support
    elif wire_mode == "pcet":
        interpretation_penalty = 1.0 - pcet_support
    else:
        interpretation_penalty = 0.0
    return 1.0 + distance / max_hop_distance + interaction_penalty + geometry_penalty + interpretation_penalty


def _empty_tables() -> Dict[str, pd.DataFrame]:
    return {
        "nodes": pd.DataFrame(columns=WIRE_NODE_COLUMNS),
        "edges": pd.DataFrame(columns=WIRE_EDGE_COLUMNS),
        "paths": pd.DataFrame(columns=WIRE_PATH_COLUMNS),
        "target_status": pd.DataFrame(columns=WIRE_TARGET_COLUMNS),
    }


def _neighbor_cells(coordinates: np.ndarray, cell_size: float) -> Iterable[Tuple[int, int, int]]:
    cell = np.floor(coordinates / cell_size).astype(int)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                yield int(cell[0] + dx), int(cell[1] + dy), int(cell[2] + dz)


def _reachable(adjacency: Mapping[str, Sequence[Tuple[str, Dict[str, Any]]]], starts: Set[str]) -> Set[str]:
    seen = set(starts)
    queue = deque(starts)
    while queue:
        node = queue.popleft()
        for neighbor, _ in adjacency.get(node, ()):
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    return seen


def _shortest_paths(
    adjacency: Mapping[str, Sequence[Tuple[str, Dict[str, Any]]]],
    source_ids: Set[str],
    target_ids: Set[str],
    target_labels_by_id: Mapping[str, Sequence[str]],
    *,
    max_hops: int,
    wire_mode: str,
    residue_key_by_node: Optional[Mapping[str, Tuple[str, str, str, str, str, str]]] = None,
    enforce_unique_residues: bool = False,
) -> List[Dict[str, Any]]:
    """Find one transparent best path for each reachable target atom."""
    heap: List[
        Tuple[float, int, Tuple[str, ...], str, Tuple[str, ...], Tuple[str, ...], Tuple[str, ...]]
    ] = []
    for source_id in sorted(source_ids):
        heapq.heappush(heap, (0.0, 0, (source_id,), source_id, (source_id,), (), ()))

    best_by_state: Dict[Any, Tuple[float, int, Tuple[str, ...]]] = {}
    paths: List[Dict[str, Any]] = []
    reached_targets: Set[str] = set()
    while heap:
        cost, hops, path_key, node, node_path, edge_types, edge_geometry = heapq.heappop(heap)
        state = (cost, hops, path_key)
        residue_path = tuple(
            residue_key_by_node[node_id] for node_id in node_path
        ) if enforce_unique_residues and residue_key_by_node is not None else ()
        state_key: Any = (node, residue_path) if enforce_unique_residues else node
        previous = best_by_state.get(state_key)
        if previous is not None and previous <= state:
            continue
        best_by_state[state_key] = state

        if node in target_ids and node not in source_ids and node not in reached_targets:
            reached_targets.add(node)
            paths.append(
                {
                    "target_label": "; ".join(target_labels_by_id.get(node, ("target",))),
                    "source_node_id": node_path[0],
                    "target_node_id": node,
                    "hops": hops,
                    "wire_mode": wire_mode,
                    "path_cost": cost,
                    "path_score": 1.0 / (1.0 + cost),
                    "node_ids": json.dumps(list(node_path)),
                    "edge_types": json.dumps(list(edge_types)),
                    "edge_geometry": json.dumps(list(edge_geometry)),
                }
            )

        if hops >= max_hops:
            continue
        for neighbor, edge in adjacency.get(node, ()):
            if neighbor in node_path:
                continue
            if enforce_unique_residues and residue_key_by_node is not None:
                neighbor_residue = residue_key_by_node[neighbor]
                if neighbor_residue in residue_path:
                    continue
            next_hops = hops + 1
            next_cost = cost + float(edge["edge_cost"])
            next_path = node_path + (neighbor,)
            next_types = edge_types + (str(edge["interaction_type"]),)
            next_geometry = edge_geometry + (str(edge["geometry_status"]),)
            heapq.heappush(
                heap,
                (next_cost, next_hops, next_path, neighbor, next_path, next_types, next_geometry),
            )
    return paths


def build_protein_wire_network(
    atoms: Sequence[Mapping[str, Any]],
    cofactor_resname: Union[str, Sequence[str]],
    targets: Sequence[Mapping[str, Any]],
    *,
    max_hop_distance: float = 3.6,
    max_hops: int = 8,
    include_water: bool = True,
    include_backbone: bool = False,
    first_model_only: bool = False,
    wire_mode: str = "generic",
    redox_hop_distance: float = 6.0,
) -> Dict[str, pd.DataFrame]:
    """Build a cofactor-to-target protein-wire graph.

    ``targets`` contains selectors such as ``{"residue": "TRP",
    "residue_number": 91, "chain": "A", "atom": "NE1"}``.  Omitting
    ``atom`` selects the eligible relay atoms in that residue.  A selector may
    include ``label`` for a stable human-readable target name.

    ``wire_mode`` is ``generic``, ``proton``, ``electron``, ``pcet``, or
    ``redox``.  The proton policy favors polar and water-mediated hops; the
    electron policy favors aromatic/redox and sulfur hops.  ``redox`` is an
    opt-in residue-level relay mode: only cofactor atoms and redox-capable
    residues are graph nodes, aromatic contacts may extend to
    ``redox_hop_distance`` (6.0 Å by default), and a path may not revisit a
    residue.  Neither policy claims a quantum or thermodynamic transfer rate:
    each is an explainable ranking of structural candidates.

    The returned ``nodes`` and ``edges`` tables describe the source-to-target
    connected subgraph.  ``paths`` contains the lowest-cost path to each
    reachable target atom, retaining the ordered node IDs, interaction
    classes, and geometry evidence as JSON strings for CSV/JSON export.
    ``target_status`` reports every requested target, including residue-name
    or atom-number mismatches that would otherwise be silently omitted.
    """
    if max_hop_distance <= 0:
        raise ValueError("max_hop_distance must be positive")
    if redox_hop_distance <= 0:
        raise ValueError("redox_hop_distance must be positive")
    if max_hops < 1:
        raise ValueError("max_hops must be at least 1")
    wire_mode = _normalized_wire_mode(wire_mode)
    if not targets:
        raise ValueError("at least one target selector is required")

    cofactor_names = _normalized_names(cofactor_resname)
    if not cofactor_names:
        raise ValueError("at least one cofactor residue name is required")

    normalized_atoms = [_normalized_atom(atom) for atom in atoms]
    if first_model_only and normalized_atoms:
        first_model = _text(normalized_atoms[0].get("model_id"))
        normalized_atoms = [atom for atom in normalized_atoms if _text(atom.get("model_id")) == first_model]

    graph_distance = redox_hop_distance if wire_mode == "redox" else max_hop_distance
    candidates = [
        atom for atom in normalized_atoms
        if _is_relay_candidate(
            atom,
            cofactor_names=cofactor_names,
            include_water=include_water,
            include_backbone=include_backbone,
            wire_mode=wire_mode,
        )
    ]
    by_id = {atom["node_id"]: atom for atom in candidates}
    source_ids = {
        atom["node_id"]
        for atom in candidates
        if _text(atom.get("residue")).upper() in cofactor_names
    }
    if not source_ids:
        raise ValueError(f"no cofactor atoms found for {sorted(cofactor_names)}")

    target_labels_by_id: Dict[str, List[str]] = defaultdict(list)
    target_status_rows: List[Dict[str, Any]] = []
    for index, selector in enumerate(targets, start=1):
        label = _target_label(selector, index)
        selector_atom = _text(selector.get("atom", selector.get("atom_name")))
        full_matches = [atom for atom in normalized_atoms if _target_matches(atom, selector)]
        matched_node_ids: List[str] = []
        for atom in full_matches:
            if atom["node_id"] not in by_id:
                # Explicit atom selectors are allowed to promote a chemically
                # unusual target atom into the endpoint set, while residue-only
                # selectors remain limited to relay-eligible atoms.
                if not selector_atom:
                    continue
                if _text(atom.get("element")).upper() in {"H", "D"}:
                    continue
                by_id[atom["node_id"]] = atom
            if label not in target_labels_by_id[atom["node_id"]]:
                target_labels_by_id[atom["node_id"]].append(label)
            matched_node_ids.append(atom["node_id"])

        observed_residues: Set[str] = set()
        if not full_matches:
            requested_number = _text(selector.get("residue_number", selector.get("resnum")))
            requested_chain = _text(selector.get("chain"))
            requested_insertion = _text(selector.get("insertion_code"))
            requested_model = _text(selector.get("model_id"))
            for atom in normalized_atoms:
                if requested_number and _text(atom.get("residue_number")) != requested_number:
                    continue
                if requested_chain and _text(atom.get("chain")) != requested_chain:
                    continue
                if requested_insertion and _text(atom.get("insertion_code")) != requested_insertion:
                    continue
                if requested_model and _text(atom.get("model_id")) != requested_model:
                    continue
                observed_residues.add(_text(atom.get("residue")).upper())

        if matched_node_ids:
            target_status = "matched"
        elif full_matches:
            target_status = "not_relay_eligible"
        elif observed_residues:
            requested_residue = _text(selector.get("residue", selector.get("resname"))).upper()
            target_status = "residue_name_mismatch" if requested_residue not in observed_residues else "atom_not_found"
        else:
            target_status = "target_not_found"
        target_status_rows.append(
            {
                "target_label": label,
                "requested_residue": _text(selector.get("residue", selector.get("resname"))).upper(),
                "requested_residue_number": _text(selector.get("residue_number", selector.get("resnum"))),
                "requested_chain": _text(selector.get("chain")),
                "requested_atom": selector_atom,
                "matched": bool(matched_node_ids),
                "status": target_status,
                "observed_residues": "; ".join(sorted(observed_residues)),
                "matched_node_ids": json.dumps(matched_node_ids),
            }
        )

    target_ids = set(target_labels_by_id)
    if not target_ids:
        unresolved = ", ".join(row["target_label"] for row in target_status_rows)
        raise ValueError(f"none of the target selectors matched eligible atoms: {unresolved}")

    atoms_by_residue: Dict[Tuple[str, str, str, str, str, str], List[Dict[str, Any]]] = defaultdict(list)
    for atom in normalized_atoms:
        atoms_by_residue[_residue_key(atom)].append(atom)
    aromatic_groups: Dict[Tuple[str, str, str, str, str, str], Tuple[np.ndarray, np.ndarray]] = {}
    for residue_key, residue_atoms in atoms_by_residue.items():
        aromatic_atoms = [atom for atom in residue_atoms if _is_aromatic(atom)]
        if len(aromatic_atoms) < 3:
            continue
        coordinates = np.asarray([atom["coordinates"] for atom in aromatic_atoms], dtype=float)
        center = coordinates.mean(axis=0)
        _, _, vh = np.linalg.svd(coordinates - center, full_matrices=False)
        normal = vh[-1]
        normal_norm = float(np.linalg.norm(normal))
        if normal_norm:
            aromatic_groups[residue_key] = (center, normal / normal_norm)

    # Spatial hashing keeps the graph construction practical for full protein
    # structures without adding scipy as a runtime dependency.
    cells: Dict[Tuple[str, int, int, int], List[str]] = defaultdict(list)
    for node_id, atom in by_id.items():
        coords = atom["coordinates"]
        cell = tuple(np.floor(coords / graph_distance).astype(int))
        cells[(_text(atom.get("model_id")), *cell)].append(node_id)

    adjacency: Dict[str, List[Tuple[str, Dict[str, Any]]]] = defaultdict(list)
    edge_records: List[Dict[str, Any]] = []
    edge_seen: Set[Tuple[str, str]] = set()
    for node_id, atom in by_id.items():
        for cell in _neighbor_cells(atom["coordinates"], graph_distance):
            cell_key = (_text(atom.get("model_id")), *cell)
            for other_id in cells.get(cell_key, ()):
                if other_id == node_id:
                    continue
                other = by_id[other_id]
                if _residue_key(atom) == _residue_key(other):
                    continue
                if _text(atom.get("model_id")) != _text(other.get("model_id")):
                    continue
                pair = tuple(sorted((node_id, other_id)))
                if pair in edge_seen:
                    continue
                distance = float(np.linalg.norm(atom["coordinates"] - other["coordinates"]))
                if distance > graph_distance:
                    continue
                edge_seen.add(pair)
                interaction_type = _interaction_type(atom, other)
                if node_id in source_ids or other_id in source_ids:
                    interaction_type = "cofactor_contact"
                geometry_status, geometry_score = _edge_geometry(
                    atom,
                    other,
                    all_atoms_by_residue=atoms_by_residue,
                    aromatic_groups=aromatic_groups,
                    distance=distance,
                    aromatic_center_distance_cutoff=6.5 if wire_mode == "redox" else 5.5,
                )
                proton_support, electron_support, pcet_support = _edge_support(
                    atom,
                    other,
                    interaction_type,
                    geometry_status,
                    cofactor_names,
                )
                edge = {
                    "edge_id": f"wire_edge_{len(edge_records) + 1}",
                    "src_node_id": pair[0],
                    "dst_node_id": pair[1],
                    "wire_mode": wire_mode,
                    "interaction_type": interaction_type,
                    "geometry_status": geometry_status,
                    "geometry_score": geometry_score,
                    "proton_support": proton_support,
                    "electron_support": electron_support,
                    "pcet_support": pcet_support,
                    "distance_A": distance,
                    "edge_cost": _edge_cost(
                        distance,
                        interaction_type,
                        geometry_status,
                        geometry_score,
                        graph_distance,
                        wire_mode,
                        proton_support,
                        electron_support,
                        pcet_support,
                    ),
                }
                edge_records.append(edge)
                adjacency[pair[0]].append((pair[1], edge))
                adjacency[pair[1]].append((pair[0], edge))

    source_reachable = _reachable(adjacency, source_ids)
    target_reachable = _reachable(adjacency, target_ids)
    network_node_ids = source_reachable.intersection(target_reachable)
    network_node_ids.update(source_ids.intersection(source_reachable))
    network_node_ids.update(target_ids.intersection(target_reachable))

    paths = _shortest_paths(
        adjacency,
        source_ids,
        target_ids,
        target_labels_by_id,
        max_hops=max_hops,
        wire_mode=wire_mode,
        residue_key_by_node={node_id: _residue_key(atom) for node_id, atom in by_id.items()},
        enforce_unique_residues=wire_mode == "redox",
    )
    for index, path in enumerate(sorted(paths, key=lambda row: (float(row["path_cost"]), row["target_label"])), start=1):
        path["path_id"] = f"wire_path_{index}"
        node_ids = json.loads(path["node_ids"])
        path["residue_path"] = json.dumps([
            {
                "residue": by_id[node_id].get("residue", ""),
                "residue_number": by_id[node_id].get("residue_number", ""),
                "chain": by_id[node_id].get("chain", ""),
                "atom": by_id[node_id].get("name", ""),
                "motif": by_id[node_id].get("motif", ""),
            }
            for node_id in node_ids
        ])

    node_rows: List[Dict[str, Any]] = []
    for node_id in sorted(network_node_ids):
        atom = by_id[node_id]
        is_source = node_id in source_ids
        is_target = node_id in target_ids
        capabilities = _capabilities(atom, cofactor_names)
        roles = []
        if is_source:
            roles.append("source")
        if is_target:
            roles.append("target")
        if not roles:
            roles.append("relay")
        coordinates = atom["coordinates"]
        node_rows.append(
            {
                "node_id": node_id,
                "node_role": "+".join(roles),
                "target_labels": "; ".join(target_labels_by_id.get(node_id, ())),
                "model_id": atom.get("model_id", ""),
                "residue_name": atom.get("residue", ""),
                "residue_number": atom.get("residue_number", ""),
                "chain": atom.get("chain", ""),
                "insertion_code": atom.get("insertion_code", ""),
                "hetero_flag": atom.get("hetero_flag", ""),
                "atom_name": atom.get("name", ""),
                "element": atom.get("element", ""),
                "motif": atom.get("motif", ""),
                **capabilities,
                "x": float(coordinates[0]),
                "y": float(coordinates[1]),
                "z": float(coordinates[2]),
            }
        )

    network_edges = [
        edge for edge in edge_records
        if edge["src_node_id"] in network_node_ids and edge["dst_node_id"] in network_node_ids
    ]
    paths.sort(key=lambda row: (float(row["path_cost"]), row["target_label"]))
    return {
        "nodes": pd.DataFrame(node_rows, columns=WIRE_NODE_COLUMNS),
        "edges": pd.DataFrame(network_edges, columns=WIRE_EDGE_COLUMNS),
        "paths": pd.DataFrame(paths, columns=WIRE_PATH_COLUMNS),
        "target_status": pd.DataFrame(target_status_rows, columns=WIRE_TARGET_COLUMNS),
    }


def analyze_protein_wires(
    structure_path: PathLike,
    cofactor_resname: Union[str, Sequence[str]],
    targets: Sequence[Mapping[str, Any]],
    **kwargs: Any,
) -> Dict[str, pd.DataFrame]:
    """Load a PDB/mmCIF structure and build its cofactor-to-target wires."""
    path = Path(structure_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Structure file not found: {path}")
    _, atoms = unpack_pdb_file(str(path))
    return build_protein_wire_network(atoms, cofactor_resname, targets, **kwargs)


__all__ = [
    "WIRE_EDGE_COLUMNS",
    "WIRE_MODES",
    "WIRE_NODE_COLUMNS",
    "WIRE_PATH_COLUMNS",
    "WIRE_TARGET_COLUMNS",
    "analyze_protein_wires",
    "build_protein_wire_network",
]
