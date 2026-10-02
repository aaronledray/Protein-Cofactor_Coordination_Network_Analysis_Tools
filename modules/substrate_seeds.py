"""Hypothetical substrate-point seeds and chains through a coordination network.

A substrate point is a user-supplied position (for example a modelled O2, H2O2,
or ligand-binding point) that is not part of the deposited structure. It is
treated as an extra seed: the network atoms it touches are reported, and
ranked chains are traced from the point through the existing coordination
network (adjacent-shell contacts) down to cofactor atoms.

Everything here is opt-in, side-effect free, and built on the tables returned
by ``analyze_structure()``. Chains are structural hypotheses about contact
connectivity from the supplied point, not binding, reactivity, or energetics.
"""

import heapq
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from .coordination_api import analyze_structure

DEFAULT_CONTACT_CUTOFF_A = 4.0
POINT_NODE = "substrate_point"

CONTACT_COLUMNS = [
    "structure_id", "site_id", "residue_name", "residue_number", "chain",
    "insertion_code", "model_id", "atom_name", "element", "motif", "shell",
    "distance_A",
]
SUMMARY_COLUMNS = [
    "structure_id", "site_id", "chain_id", "rank", "entry_residue", "entry_atom",
    "entry_shell", "terminal_atom", "n_hops", "total_length_A", "shells_traversed",
]
STEP_COLUMNS = [
    "structure_id", "site_id", "chain_id", "step", "node_type", "residue_name",
    "residue_number", "chain", "insertion_code", "model_id", "atom_name", "motif",
    "shell", "edge_from_previous", "edge_distance_A",
]

Point = Union[Sequence[float], Mapping[str, Any]]


def _node(row: Mapping[str, Any]) -> Tuple:
    return (
        str(row["residue_name"]), str(row["residue_number"]), str(row["chain"]),
        str(row.get("insertion_code", "") or "").strip(), str(row["model_id"]),
        str(row["atom_name"]),
    )


def _contact_node(row: Mapping[str, Any], prefix: str) -> Tuple:
    return (
        str(row[f"{prefix}_resname"]), str(row[f"{prefix}_resnum"]),
        str(row[f"{prefix}_chain"]),
        str(row.get(f"{prefix}_insertion_code", "") or "").strip(),
        str(row[f"{prefix}_model_id"]), str(row[f"{prefix}_atom"]),
    )


def resolve_substrate_point(site_atoms: pd.DataFrame, point: Point) -> np.ndarray:
    """Return xyz for an explicit point or a cofactor-atom-relative one.

    ``point`` is either ``(x, y, z)`` or ``{"cofactor_atom": "FE",
    "offset": (dx, dy, dz)}``; the latter is resolved within one site.
    """
    if isinstance(point, Mapping):
        name = point.get("cofactor_atom")
        offset = point.get("offset")
        if not name or offset is None or len(offset) != 3:
            raise ValueError("relative point needs 'cofactor_atom' and a 3-vector 'offset'")
        anchors = site_atoms[
            (site_atoms["shell"] == "Cofactor") & (site_atoms["atom_name"] == name)
        ]
        if anchors.empty:
            raise ValueError(f"cofactor atom {name!r} not found in this site")
        base = anchors.iloc[0][["x", "y", "z"]].to_numpy(dtype=float)
        return base + np.asarray(offset, dtype=float)
    values = np.asarray(point, dtype=float)
    if values.shape != (3,) or not np.isfinite(values).all():
        raise ValueError("point must be three finite coordinates")
    return values


def _adjacency(contacts: pd.DataFrame) -> Dict[Tuple, List[Tuple[Tuple, float, str]]]:
    graph: Dict[Tuple, List[Tuple[Tuple, float, str]]] = {}
    for row in contacts.to_dict("records"):
        src, dst = _contact_node(row, "src"), _contact_node(row, "dst")
        kind = "direct_coordination" if row.get("direct_coordination") else "shell_contact"
        distance = float(row["distance_A"])
        graph.setdefault(src, []).append((dst, distance, kind))
        graph.setdefault(dst, []).append((src, distance, kind))
    return graph


def _best_chain(graph, start: Tuple, cofactor_nodes: set, first_distance: float):
    """Fewest-hops (then shortest) path from ``start`` to any cofactor atom."""
    queue = [(1, first_distance, start, None)]
    best: Dict[Tuple, Tuple[int, float]] = {}
    parents: Dict[Tuple, Tuple[Optional[Tuple], str, float]] = {}
    while queue:
        hops, length, node, parent = heapq.heappop(queue)
        if node in best:
            continue
        best[node] = (hops, length)
        if parent is not None:
            parents[node] = parent
        if node in cofactor_nodes:
            path = [node]
            while path[-1] in parents and parents[path[-1]][0] is not None:
                path.append(parents[path[-1]][0])
            path.reverse()
            return path, parents, hops, length
        for neighbor, distance, kind in sorted(graph.get(node, ()), key=lambda e: (e[1], e[0])):
            if neighbor not in best:
                heapq.heappush(queue, (hops + 1, length + distance, neighbor, (node, kind, distance)))
    return None


def substrate_chains(
    tables: Mapping[str, pd.DataFrame],
    point: Point,
    *,
    contact_cutoff: float = DEFAULT_CONTACT_CUTOFF_A,
    max_chains: int = 10,
) -> Dict[str, pd.DataFrame]:
    """Trace ranked chains from a hypothetical substrate point to the cofactor.

    Per site, network atoms within ``contact_cutoff`` of the point are the
    point's contacts. For each contacted residue, the best chain (fewest hops,
    then shortest summed distance) runs from its nearest contacted atom along
    shell-adjacent contacts to a cofactor atom; chains are ranked globally per
    site and capped at ``max_chains``. Returns ``seed``, ``contacts``,
    ``chains`` (one row per chain) and ``steps`` (one row per node).
    """
    if contact_cutoff <= 0:
        raise ValueError("contact_cutoff must be positive")
    if max_chains < 1:
        raise ValueError("max_chains must be at least 1")

    atoms, all_contacts = tables["atoms"], tables["contacts"]
    seeds, contact_rows, summaries, steps = [], [], [], []
    for site_id, site_atoms in atoms.groupby("site_id", sort=True):
        xyz = resolve_substrate_point(site_atoms, point)
        structure_id = str(site_atoms["structure_id"].iloc[0])
        seeds.append({"structure_id": structure_id, "site_id": site_id,
                      "x": xyz[0], "y": xyz[1], "z": xyz[2], "contact_cutoff_A": contact_cutoff})

        coordinates = site_atoms[["x", "y", "z"]].to_numpy(dtype=float)
        distances = np.linalg.norm(coordinates - xyz, axis=1)
        near = site_atoms[distances <= contact_cutoff].assign(distance_A=distances[distances <= contact_cutoff])
        near = near.sort_values(["distance_A", "residue_number", "atom_name"], kind="mergesort")
        for row in near.to_dict("records"):
            contact_rows.append({column: row.get(column) for column in CONTACT_COLUMNS})
        if near.empty:
            continue

        site_contacts = all_contacts[all_contacts["site_id"] == site_id]
        graph = _adjacency(site_contacts)
        cofactor_nodes = {_node(r) for r in site_atoms[site_atoms["shell"] == "Cofactor"].to_dict("records")}
        info = {_node(r): r for r in site_atoms.to_dict("records")}

        seen_residues, candidates = set(), []
        for row in near.to_dict("records"):
            residue = (row["residue_name"], str(row["residue_number"]), str(row["chain"]),
                       str(row["model_id"]))
            if residue in seen_residues:
                continue  # nearest atom of each residue is its entry point
            seen_residues.add(residue)
            result = _best_chain(graph, _node(row), cofactor_nodes, float(row["distance_A"]))
            if result is not None:
                candidates.append((result, row))
        candidates.sort(key=lambda item: (item[0][2], round(item[0][3], 6), item[1]["residue_number"]))

        for rank, ((path, parents, hops, length), entry) in enumerate(candidates[:max_chains], start=1):
            chain_id = f"{site_id}:chain_{rank}"
            shells = [info[n]["shell"] for n in path]
            summaries.append({
                "structure_id": structure_id, "site_id": site_id, "chain_id": chain_id,
                "rank": rank, "entry_residue": f'{entry["residue_name"]}{entry["residue_number"]}',
                "entry_atom": entry["atom_name"], "entry_shell": entry["shell"],
                "terminal_atom": path[-1][-1], "n_hops": hops,
                "total_length_A": round(length, 3),
                "shells_traversed": ">".join(dict.fromkeys(shells)),
            })
            steps.append({"structure_id": structure_id, "site_id": site_id, "chain_id": chain_id,
                          "step": 0, "node_type": POINT_NODE, "atom_name": "SUBSTRATE_POINT",
                          "x": xyz[0], "y": xyz[1], "z": xyz[2]})
            for step, node in enumerate(path, start=1):
                data = info[node]
                if step == 1:
                    kind, distance = "point_contact", float(entry["distance_A"])
                else:
                    _, kind, distance = parents[node]
                steps.append({
                    "structure_id": structure_id, "site_id": site_id, "chain_id": chain_id,
                    "step": step,
                    "node_type": "cofactor" if node in cofactor_nodes else "network",
                    "residue_name": data["residue_name"], "residue_number": data["residue_number"],
                    "chain": data["chain"], "insertion_code": data["insertion_code"],
                    "model_id": data["model_id"], "atom_name": data["atom_name"],
                    "motif": data["motif"], "shell": data["shell"],
                    "edge_from_previous": kind, "edge_distance_A": round(distance, 3),
                })

    steps_frame = pd.DataFrame(steps, columns=STEP_COLUMNS + ["x", "y", "z"])
    # The point row has no residue, which would otherwise turn numbers into floats.
    for column in ("residue_number", "step"):
        try:
            steps_frame[column] = steps_frame[column].astype("Int64")
        except (TypeError, ValueError):
            pass
    return {
        "seed": pd.DataFrame(seeds, columns=["structure_id", "site_id", "x", "y", "z", "contact_cutoff_A"]),
        "contacts": pd.DataFrame(contact_rows, columns=CONTACT_COLUMNS),
        "chains": pd.DataFrame(summaries, columns=SUMMARY_COLUMNS),
        "steps": steps_frame,
    }


def analyze_substrate_seed(
    structure_path,
    cofactor_resname: Union[str, Sequence[str]],
    point: Point,
    *,
    contact_cutoff: float = DEFAULT_CONTACT_CUTOFF_A,
    max_chains: int = 10,
    shells: int = 3,
    **analysis_options: Any,
) -> Dict[str, pd.DataFrame]:
    """Run the coordination analysis, then trace chains from ``point``.

    Returns the usual ``residues``/``atoms``/``links``/``contacts`` tables plus
    ``seed``, ``substrate_contacts``, ``chains`` and ``steps``. Per-site
    analysis is the default so a point is resolved against each cofactor copy
    separately; pass ``site_mode="union"`` to treat the structure as one site.
    """
    options = {"site_mode": "per-site", "direct_coordination": True, **analysis_options}
    tables = analyze_structure(structure_path, cofactor_resname, shells=shells, **options)
    result = substrate_chains(
        tables, point, contact_cutoff=contact_cutoff, max_chains=max_chains
    )
    return {
        **tables,
        "seed": result["seed"],
        "substrate_contacts": result["contacts"],
        "chains": result["chains"],
        "steps": result["steps"],
    }
