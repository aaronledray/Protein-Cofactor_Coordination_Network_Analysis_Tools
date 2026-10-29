"""Water-aware hydrogen-bond wires across a set of structures, voted by sequence.

:func:`conserved_interactions.relay_steps` asks how many links separate a source (for example the
propionate oxygens of a heme) from each residue. This module asks which *wires* recur: the specific
path through residues and waters, and how many independent paths there are.

For each structure a directed relay graph is built from the tidy tables of
:mod:`modules.conserved_interactions`. Nodes are the source, each modeled water, and each residue
(a shared ``position``). An edge joins two nodes whose N/O atoms (or the source atoms, or a water
oxygen) are within ``link``; a residue can be reached through any N/O atom but passes the wire onward
only if it was reached through a side-chain N/O atom (``("r", position)``); a residue reached only
through its backbone is a terminal node (``("rb", position)``), as in
:func:`conserved_interactions.relay_steps`. Three quantities are voted across sequences with the
one-vote-per-sequence rule:

* **reach**: the fewest links from the source to a target position;
* **path signature**: the shortest paths as tuples of positions with ``"w"`` for a water, for example
  ``("w", 404)`` is source to water to residue 404; a signature is conserved when most sequences of every
  group show it among their shortest paths;
* **independent paths**: the number of vertex-disjoint paths from the source to the target (each residue
  and each water used at most once), a measure of how much the connection depends on a single water
  or side chain.

The graph is a hydrogen-bond-capable geometry on heavy atoms, not a proton path: no hydrogens,
donor/acceptor roles, or protonation states are used. A wire that is conserved here is a conserved
connectivity.
"""

from collections import deque
from typing import Dict, Hashable, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from .conserved_interactions import BACKBONE, _close_pairs, _heavy, sequence_vote

SOURCE = ("src",)
MAX_SIGNATURES = 200


def relay_graph(atoms: pd.DataFrame, sources: pd.DataFrame, waters: Optional[pd.DataFrame] = None, *, link: float = 3.5) -> Dict[Hashable, Set[Hashable]]:
    """Directed relay graph of one structure: ``{node: {successors}}``.

    ``atoms``, ``sources`` and ``waters`` hold a single structure. Nodes are ``("src",)``,
    ``("w", i)``, ``("r", position)`` (reached through a side chain; relays onward) and
    ``("rb", position)`` (reached only through the backbone; terminal).
    """
    polar = _heavy(atoms)
    polar = polar[polar["element"].str.upper().isin(["N", "O"])].reset_index(drop=True)
    adj: Dict[Hashable, Set[Hashable]] = {SOURCE: set()}
    if polar.empty:
        return adj
    pxyz = polar[["x", "y", "z"]].to_numpy(float)
    ppos = polar["position"].to_numpy()
    side = np.array([n not in BACKBONE for n in polar["atom"]])
    wxyz = waters[["x", "y", "z"]].to_numpy(float) if waters is not None and len(waters) else np.empty((0, 3))
    sxyz = sources[["x", "y", "z"]].to_numpy(float)
    for p in np.unique(ppos):
        adj[("r", p)] = set()
        adj[("rb", p)] = set()
    for k in range(len(wxyz)):
        adj[("w", k)] = set()

    def enter(b):
        return ("r", ppos[b]) if side[b] else ("rb", ppos[b])

    def add(u, v):
        adj.setdefault(u, set()).add(v)

    i, j, _ = _close_pairs(sxyz, pxyz, link)
    for b in np.unique(j):
        add(SOURCE, enter(b))
    if len(wxyz):
        i, j, _ = _close_pairs(sxyz, wxyz, link)
        for k in np.unique(j):
            add(SOURCE, ("w", int(k)))
        i, j, _ = _close_pairs(wxyz, wxyz, link)
        for a, b in zip(i, j):
            if a != b:
                add(("w", int(a)), ("w", int(b)))
        i, j, _ = _close_pairs(wxyz, pxyz, link)
        for a, b in zip(i, j):
            add(("w", int(a)), enter(b))
            if side[b]:
                add(("r", ppos[b]), ("w", int(a)))
    i, j, _ = _close_pairs(pxyz[side], pxyz, link)
    side_pos = ppos[side]
    for a, b in zip(i, j):
        if side_pos[a] != ppos[b]:
            add(("r", side_pos[a]), enter(b))
    return adj


def _shortest(adj, targets: Set, max_steps: int):
    """Layered BFS from the source: distance and shortest-path predecessors of every reached node."""
    dist, pred = {SOURCE: 0}, {SOURCE: []}
    queue = deque([SOURCE])
    while queue:
        u = queue.popleft()
        if dist[u] >= max_steps:
            continue
        for v in adj.get(u, ()):
            if v not in dist:
                dist[v], pred[v] = dist[u] + 1, [u]
                queue.append(v)
            elif dist[v] == dist[u] + 1:
                pred[v].append(u)
    return dist, pred


def shortest_signatures(adj, targets: Iterable, *, max_steps: int = 6) -> Dict[object, Tuple[int, List[tuple]]]:
    """For each target position, the fewest links and the shortest paths as signatures.

    A signature lists the nodes after the source, with ``"w"`` for a water and the position for a
    residue. Up to ``MAX_SIGNATURES`` signatures are kept per target.
    """
    dist, pred = _shortest(adj, set(), max_steps)
    out: Dict[object, Tuple[int, List[tuple]]] = {}
    for t in targets:
        reached = [n for n in (("r", t), ("rb", t)) if n in dist]
        if not reached:
            continue
        best = min(dist[n] for n in reached)
        starts = [n for n in reached if dist[n] == best]
        sigs: List[tuple] = []

        def walk(v, tail):
            if len(sigs) >= MAX_SIGNATURES:
                return
            if v == SOURCE:
                sigs.append(tuple(reversed(tail)))
                return
            for u in pred[v]:
                walk(u, tail + ["w" if v[0] == "w" else v[1]])

        for node in starts:
            walk(node, [])
        out[t] = (best, sorted(set(sigs), key=lambda s: (len(s), str(s))))
    return out


def disjoint_path_count(adj, target, *, cap: int = 8) -> int:
    """Vertex-disjoint paths from the source to ``("r", target)`` (each residue and water used once)."""
    goal = ("goal", target)
    capacity: Dict[tuple, int] = {}
    graph: Dict[tuple, Set[tuple]] = {}

    def arc(a, b, c):
        capacity[(a, b)] = capacity.get((a, b), 0) + c
        capacity.setdefault((b, a), 0)
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)

    # node splitting: (v, 0) -> (v, 1) with capacity 1 (unlimited for the source and the target residue); the two residue nodes
    # of the target (side-chain and backbone entry) both feed the goal
    for u, succ in adj.items():
        arc((u, 0), (u, 1), cap if u in (SOURCE, ("r", target), ("rb", target)) else 1)
        for v in succ:
            arc((u, 1), (v, 0), 1)
    for node in (("r", target), ("rb", target)):
        if (node, 1) in graph:
            arc((node, 1), (goal, 0), cap)
    arc((goal, 0), (goal, 1), cap)
    if not graph.get((goal, 0), set()) - {(goal, 1)}:
        return 0
    flow = 0
    while flow < cap:
        parent = {(SOURCE, 0): None}
        queue = deque([(SOURCE, 0)])
        while queue and (goal, 1) not in parent:
            a = queue.popleft()
            for b in graph.get(a, ()):
                if b not in parent and capacity[(a, b)] > 0:
                    parent[b] = a
                    queue.append(b)
        if (goal, 1) not in parent:
            break
        b = (goal, 1)
        while parent[b] is not None:
            a = parent[b]
            capacity[(a, b)] -= 1
            capacity[(b, a)] += 1
            b = a
        flow += 1
    return flow


def conserved_wires(atoms: pd.DataFrame, structures: pd.DataFrame, sources: pd.DataFrame, waters: Optional[pd.DataFrame],
                    targets: Sequence, *, link: float = 3.5, max_steps: int = 6, groups: Optional[Sequence[str]] = None,
                    support: float = 0.90, min_sequences: int = 5, structure_fraction: float = 0.5):
    """Reach, shortest-path signatures, and independent paths from a source to target positions.

    Returns ``(reach, signatures, paths)``: ``reach`` has per target and group the fraction of
    sequences reaching it within each number of links (``sequence_vote`` columns for flags
    ``within_k``); ``signatures`` has per target and signature the group fractions and
    ``conserved_present``; ``paths`` has per structure and target the fewest links and the
    vertex-disjoint path count. A structure is evaluated only when it maps the target position.
    """
    meta = structures[["structure_id", "group", "sequence_id"]].drop_duplicates("structure_id")
    steps_rows, sig_rows = [], []
    for sid in meta["structure_id"]:
        a = atoms[atoms["structure_id"] == sid]
        present = set(a["position"])
        src = sources[sources["structure_id"] == sid]
        if a.empty or src.empty:
            continue
        w = waters[waters["structure_id"] == sid] if waters is not None else None
        adj = relay_graph(a, src, w, link=link)
        found = shortest_signatures(adj, [t for t in targets if t in present], max_steps=max_steps)
        for t in targets:
            if t not in present:
                continue
            steps, sigs = found.get(t, (np.nan, []))
            steps_rows.append({"structure_id": sid, "target": t, "steps": steps, "disjoint_paths": disjoint_path_count(adj, t) if not np.isnan(steps) else 0})
            for s in sigs:
                sig_rows.append({"structure_id": sid, "target": t, "signature": s})
    paths = pd.DataFrame(steps_rows, columns=["structure_id", "target", "steps", "disjoint_paths"])
    for k in range(1, max_steps + 1):
        paths[f"within_{k}"] = paths["steps"] <= k
    reach = sequence_vote(paths, structures, ["target"], [f"within_{k}" for k in range(1, max_steps + 1)], groups=groups, support=support,
                          min_sequences=min_sequences, structure_fraction=structure_fraction)
    observed = pd.DataFrame(sig_rows, columns=["structure_id", "target", "signature"])
    if observed.empty:
        return reach, pd.DataFrame(), paths
    candidates = observed[["target", "signature"]].drop_duplicates()
    have = set(zip(observed["structure_id"], observed["target"], observed["signature"]))
    cand_by_target = candidates.groupby("target")["signature"].apply(list).to_dict()
    rows = [{"structure_id": sid, "target": t, "signature": s, "present": (sid, t, s) in have}
            for sid, t in set(zip(paths["structure_id"], paths["target"])) for s in cand_by_target.get(t, [])]
    sig_obs = pd.DataFrame(rows)
    sig_obs["signature"] = sig_obs["signature"].astype(str)
    signatures = sequence_vote(sig_obs, structures, ["target", "signature"], ["present"], groups=groups, support=support,
                               min_sequences=min_sequences, structure_fraction=structure_fraction)
    return reach, signatures, paths
