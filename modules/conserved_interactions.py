"""Conserved interactions across a set of structures.

Given structures that share a residue correspondence (``position`` labels from an
alignment or a structural superposition), this module asks which interactions recur
in nearly every sequence of every group, and what chemistry they have. It is the
reusable form of the analysis in the heme-copper oxidase comparison and shares its
rules, so results can be checked against it.

The pieces, all pure functions on tidy tables (:func:`load_atoms` builds the atom table from a
structure file and a position map):

* :func:`residue_contacts` -- residue-residue contacts per structure with typed
  chemistry (closest-atom class, polar-capable N/O pair, hydrophobic, ring-atom,
  ionic contact);
* :func:`cofactor_contacts` and :func:`water_bridges` -- contacts between residues and
  named cofactor features, and residues bridged to a feature by one modeled water;
* :func:`relay_steps` -- fewest hydrogen-bond links from a feature to each residue
  through waters and side-chain N/O atoms;
* :func:`sequence_vote` -- the conservation rule: one vote per sequence, then support
  required in every group;
* :func:`anchor_audit` -- flags structures whose correspondence puts a known anchor
  residue far from its site (a misaligned region), which would otherwise count as
  missing contacts.

Input tables
------------
``atoms``: ``structure_id, position, seq_index, resname, atom, element, x, y, z`` (one
row per heavy atom; ``seq_index`` is the residue's order in its own chain, used to
exclude sequence neighbors). ``features``: ``structure_id, feature, atom, element, kind,
x, y, z`` with ``kind`` in ``{"metal", "oxygen"}``. ``waters``: ``structure_id, x, y, z``.
``structures``: ``structure_id, group, sequence_id``.

Geometry is judged from heavy atoms; no hydrogens are used, so "polar" means a
hydrogen-bond-capable N/O pair, not a hydrogen bond (see :mod:`modules.protonation`).
"""

from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

BACKBONE = frozenset({"N", "CA", "C", "O", "OXT"})
RING_ATOMS: Dict[str, frozenset] = {
    "PHE": frozenset({"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
    "TYR": frozenset({"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
    "TRP": frozenset({"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}),
    "HIS": frozenset({"CG", "ND1", "CD2", "CE1", "NE2"}),
}
POSITIVE_ATOMS = frozenset({("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2")})
NEGATIVE_ATOMS = frozenset({("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")})
RESIDUE_CLASS = {
    **{a: "hydrophobic" for a in "AVLIMPGC"},
    **{a: "aromatic" for a in "FWY"},
    **{a: "polar" for a in "STNQH"},
    **{a: "positive" for a in "KR"},
    **{a: "negative" for a in "DE"},
}
CONTACT_CLASSES = ("bb_bb", "bb_sc", "sc_bb", "sc_sc")
CONTACT_FLAGS = ("polar", "polar_sc", "polar_bb_bb", "hydrophobic", "aromatic", "ionic")
FEATURE_CLASSES = ("coord", "polar", "ionic", "contact")
_BLOCK = 4_000_000


def _close_pairs(left: np.ndarray, right: np.ndarray, cutoff: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Index pairs ``(i, j)`` with ``|left[i] - right[j]| <= cutoff`` and their distances."""
    if len(left) == 0 or len(right) == 0:
        empty = np.empty(0, dtype=int)
        return empty, empty, np.empty(0)
    step = max(1, _BLOCK // len(right))
    ii, jj, dd = [], [], []
    for start in range(0, len(left), step):
        block = left[start:start + step]
        d = np.sqrt(((block[:, None, :] - right[None, :, :]) ** 2).sum(-1))
        i, j = np.nonzero(d <= cutoff)
        ii.append(i + start)
        jj.append(j)
        dd.append(d[i, j])
    return np.concatenate(ii), np.concatenate(jj), np.concatenate(dd)


def _heavy(atoms: pd.DataFrame) -> pd.DataFrame:
    return atoms[~atoms["element"].str.upper().isin(["H", "D"])]


def _per_structure(table: pd.DataFrame):
    for sid, frame in table.groupby("structure_id", sort=False):
        yield sid, frame


# ----------------------------------------------------------------------------- loading
def load_atoms(path: str, chain_id: str, position_map: Optional[Mapping[Tuple[int, str], object]] = None, *,
               structure_id: Optional[str] = None, model_index: int = 0) -> pd.DataFrame:
    """Heavy atoms of one chain as the tidy ``atoms`` table.

    ``position_map`` maps ``(residue_number, insertion_code)`` to the shared position label
    (for example an alignment column); residues without an entry are skipped, and with no
    map the residue number is the position. ``seq_index`` is the residue's order among the
    chain's amino-acid residues. Altlocs other than blank or ``A`` and hydrogens are dropped.
    """
    from .io_utils import unpack_pdb_file

    structure, _ = unpack_pdb_file(str(path))
    model = list(structure)[model_index]
    sid = structure_id or Path(path).stem
    rows = []
    order = 0
    for residue in model[chain_id]:
        if residue.id[0] != " ":
            continue
        order += 1
        key = (int(residue.id[1]), residue.id[2].strip())
        position = key[0] if position_map is None else position_map.get(key)
        if position is None:
            continue
        # get_unpacked_list() yields each conformer of a disordered atom; iterating the residue
        # would give one DisorderedAtom whose altloc is that of the highest-occupancy conformer,
        # so the atom would vanish whenever B outweighs A.
        for atom in residue.get_unpacked_list():
            if atom.element.upper() in ("H", "D") or atom.get_altloc() not in (" ", "A"):
                continue
            x, y, z = (float(v) for v in atom.coord)
            rows.append((sid, position, order, residue.resname.strip().upper(), atom.get_id(), atom.element.upper(), x, y, z))
    return pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])


# ----------------------------------------------------------------------------- contacts
def residue_contacts(atoms: pd.DataFrame, *, cutoff: float = 4.0, min_distance: float = 1.9,
                     min_separation: int = 2, chemistry: bool = True) -> pd.DataFrame:
    """Residue-residue contacts per structure with typed chemistry.

    A pair of residues is in contact when its closest heavy atoms are between
    ``min_distance`` (closer pairs are covalent candidates) and ``cutoff`` and the
    residues are at least ``min_separation`` apart in ``seq_index``. Chemistry flags use
    the thresholds below, capped at ``cutoff``: ``polar`` (N/O pair within 3.5 Å),
    ``polar_sc`` (at least one side-chain atom), ``polar_bb_bb`` (both backbone),
    ``hydrophobic`` (side-chain C/S pair within 4.0 Å), ``aromatic`` (a ring atom of
    either residue within 4.0 Å of any heavy atom of the other), ``ionic`` (Lys/Arg
    N to Asp/Glu O within 4.0 Å); ``closest`` is the backbone/side-chain class of the
    single closest atom pair.
    """
    columns = ["structure_id", "position_a", "position_b", "resname_a", "resname_b", "min_distance", "n_contacts"]
    if chemistry:
        columns += ["closest", *CONTACT_FLAGS]
    frames: List[pd.DataFrame] = []
    for sid, frame in _per_structure(_heavy(atoms)):
        frame = frame.reset_index(drop=True)
        xyz = frame[["x", "y", "z"]].to_numpy(float)
        i, j, d = _close_pairs(xyz, xyz, cutoff)
        pos = frame["position"].to_numpy()
        seq = frame["seq_index"].to_numpy()
        keep = (i < j) & (pos[i] != pos[j]) & (np.abs(seq[i] - seq[j]) >= min_separation)
        i, j, d = i[keep], j[keep], d[keep]
        if len(i) == 0:
            continue
        swap = pos[i] > pos[j]
        a, b = np.where(swap, j, i), np.where(swap, i, j)
        name = frame["atom"].to_numpy()
        elem = frame["element"].str.upper().to_numpy()
        res = frame["resname"].str.upper().to_numpy()
        bb = np.array([n in BACKBONE for n in name])
        pair = pd.DataFrame({"position_a": pos[a], "position_b": pos[b], "d": d, "closest": np.where(bb[a], "bb", "sc").astype(object) + "_" + np.where(bb[b], "bb", "sc").astype(object)})
        if chemistry:
            polar = np.isin(elem[a], ["N", "O"]) & np.isin(elem[b], ["N", "O"]) & (d <= min(3.5, cutoff))
            ring = np.array([n in RING_ATOMS.get(r, ()) for n, r in zip(name, res)])
            posi = np.array([(r, n) in POSITIVE_ATOMS for n, r in zip(name, res)])
            nega = np.array([(r, n) in NEGATIVE_ATOMS for n, r in zip(name, res)])
            within = d <= min(4.0, cutoff)
            pair["polar"] = polar
            pair["polar_sc"] = polar & (~bb[a] | ~bb[b])
            pair["polar_bb_bb"] = polar & bb[a] & bb[b]
            pair["hydrophobic"] = np.isin(elem[a], ["C", "S"]) & np.isin(elem[b], ["C", "S"]) & ~bb[a] & ~bb[b] & within
            pair["aromatic"] = (ring[a] | ring[b]) & within
            pair["ionic"] = ((posi[a] & nega[b]) | (nega[a] & posi[b])) & within
        pair["ra"], pair["rb"] = res[a], res[b]
        keys = ["position_a", "position_b"]
        flags = list(CONTACT_FLAGS) if chemistry else []
        closest = pair.sort_values("d", kind="stable").drop_duplicates(keys)[keys + ["d", "closest", "ra", "rb"]]
        agg = pair.groupby(keys, sort=False).agg(n_contacts=("d", "size"), **{f: (f, "any") for f in flags}).reset_index()
        out = closest.merge(agg, on=keys)
        out = out[out["d"] >= min_distance].rename(columns={"d": "min_distance", "ra": "resname_a", "rb": "resname_b"})
        out.insert(0, "structure_id", sid)
        frames.append(out)
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True)[columns]


# ----------------------------------------------------------------------------- cofactor features
def cofactor_contacts(atoms: pd.DataFrame, features: pd.DataFrame, *, coordination_distance: float = 2.8,
                      polar_distance: float = 3.5, ionic_distance: float = 4.0, contact_distance: float = 4.0,
                      record_distance: float = 8.0) -> pd.DataFrame:
    """Contacts between residues and cofactor features, one row per structure x feature x position.

    Classes: ``coord`` (a ``metal`` feature atom to a residue N, O, or S within
    ``coordination_distance``), ``polar`` (an ``oxygen`` feature atom to a residue N or O
    within ``polar_distance``), ``ionic`` (an ``oxygen`` feature atom to Lys NZ or Arg
    NE/NH1/NH2 within ``ionic_distance``), ``contact`` (any heavy-atom pair within
    ``contact_distance``). A feature groups all its atoms; residues farther than
    ``record_distance`` from every atom of the feature are not recorded (absent).
    """
    rows: List[dict] = []
    heavy = _heavy(atoms)
    for sid, frame in _per_structure(heavy):
        feats = features[features["structure_id"] == sid]
        if feats.empty:
            continue
        frame = frame.reset_index(drop=True)
        xyz = frame[["x", "y", "z"]].to_numpy(float)
        elem = frame["element"].str.upper().to_numpy()
        name = frame["atom"].to_numpy()
        res = frame["resname"].str.upper().to_numpy()
        pos = frame["position"].to_numpy()
        ionic_atom = np.array([(r, n) in POSITIVE_ATOMS for r, n in zip(res, name)])
        for feature, fgrp in feats.groupby("feature", sort=False):
            fxyz = fgrp[["x", "y", "z"]].to_numpy(float)
            metal = (fgrp["kind"] == "metal").to_numpy()
            i, j, d = _close_pairs(xyz, fxyz, record_distance)
            if len(i) == 0:
                continue
            table = pd.DataFrame({"i": i, "j": j, "d": d, "pos": pos[i]})
            for p, grp in table.groupby("pos", sort=False):
                ii, jj, dd = grp["i"].to_numpy(), grp["j"].to_numpy(), grp["d"].to_numpy()
                is_metal = metal[jj]
                row = {"structure_id": sid, "feature": feature, "position": p, "resname": res[ii[0]], "min_distance": float(dd.min()),
                       "coord": bool((is_metal & np.isin(elem[ii], ["N", "O", "S"]) & (dd <= coordination_distance)).any()),
                       "polar": bool((~is_metal & np.isin(elem[ii], ["N", "O"]) & (dd <= polar_distance)).any()),
                       "ionic": bool((~is_metal & ionic_atom[ii] & (dd <= ionic_distance)).any()),
                       "contact": bool((dd <= contact_distance).any())}
                rows.append(row)
    return pd.DataFrame(rows, columns=["structure_id", "feature", "position", "resname", "min_distance", *FEATURE_CLASSES])


def water_bridges(atoms: pd.DataFrame, features: pd.DataFrame, waters: pd.DataFrame, *, bridge: float = 3.5,
                  direct: float = 3.5) -> pd.DataFrame:
    """Residues bridged to an oxygen feature by one modeled water.

    A residue is bridged when a water oxygen lies within ``bridge`` of a feature atom and
    within ``bridge`` of a residue N or O atom; ``direct`` is a residue N/O within
    ``direct`` of the feature; ``indirect`` is bridged and not direct.
    """
    rows: List[dict] = []
    heavy = _heavy(atoms)
    for sid, frame in _per_structure(heavy):
        w = waters[waters["structure_id"] == sid][["x", "y", "z"]].to_numpy(float)
        feats = features[(features["structure_id"] == sid) & (features["kind"] == "oxygen")]
        polar = frame[frame["element"].str.upper().isin(["N", "O"])].reset_index(drop=True)
        if polar.empty:
            continue
        pxyz = polar[["x", "y", "z"]].to_numpy(float)
        ppos = polar["position"].to_numpy()
        for feature, fgrp in feats.groupby("feature", sort=False):
            fxyz = fgrp[["x", "y", "z"]].to_numpy(float)
            near = np.zeros(len(w), bool)
            if len(w):
                _, jj, _ = _close_pairs(fxyz, w, bridge)
                near[np.unique(jj)] = True
            bridged, direct_set = set(), set()
            if near.any():
                _, pj, _ = _close_pairs(w[near], pxyz, bridge)
                bridged = set(ppos[np.unique(pj)])
            _, dj, _ = _close_pairs(fxyz, pxyz, direct)
            direct_set = set(ppos[np.unique(dj)])
            for p in np.unique(ppos):
                rows.append({"structure_id": sid, "feature": feature, "position": p, "bridged": p in bridged,
                             "direct": p in direct_set, "indirect": (p in bridged) and (p not in direct_set)})
    return pd.DataFrame(rows, columns=["structure_id", "feature", "position", "bridged", "direct", "indirect"])


def relay_steps(atoms: pd.DataFrame, sources: pd.DataFrame, waters: Optional[pd.DataFrame] = None, *,
                link: float = 3.5, max_steps: int = 6) -> pd.DataFrame:
    """Fewest hydrogen-bond links from source atoms to each residue.

    Nodes are residue N/O atoms and (optionally) modeled waters; two nodes are linked
    when closer than ``link``. A residue reached through a side-chain atom relays
    onward; one reached only through its backbone N or O is a terminal. ``sources`` has
    ``structure_id, x, y, z`` (for example the propionate oxygens). Returns
    ``structure_id, position, steps`` for the residues reached within ``max_steps``.
    """
    out: List[dict] = []
    heavy = _heavy(atoms)
    for sid, frame in _per_structure(heavy):
        src = sources[sources["structure_id"] == sid][["x", "y", "z"]].to_numpy(float)
        polar = frame[frame["element"].str.upper().isin(["N", "O"])].reset_index(drop=True)
        if polar.empty or len(src) == 0:
            continue
        axyz = polar[["x", "y", "z"]].to_numpy(float)
        owner = polar["position"].to_numpy()
        side = np.array([n not in BACKBONE for n in polar["atom"]])
        wxyz = waters[waters["structure_id"] == sid][["x", "y", "z"]].to_numpy(float) if waters is not None else np.empty((0, 3))
        by_pos: Dict[object, List[int]] = {}
        for k, (p, s) in enumerate(zip(owner, side)):
            if s:
                by_pos.setdefault(p, []).append(k)

        def near(point):
            ai, _, _ = _close_pairs(axyz, point[None, :], link)
            wi, _, _ = _close_pairs(wxyz, point[None, :], link) if len(wxyz) else (np.empty(0, int), None, None)
            return {("a", int(x)) for x in ai} | {("w", int(x)) for x in wi}

        level = set()
        for s in src:
            level |= near(s)
        best: Dict[object, int] = {}
        done_w, done_r = set(), set()
        for step in range(1, max_steps + 1):
            nxt = set()
            for kind, k in level:
                if kind == "w":
                    if k not in done_w:
                        done_w.add(k)
                        nxt |= near(wxyz[k])
                else:
                    p = owner[k]
                    best.setdefault(p, step)
                    if side[k] and p not in done_r:
                        done_r.add(p)
                        for m in by_pos[p]:
                            nxt |= near(axyz[m])
            level = nxt
            if not level:
                break
        out.extend({"structure_id": sid, "position": p, "steps": s} for p, s in best.items())
    return pd.DataFrame(out, columns=["structure_id", "position", "steps"])


# ----------------------------------------------------------------------------- the vote
def sequence_vote(observations: pd.DataFrame, structures: pd.DataFrame, keys: Sequence[str], flags: Sequence[str], *,
                  condition: Optional[str] = None, structure_fraction: float = 0.5, support: float = 0.90,
                  min_sequences: int = 5, groups: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """One vote per sequence, then support in every group.

    ``observations`` has ``structure_id``, the ``keys`` columns, and boolean ``flags``
    columns. Within a (key, sequence) the vote is positive when at least
    ``structure_fraction`` of the sequence's structures (those present in ``observations``
    for that key, and where ``condition`` is true when given) are positive. A key is
    ``conserved_<flag>`` when every group has at least ``min_sequences`` voting sequences
    and at least ``support`` of them positive. Returns one row per key with
    ``n_<group>`` and ``frac_<group>`` for each flag (``<flag>_n_<group>``, ``<flag>_<group>``)
    and ``conserved_<flag>``.
    """
    meta = structures[["structure_id", "group", "sequence_id"]].drop_duplicates("structure_id")
    data = observations.merge(meta, on="structure_id", how="inner")
    groups = list(groups) if groups is not None else sorted(meta["group"].unique())
    keys = list(keys)
    result: Optional[pd.DataFrame] = None
    for flag in flags:
        sub = data if condition is None else data[data[condition].astype(bool)]
        share = sub.assign(_flag=sub[flag].astype(float)).groupby(keys + ["group", "sequence_id"], sort=False)["_flag"].mean()
        per_seq = (share >= structure_fraction).rename("vote").reset_index()
        agg = per_seq.groupby(keys + ["group"], sort=False)["vote"].agg(["sum", "count"]).reset_index()
        wide = agg.pivot_table(index=keys, columns="group", values=["sum", "count"], aggfunc="first")
        frame = pd.DataFrame(index=wide.index)
        conserved = pd.Series(True, index=wide.index)
        for g in groups:
            n = wide[("count", g)] if ("count", g) in wide.columns else pd.Series(0, index=wide.index)
            pos = wide[("sum", g)] if ("sum", g) in wide.columns else pd.Series(0, index=wide.index)
            n, pos = n.fillna(0).astype(int), pos.fillna(0)
            frame[f"{flag}_n_{g}"] = n
            frame[f"{flag}_{g}"] = (pos / n.where(n > 0)).round(3)
            conserved &= (n >= min_sequences) & (pos / n.where(n > 0) >= support)
        frame[f"conserved_{flag}"] = conserved.astype(bool)
        frame = frame.reset_index()
        result = frame if result is None else result.merge(frame, on=keys, how="outer")
    return result if result is not None else pd.DataFrame(columns=keys)


def conserved_contacts(atoms: pd.DataFrame, structures: pd.DataFrame, *, groups: Optional[Sequence[str]] = None,
                       cutoff: float = 4.0, min_separation: int = 2, support: float = 0.90, min_sequences: int = 5,
                       structure_fraction: float = 0.5) -> pd.DataFrame:
    """Contacts passing the sequence vote, with the conserved chemistry of each.

    Every pair of positions mapped in a structure is a candidate; the contact vote counts
    a structure as negative when both positions are present but not in contact. Chemistry
    flags are voted over contact-positive structures only (``closest`` classes and the
    flags in :data:`CONTACT_FLAGS`).
    """
    contacts = residue_contacts(atoms, cutoff=cutoff, min_separation=min_separation)
    present = _heavy(atoms).groupby("structure_id")["position"].apply(set)
    # candidate pairs: contacts anywhere in the set, evaluated in every structure that maps both positions
    candidates = list(contacts[["position_a", "position_b"]].drop_duplicates().itertuples(index=False, name=None))
    in_contact = set(zip(contacts["structure_id"], contacts["position_a"], contacts["position_b"]))
    rows = [(sid, a, b, (sid, a, b) in in_contact) for sid, have in present.items() for a, b in candidates if a in have and b in have]
    obs = pd.DataFrame(rows, columns=["structure_id", "position_a", "position_b", "contact"])
    obs = obs.merge(contacts, on=["structure_id", "position_a", "position_b"], how="left")
    obs["contact"] = obs["contact"].astype(bool)
    for c in CONTACT_CLASSES:
        obs[c] = obs["closest"] == c
    keys = ["position_a", "position_b"]
    base = sequence_vote(obs, structures, keys, ["contact"], groups=groups, support=support,
                         min_sequences=min_sequences, structure_fraction=structure_fraction)
    chem = sequence_vote(obs, structures, keys, list(CONTACT_CLASSES) + list(CONTACT_FLAGS), condition="contact", groups=groups,
                         support=support, min_sequences=min_sequences, structure_fraction=structure_fraction)
    return base.merge(chem, on=keys, how="left")


# ----------------------------------------------------------------------------- audits
def anchor_audit(atoms: pd.DataFrame, sites: pd.DataFrame, anchor_position, structures: Optional[pd.DataFrame] = None,
                 *, max_distance: float = 8.0) -> pd.DataFrame:
    """Flag structures whose correspondence misplaces an anchor residue.

    ``sites`` has ``structure_id, x, y, z`` (the site the anchor should be next to, for
    example the heme iron whose axial histidine defines the anchor position). A structure
    fails when its residue at ``anchor_position`` is missing or has no atom within
    ``max_distance`` of the site. Failing sequences would otherwise count as negative for
    every contact in the misaligned region.
    """
    rows = []
    heavy = _heavy(atoms)
    for sid, site in sites.groupby("structure_id", sort=False):
        s = site[["x", "y", "z"]].to_numpy(float)
        res = heavy[(heavy["structure_id"] == sid) & (heavy["position"] == anchor_position)]
        if res.empty:
            rows.append({"structure_id": sid, "anchor_resname": None, "distance": np.nan, "ok": False})
            continue
        _, _, d = _close_pairs(res[["x", "y", "z"]].to_numpy(float), s, 1e9)
        rows.append({"structure_id": sid, "anchor_resname": res["resname"].iloc[0], "distance": float(d.min()), "ok": bool(d.min() <= max_distance)})
    table = pd.DataFrame(rows, columns=["structure_id", "anchor_resname", "distance", "ok"])
    if structures is not None:
        table = table.merge(structures[["structure_id", "group", "sequence_id"]].drop_duplicates("structure_id"), on="structure_id", how="left")
    return table


def misaligned_sequences(audit: pd.DataFrame) -> List[str]:
    """Sequence ids with at least one failing structure in an :func:`anchor_audit` table."""
    if "sequence_id" not in audit.columns:
        raise ValueError("pass structures to anchor_audit so sequence ids are available")
    return sorted(audit.loc[~audit["ok"], "sequence_id"].dropna().unique().tolist())


def chemistry_summary(votes: pd.DataFrame, *, local_separation: Optional[pd.Series] = None) -> pd.DataFrame:
    """Counts of conserved contacts by conserved chemistry class (optionally split local / distant).

    ``votes`` is the output of :func:`conserved_contacts`; ``local_separation`` is an
    optional boolean series (True = local) aligned to ``votes`` that adds the split.
    """
    cols = [c for c in votes.columns if c.startswith("conserved_") and c != "conserved_contact"]
    passing = votes[votes["conserved_contact"]]
    rows = []
    for c in cols:
        row = {"class": c.replace("conserved_", ""), "all": int(passing[c].fillna(False).sum())}
        if local_separation is not None:
            loc = local_separation.loc[passing.index].astype(bool)
            row["local"] = int(passing.loc[loc, c].fillna(False).sum())
            row["distant"] = int(passing.loc[~loc, c].fillna(False).sum())
        rows.append(row)
    return pd.DataFrame(rows)
