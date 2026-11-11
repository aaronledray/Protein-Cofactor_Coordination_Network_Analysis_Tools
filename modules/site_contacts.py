"""Atom-level contacts of chosen source atoms, with every alternate conformation kept (experimental, opt-in).

``analyze_structure`` builds coordination *shells*: a pruned list (one atom per residue moiety, no carbon seeds) in which only the
highest-occupancy conformer of an atom with alternate locations is visible, and the conformer is not labelled. This module is the
atom-level companion for the question "what lies around these atoms?": it lists every heavy atom (protein, water, other hetero
groups) within a neighbor cutoff of each selected source atom, with the alternate-location label and occupancy of both atoms.
Nothing in the shell analysis, its tables or its defaults is changed or called.

Inputs: a structure file (PDB or mmCIF, ``.gz`` accepted) and ``sources``, a list of selectors. Each selector is a mapping:

``label``    the ``source_group`` name written to the tables (required);
``resname``  residue name, compared upper-cased (give this, ``resnum``, or both);
``atoms``    explicit atom names within the matched residues (default: every heavy atom of them);
``chain``, ``resnum``, ``icode``  optional filters in the file's own numbering (verbatim; no renumbering is ever applied);
``charge``   optional ``"negative"``/``"positive"`` to override the atom's built-in formal-charge class used for ``ionic_candidate``.

Definitions (heavy atoms only; hydrogens and deuteriums are ignored and counted in the manifest):

* *neighbor*: any heavy atom within ``neighbor_cutoff`` (4.5 A) of a source atom, whatever its chemistry; atoms of the source atom's own
  residue are skipped (``include_same_residue=True`` keeps them and sets ``same_residue``);
* ``polar_candidate``: source and partner are both N or O, within ``polar_cutoff`` (3.5 A). It is a heavy-atom N/O pair, **not** a
  hydrogen bond: no hydrogen is placed and no angle is judged;
* ``close_polar``: a ``polar_candidate`` within ``close_polar_cutoff`` (3.2 A);
* ``ionic_candidate``: the two atoms have opposite built-in charge classes within ``ionic_cutoff`` (4.0 A). Negative: Asp OD1/OD2, Glu
  OE1/OE2, HEM propionate oxygens; positive: Lys NZ, Arg NE/NH1/NH2. Histidine protonation is never assumed;
* ``water_bridge``: one modelled water oxygen within ``water_bridge_cutoff`` (3.5 A) of a polar (N/O) source atom and of an N/O atom of
  a different, non-water residue.

Alternate locations: ``altloc_policy="all"`` keeps every conformer of every atom and writes one row per (source conformer, partner
conformer) pair; ``"highest_occupancy"`` keeps the conformer Biopython selects (highest occupancy, the first listed on a tie), which is
what ``analyze_structure`` sees. ``altloc_compatible`` is False when both atoms carry different non-blank location labels, a pairing a
single physical state cannot contain (the labels are only guaranteed consistent within a residue, so the row is kept, not dropped).
With ``compare_highest_occupancy=True`` the residue summary gains ``in_highest_occupancy`` and ``min_distance_highest_occupancy_A``.

Only the deposited asymmetric unit is analysed; symmetry mates and biological-assembly transforms are never applied (see
:mod:`modules.assembly_policy`). Residue names, chains, numbers, insertion codes and hetero flags are reported exactly as read. A table
of contacts describes a crystal structure; it does not choose a quantum region, protonation states, or a model of another species.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import platform
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from .conserved_interactions import NEGATIVE_ATOMS, POSITIVE_ATOMS
from .io_utils import _choose_parser
from .motif_registry import HEME_PROPIONATE_OXYGENS

SITE_CONTACTS_VERSION = "1"
ALTLOC_POLICIES = ("all", "highest_occupancy")
WATER_RESNAMES = frozenset({"HOH", "WAT", "H2O", "DOD", "SOL"})
HYDROGENS = frozenset({"H", "D"})
POLAR_ELEMENTS = frozenset({"N", "O"})
SYMMETRY_NOTE = "asymmetric unit only: crystal-symmetry mates and biological-assembly transforms are not applied"

CONTACT_COLUMNS = [
    "structure_id", "input_sha256", "model_id", "altloc_policy", "source_group",
    "source_chain", "source_resname", "source_resnum", "source_icode", "source_hetero_flag", "source_atom",
    "source_altloc", "source_occupancy",
    "partner_chain", "partner_resname", "partner_resnum", "partner_icode", "partner_hetero_flag", "partner_atom",
    "partner_element", "partner_altloc", "partner_occupancy", "partner_kind",
    "distance_A", "polar_candidate", "close_polar", "ionic_candidate", "same_residue",
    "altloc_compatible", "seq_separation",
]
SUMMARY_COLUMNS = [
    "structure_id", "model_id", "altloc_policy", "source_group", "source_chain", "source_resname", "source_resnum", "source_icode",
    "partner_chain", "partner_resname", "partner_resnum", "partner_icode", "partner_hetero_flag", "partner_kind",
    "min_distance_A", "n_atom_pairs", "source_atoms", "partner_atoms", "partner_altlocs_present", "classes",
]
SUMMARY_COMPARISON_COLUMNS = ["in_highest_occupancy", "min_distance_highest_occupancy_A"]
BRIDGE_COLUMNS = [
    "structure_id", "input_sha256", "model_id", "altloc_policy", "source_group",
    "source_chain", "source_resname", "source_resnum", "source_icode", "source_atom", "source_altloc", "source_occupancy",
    "water_chain", "water_resname", "water_resnum", "water_icode", "water_atom", "water_altloc", "water_occupancy",
    "source_water_distance_A",
    "partner_chain", "partner_resname", "partner_resnum", "partner_icode", "partner_hetero_flag", "partner_atom",
    "partner_element", "partner_altloc", "partner_occupancy", "partner_kind",
    "water_partner_distance_A", "source_partner_distance_A", "direct_polar_candidate", "altloc_compatible",
]
STATUS_COLUMNS = [
    "structure_id", "input_sha256", "source_group", "requested_resname", "requested_atoms", "requested_chain",
    "requested_resnum", "requested_icode", "status", "n_residues_matched", "n_atom_conformers", "matched_residues",
    "missing_atoms", "note",
]
TABLE_NAMES = {
    "site_atom_contacts": CONTACT_COLUMNS,
    "site_residue_summary": SUMMARY_COLUMNS + SUMMARY_COMPARISON_COLUMNS,
    "site_water_bridges": BRIDGE_COLUMNS,
    "site_status": STATUS_COLUMNS,
}


# ----------------------------------------------------------------------------- helpers
def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _structure_id(path: Path) -> str:
    name = path.name
    for suffix in (".gz", ".mmcif", ".cif", ".pdb"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return name


def _element(atom) -> str:
    element = str(getattr(atom, "element", "") or "").strip().upper()
    if element:
        return element
    letters = re.sub(r"[^A-Za-z]", "", atom.get_name())
    return letters[:1].upper()


def _charge_class(resname: str, atom: str) -> str:
    resname = resname.upper()
    if (resname, atom) in NEGATIVE_ATOMS or (resname == "HEM" and atom in HEME_PROPIONATE_OXYGENS):
        return "negative"
    if (resname, atom) in POSITIVE_ATOMS:
        return "positive"
    return ""


def _kind(hetero_id: str, resname: str) -> str:
    if hetero_id == "W" or resname.upper() in WATER_RESNAMES:
        return "water"
    return "protein" if hetero_id == " " else "hetero"


def _opened(path: Path):
    """The Biopython structure of ``path`` (PDB or mmCIF, optionally gzipped)."""
    import gzip
    import tempfile

    name = str(path)
    if name.lower().endswith(".gz"):
        # io_utils' parser choice is by extension; hand it the decompressed text under the inner name.
        inner = name[:-3]
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder) / Path(inner).name
            with gzip.open(path, "rb") as source, open(target, "wb") as sink:
                sink.write(source.read())
            return _choose_parser(str(target)).get_structure("structure", str(target))
    return _choose_parser(name).get_structure("structure", name)


def _load_atoms(structure, models: str) -> tuple:
    """(heavy-atom conformer table, counts). One row per atom conformer; ``selected`` marks Biopython's highest-occupancy choice."""
    rows: List[tuple] = []
    counts = {"n_models": 0, "n_chains": 0, "n_hydrogens_ignored": 0, "n_atoms_with_altlocs": 0, "n_disordered_residues": 0}
    for model_index, model in enumerate(structure):
        if models == "first" and model_index > 0:
            break
        counts["n_models"] += 1
        for chain in model:
            counts["n_chains"] += 1
            for residue in chain:
                disordered_residue = residue.is_disordered() == 2
                if disordered_residue:
                    counts["n_disordered_residues"] += 1
                    variants = [(child, child is residue.selected_child) for child in residue.disordered_get_list()]
                else:
                    variants = [(residue, True)]
                for variant, variant_selected in variants:
                    hetero_id, number, icode = variant.get_id()
                    resname = variant.get_resname().strip()
                    kind = _kind(hetero_id, resname)
                    for atom in variant.get_list():
                        if atom.is_disordered():
                            conformers = atom.disordered_get_list()
                            chosen = atom.selected_child
                            if len(conformers) > 1:
                                counts["n_atoms_with_altlocs"] += 1
                        else:
                            conformers, chosen = [atom], atom
                        for conformer in conformers:
                            element = _element(conformer)
                            if element in HYDROGENS:
                                counts["n_hydrogens_ignored"] += 1
                                continue
                            occupancy = conformer.get_occupancy()
                            x, y, z = (float(v) for v in conformer.coord)
                            rows.append((
                                model.id, chain.id, resname, int(number), str(icode or "").strip(), hetero_id.strip(), kind,
                                conformer.get_name().strip(), element, conformer.get_altloc().strip(),
                                np.nan if occupancy is None else float(occupancy), x, y, z,
                                bool(variant_selected and conformer is chosen),
                            ))
    columns = ["model_id", "chain", "resname", "resnum", "icode", "hetero_flag", "kind", "atom", "element", "altloc",
               "occupancy", "x", "y", "z", "selected"]
    return pd.DataFrame(rows, columns=columns), counts


def _residue_key(frame: pd.DataFrame) -> pd.Series:
    """Residue identity (model, chain, number, insertion code, hetero flag); a point-mutation variant shares it with its sibling."""
    return (frame["model_id"].astype(str) + "|" + frame["chain"].astype(str) + "|" + frame["resnum"].astype(str) + "|"
            + frame["icode"].astype(str) + "|" + frame["hetero_flag"].astype(str))


def _match_sources(atoms: pd.DataFrame, sources: Sequence[Mapping[str, Any]], structure_id: str, input_sha256: str):
    """(source atom table with ``group``/``charge`` columns, status rows)."""
    pieces: List[pd.DataFrame] = []
    status: List[dict] = []
    for spec in sources:
        label = spec.get("label")
        if not label:
            raise ValueError("every source needs a label")
        resname, chain, resnum, icode = spec.get("resname"), spec.get("chain"), spec.get("resnum"), spec.get("icode")
        wanted_atoms = [str(a).strip() for a in spec["atoms"]] if spec.get("atoms") else None
        charge_override = spec.get("charge")
        if charge_override not in (None, "negative", "positive", "none"):
            raise ValueError("charge must be 'negative', 'positive' or 'none'")
        if resname is None and resnum is None:
            raise ValueError(f"source {label!r} needs a resname, a resnum, or both")
        by_number = pd.Series(True, index=atoms.index)
        if chain is not None:
            by_number &= atoms["chain"] == str(chain)
        if resnum is not None:
            by_number &= atoms["resnum"] == int(resnum)
        if icode is not None:
            by_number &= atoms["icode"] == str(icode).strip()
        mask = by_number.copy()
        if resname is not None:
            mask &= atoms["resname"].str.upper() == str(resname).strip().upper()
        matched = atoms[mask]
        if wanted_atoms is not None:
            matched = matched[matched["atom"].isin(wanted_atoms)]
        residues = sorted({(r.chain, r.resname, r.resnum, r.icode) for r in matched.itertuples()})
        found_atoms = set(matched["atom"])
        missing = [a for a in (wanted_atoms or []) if a not in found_atoms]
        note = ""
        if not len(matched):
            state = "not_found"
            if resnum is not None and resname is not None and by_number.any():
                seen = sorted({r.resname for r in atoms[by_number].itertuples()})
                note = f"residue number matched but the residue name differs: {'|'.join(seen)}"
            else:
                note = "no matching residue in this file"
        elif missing:
            state = "partial"
        else:
            state = "matched"
        status.append({
            "structure_id": structure_id, "input_sha256": input_sha256, "source_group": label,
            "requested_resname": "" if resname is None else str(resname), "requested_atoms": "|".join(wanted_atoms or []),
            "requested_chain": "" if chain is None else str(chain), "requested_resnum": "" if resnum is None else int(resnum),
            "requested_icode": "" if icode is None else str(icode), "status": state,
            "n_residues_matched": len(residues), "n_atom_conformers": int(len(matched)),
            "matched_residues": "|".join(f"{c}:{n}{i}:{r}" for c, r, n, i in residues),
            "missing_atoms": "|".join(missing), "note": note,
        })
        if len(matched):
            piece = matched.copy()
            piece["group"] = label
            piece["charge"] = [
                (charge_override if charge_override else _charge_class(r.resname, r.atom)) if charge_override != "none" else ""
                for r in piece.itertuples()
            ]
            pieces.append(piece)
    sources_table = pd.concat(pieces, ignore_index=True) if pieces else atoms.iloc[0:0].assign(group="", charge="")
    return sources_table, pd.DataFrame(status, columns=STATUS_COLUMNS)


def _compatible(a: Any, b: Any) -> bool:
    return not (a and b and a != b)


def _contact_rows(atoms: pd.DataFrame, sources: pd.DataFrame, *, structure_id: str, input_sha256: str, policy: str,
                  neighbor: float, polar: float, close: float, ionic: float, include_same: bool) -> pd.DataFrame:
    rows: List[dict] = []
    partners_all = atoms.copy()
    partners_all["residue_key"] = _residue_key(partners_all)
    partners_all["charge"] = [_charge_class(r.resname, r.atom) for r in partners_all.itertuples()]
    sources = sources.copy()
    sources["residue_key"] = _residue_key(sources)
    for model_id, src_model in sources.groupby("model_id", sort=False):
        partners = partners_all[partners_all["model_id"] == model_id].reset_index(drop=True)
        pxyz = partners[["x", "y", "z"]].to_numpy(float)
        for s in src_model.itertuples():
            distances = np.linalg.norm(pxyz - np.array([s.x, s.y, s.z]), axis=1)
            for index in np.nonzero(distances <= neighbor)[0]:
                p = partners.iloc[index]
                same = p["residue_key"] == s.residue_key
                if same and not include_same:
                    continue
                if same and p["atom"] == s.atom:
                    continue                                              # the source atom itself or another conformer of it
                d = float(distances[index])
                is_polar = s.element in POLAR_ELEMENTS and p["element"] in POLAR_ELEMENTS and d <= polar
                opposite = {s.charge, p["charge"]} == {"negative", "positive"}
                separation = ""
                if s.kind == "protein" and p["kind"] == "protein" and s.chain == p["chain"]:
                    separation = abs(int(p["resnum"]) - int(s.resnum))
                rows.append({
                    "structure_id": structure_id, "input_sha256": input_sha256, "model_id": model_id, "altloc_policy": policy,
                    "source_group": s.group, "source_chain": s.chain, "source_resname": s.resname, "source_resnum": s.resnum,
                    "source_icode": s.icode, "source_hetero_flag": s.hetero_flag, "source_atom": s.atom,
                    "source_altloc": s.altloc, "source_occupancy": s.occupancy,
                    "partner_chain": p["chain"], "partner_resname": p["resname"], "partner_resnum": p["resnum"],
                    "partner_icode": p["icode"], "partner_hetero_flag": p["hetero_flag"], "partner_atom": p["atom"],
                    "partner_element": p["element"], "partner_altloc": p["altloc"], "partner_occupancy": p["occupancy"],
                    "partner_kind": p["kind"], "distance_A": round(d, 3), "polar_candidate": bool(is_polar),
                    "close_polar": bool(is_polar and d <= close), "ionic_candidate": bool(opposite and d <= ionic),
                    "same_residue": bool(same), "altloc_compatible": _compatible(s.altloc, p["altloc"]),
                    "seq_separation": separation,
                })
    frame = pd.DataFrame(rows, columns=CONTACT_COLUMNS)
    if len(frame):
        frame = frame.sort_values(
            ["model_id", "source_group", "source_chain", "source_resnum", "source_icode", "source_atom", "source_altloc",
             "distance_A", "partner_chain", "partner_resnum", "partner_atom", "partner_altloc"],
            kind="mergesort", ignore_index=True)
        frame["seq_separation"] = frame["seq_separation"].replace("", pd.NA).astype("Int64")
    return frame


def _class_labels(frame: pd.DataFrame) -> str:
    labels = []
    for flag in ("close_polar", "polar_candidate", "ionic_candidate"):
        if frame[flag].any():
            labels.append(flag)
    return "|".join(labels) if labels else "neighbor"


def _summary(contacts: pd.DataFrame, policy: str) -> pd.DataFrame:
    keys = ["structure_id", "model_id", "source_group", "source_chain", "source_resname", "source_resnum", "source_icode",
            "partner_chain", "partner_resname", "partner_resnum", "partner_icode", "partner_hetero_flag", "partner_kind"]
    rows = []
    for key, group in contacts.groupby(keys, sort=False, dropna=False):
        row = dict(zip(keys, key))
        row.update({
            "altloc_policy": policy, "min_distance_A": float(group["distance_A"].min()), "n_atom_pairs": int(len(group)),
            "source_atoms": "|".join(sorted(set(group["source_atom"]))),
            "partner_atoms": "|".join(sorted(set(group["partner_atom"]))),
            "partner_altlocs_present": "|".join(sorted({a for a in group["partner_altloc"] if a})),
            "classes": _class_labels(group),
        })
        rows.append(row)
    frame = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
    if len(frame):
        frame = frame.sort_values(["model_id", "source_group", "source_chain", "source_resnum", "min_distance_A",
                                   "partner_chain", "partner_resnum"], kind="mergesort", ignore_index=True)
    return frame


def _water_bridges(atoms: pd.DataFrame, sources: pd.DataFrame, *, structure_id: str, input_sha256: str, policy: str,
                   bridge: float, polar: float) -> pd.DataFrame:
    rows: List[dict] = []
    atoms = atoms.copy()
    atoms["residue_key"] = _residue_key(atoms)
    sources = sources.copy()
    sources["residue_key"] = _residue_key(sources)
    sources = sources[sources["element"].isin(POLAR_ELEMENTS)]
    for model_id, src_model in sources.groupby("model_id", sort=False):
        model_atoms = atoms[atoms["model_id"] == model_id]
        waters = model_atoms[(model_atoms["kind"] == "water") & (model_atoms["element"] == "O")].reset_index(drop=True)
        partners = model_atoms[(model_atoms["kind"] != "water") & model_atoms["element"].isin(POLAR_ELEMENTS)].reset_index(drop=True)
        if waters.empty or partners.empty:
            continue
        wxyz = waters[["x", "y", "z"]].to_numpy(float)
        pxyz = partners[["x", "y", "z"]].to_numpy(float)
        for s in src_model.itertuples():
            sxyz = np.array([s.x, s.y, s.z])
            for wi in np.nonzero(np.linalg.norm(wxyz - sxyz, axis=1) <= bridge)[0]:
                w = waters.iloc[wi]
                to_partner = np.linalg.norm(pxyz - wxyz[wi], axis=1)
                for pi in np.nonzero(to_partner <= bridge)[0]:
                    p = partners.iloc[pi]
                    if p["residue_key"] == s.residue_key:
                        continue
                    direct = float(np.linalg.norm(pxyz[pi] - sxyz))
                    rows.append({
                        "structure_id": structure_id, "input_sha256": input_sha256, "model_id": model_id, "altloc_policy": policy,
                        "source_group": s.group, "source_chain": s.chain, "source_resname": s.resname, "source_resnum": s.resnum,
                        "source_icode": s.icode, "source_atom": s.atom, "source_altloc": s.altloc, "source_occupancy": s.occupancy,
                        "water_chain": w["chain"], "water_resname": w["resname"], "water_resnum": w["resnum"],
                        "water_icode": w["icode"], "water_atom": w["atom"], "water_altloc": w["altloc"],
                        "water_occupancy": w["occupancy"],
                        "source_water_distance_A": round(float(np.linalg.norm(wxyz[wi] - sxyz)), 3),
                        "partner_chain": p["chain"], "partner_resname": p["resname"], "partner_resnum": p["resnum"],
                        "partner_icode": p["icode"], "partner_hetero_flag": p["hetero_flag"], "partner_atom": p["atom"],
                        "partner_element": p["element"], "partner_altloc": p["altloc"], "partner_occupancy": p["occupancy"],
                        "partner_kind": p["kind"], "water_partner_distance_A": round(float(to_partner[pi]), 3),
                        "source_partner_distance_A": round(direct, 3), "direct_polar_candidate": bool(direct <= polar),
                        "altloc_compatible": bool(_compatible(s.altloc, w["altloc"]) and _compatible(w["altloc"], p["altloc"])
                                                  and _compatible(s.altloc, p["altloc"])),
                    })
    frame = pd.DataFrame(rows, columns=BRIDGE_COLUMNS)
    if len(frame):
        frame = frame.sort_values(["model_id", "source_group", "source_chain", "source_resnum", "source_atom", "source_altloc",
                                   "source_water_distance_A", "water_resnum", "partner_resnum", "partner_atom"],
                                  kind="mergesort", ignore_index=True)
    return frame


def _repository_context() -> Dict[str, Any]:
    root = Path(__file__).resolve().parents[1]
    version = None
    try:
        match = re.search(r'^version\s*=\s*"([^"]+)"', (root / "pyproject.toml").read_text(), re.MULTILINE)
        version = match.group(1) if match else None
    except OSError:
        pass

    def git(*args):
        try:
            done = subprocess.run(["git", *args], capture_output=True, text=True, cwd=root, timeout=10)
            return done.stdout.strip() if done.returncode == 0 else None
        except Exception:
            return None

    status = git("status", "--porcelain")
    return {
        "package_version": version, "git_commit": git("rev-parse", "HEAD"),
        "git_working_tree_clean": None if status is None else status == "",
        "module_sha256": _sha256(Path(__file__).resolve()),
    }


# ----------------------------------------------------------------------------- public API
def site_contacts(
    structure_path,
    sources: Sequence[Mapping[str, Any]],
    *,
    neighbor_cutoff: float = 4.5,
    polar_cutoff: float = 3.5,
    close_polar_cutoff: float = 3.2,
    ionic_cutoff: float = 4.0,
    water_bridge_cutoff: float = 3.5,
    altloc_policy: str = "all",
    models: str = "all",
    include_same_residue: bool = False,
    compare_highest_occupancy: bool = True,
) -> Dict[str, Any]:
    """Atom-level neighbors, residue summary and single-water bridges of the selected source atoms of one structure.

    Returns a dict: ``site_atom_contacts``, ``site_residue_summary``, ``site_water_bridges`` and ``site_status`` (DataFrames; columns
    in ``CONTACT_COLUMNS``, ``SUMMARY_COLUMNS`` (+ ``SUMMARY_COMPARISON_COLUMNS``), ``BRIDGE_COLUMNS``, ``STATUS_COLUMNS``) and
    ``manifest`` (a JSON-ready dict: input SHA-256, tool version, parameters, counts). A source that matches nothing is a ``site_status``
    row (``not_found``), never an error. See the module docstring for the selectors, definitions and altloc policy.
    """
    if altloc_policy not in ALTLOC_POLICIES:
        raise ValueError(f"altloc_policy must be one of {ALTLOC_POLICIES}")
    if models not in ("all", "first"):
        raise ValueError("models must be 'all' or 'first'")
    if not 0 < close_polar_cutoff <= polar_cutoff <= neighbor_cutoff or not 0 < ionic_cutoff <= neighbor_cutoff:
        raise ValueError("cutoffs must satisfy 0 < close_polar <= polar <= neighbor and 0 < ionic <= neighbor")
    if water_bridge_cutoff <= 0:
        raise ValueError("water_bridge_cutoff must be positive")
    if not sources:
        raise ValueError("at least one source selector is required")

    path = Path(structure_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Structure file not found: {path}")
    digest = _sha256(path)
    structure_id = _structure_id(path)
    atoms, counts = _load_atoms(_opened(path), models)
    source_table, status = _match_sources(atoms, sources, structure_id, digest)

    def run(policy: str):
        pool = atoms if policy == "all" else atoms[atoms["selected"]]
        chosen = source_table if policy == "all" else source_table[source_table["selected"]]
        contacts = _contact_rows(pool, chosen, structure_id=structure_id, input_sha256=digest, policy=policy,
                                 neighbor=neighbor_cutoff, polar=polar_cutoff, close=close_polar_cutoff,
                                 ionic=ionic_cutoff, include_same=include_same_residue)
        return pool, chosen, contacts

    pool, chosen, contacts = run(altloc_policy)
    summary = _summary(contacts, altloc_policy)
    if compare_highest_occupancy:
        if altloc_policy == "highest_occupancy":
            reference = contacts
        else:
            _, _, reference = run("highest_occupancy")
        ref = _summary(reference, "highest_occupancy")
        keys = [c for c in SUMMARY_COLUMNS if c.startswith(("structure_id", "model_id", "source_", "partner_")) and
                c not in ("source_atoms", "partner_atoms", "partner_altlocs_present")]
        ref = ref[keys + ["min_distance_A"]].rename(columns={"min_distance_A": "min_distance_highest_occupancy_A"})
        if len(summary) and len(ref):
            summary = summary.merge(ref, on=keys, how="left")
        else:
            summary["min_distance_highest_occupancy_A"] = np.nan
        summary["in_highest_occupancy"] = summary["min_distance_highest_occupancy_A"].notna()
    else:
        summary["in_highest_occupancy"] = pd.NA
        summary["min_distance_highest_occupancy_A"] = np.nan
    summary = summary.reindex(columns=SUMMARY_COLUMNS + SUMMARY_COMPARISON_COLUMNS)
    bridges = _water_bridges(pool, chosen, structure_id=structure_id, input_sha256=digest, policy=altloc_policy,
                             bridge=water_bridge_cutoff, polar=polar_cutoff)

    import Bio

    manifest = {
        "tool": "modules/site_contacts.py", "site_contacts_version": SITE_CONTACTS_VERSION,
        **_repository_context(),
        "python": platform.python_version(), "biopython": Bio.__version__, "numpy": np.__version__, "pandas": pd.__version__,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat(),
        "input": {"path": str(path), "size_bytes": path.stat().st_size, "sha256": digest, "structure_id": structure_id},
        "parameters": {
            "neighbor_cutoff_A": neighbor_cutoff, "polar_cutoff_A": polar_cutoff, "close_polar_cutoff_A": close_polar_cutoff,
            "ionic_cutoff_A": ionic_cutoff, "water_bridge_cutoff_A": water_bridge_cutoff, "altloc_policy": altloc_policy,
            "models": models, "include_same_residue": include_same_residue, "heavy_atoms_only": True,
            "compare_highest_occupancy": compare_highest_occupancy,
        },
        "sources": [dict(spec) for spec in sources],
        "counts": {**counts, "n_heavy_atom_conformers": int(len(atoms)), "n_source_atom_conformers": int(len(source_table)),
                   "n_atom_contact_rows": int(len(contacts)), "n_residue_summary_rows": int(len(summary)),
                   "n_water_bridge_rows": int(len(bridges))},
        "symmetry": SYMMETRY_NOTE,
        "caveats": [
            "polar_candidate is a heavy-atom N/O pair within the cutoff, not a hydrogen bond: no hydrogen is placed and no angle is judged",
            "ionic_candidate uses fixed charge classes (Asp, Glu, HEM propionate negative; Lys, Arg positive); histidine protonation is not assumed",
            "altloc labels are only guaranteed consistent within a residue; see the altloc_compatible column",
            "a contact table describes a crystal structure; it does not choose a quantum region, protonation states, or a model of another species",
        ],
    }
    return {"site_atom_contacts": contacts, "site_residue_summary": summary, "site_water_bridges": bridges,
            "site_status": status, "manifest": manifest}


def write_site_contacts(result: Mapping[str, Any], output_dir, *, suffix: str = "") -> List[Path]:
    """Write the tables of :func:`site_contacts` as CSV plus ``site_manifest{suffix}.json`` (with the SHA-256 of each CSV).

    ``output_dir`` is created if missing; nothing else is written. Returns the paths written, manifest last.
    """
    folder = Path(output_dir)
    folder.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    files: Dict[str, Any] = {}
    for name in TABLE_NAMES:
        target = folder / f"{name}{suffix}.csv"
        result[name].to_csv(target, index=False, lineterminator="\n")
        written.append(target)
        files[target.name] = {"sha256": _sha256(target), "rows": int(len(result[name]))}
    manifest = dict(result["manifest"])
    manifest["output_files"] = files
    manifest_path = folder / f"site_manifest{suffix}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=False, default=str) + "\n")
    written.append(manifest_path)
    return written


__all__ = ["site_contacts", "write_site_contacts", "CONTACT_COLUMNS", "SUMMARY_COLUMNS", "SUMMARY_COMPARISON_COLUMNS",
           "BRIDGE_COLUMNS", "STATUS_COLUMNS", "ALTLOC_POLICIES", "SITE_CONTACTS_VERSION"]
