"""Protonation and hydrogen-bond geometry from optional external tools.

Crystal structures rarely show protons, so whether an interaction is a hydrogen bond, and
whether a carboxylate or histidine can carry a proton, has to be inferred. This module
wraps two optional tools and adds geometry helpers that need neither:

* :func:`run_propka` -- predicted pKa of every titratable group (``pip install propka``);
* :func:`add_hydrogens` -- hydrogens added with the ``reduce`` executable (flips Asn, Gln,
  and His and optimizes hydroxyl rotamers);
* :func:`directional_hbonds` -- donors whose hydrogen points at an acceptor
  (H...A <= 2.5 Å, D-H...A >= 120°, D...A <= 3.5 Å), with an option to test a histidine as
  doubly protonated by placing a virtual hydrogen on each ring nitrogen that lacks one;
* :func:`backbone_amide_hydrogen`, :func:`ring_nh_hydrogen`, :func:`amide_nh2_hydrogens` and
  :func:`hydrogen_bond_geometry` -- hydrogens placed from heavy atoms for the planar donors
  (backbone amide, Trp NE1, Arg, Asn/Gln amide, a His ring N as if protonated);
* :func:`carboxylate_distances` -- closest carboxylate oxygen to an acceptor, since a
  carboxylate cannot donate to a carboxylate unless one is protonated and a short O...O
  contact (about 2.7 Å or less) is the geometry of a protonated pair.

PROPKA is empirical and built on heavy-atom structures. It does not know the redox,
ligation, or charge state of metal centers, so the histidines that ligate metals return
large negative values and the result depends on whether cofactors are in the input; use
it as an environment descriptor and check the sensitivity (with and without cofactors).
"""

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

TITRATABLE = ("ASP", "GLU", "HIS", "LYS", "ARG", "TYR", "CYS")
CARBOXYL_OXYGENS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
HIS_RING = ("CG", "ND1", "CD2", "CE1", "NE2")


def run_propka(pdb_path: str, *, residue_types: Sequence[str] = TITRATABLE) -> pd.DataFrame:
    """Predicted pKa for each titratable group of a PDB file (PROPKA 3, average conformation).

    Columns: ``chain, residue_number, insertion_code, residue_name, pka, model_pka, shift,
    buried`` (burial fraction). Raises ``ImportError`` with an install hint if PROPKA is absent.
    """
    try:
        import propka.run as propka_run
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError("PROPKA is not installed; run `pip install propka`") from exc
    molecule = propka_run.single(str(pdb_path), optargs=["--quiet"], stream=None, write_pka=False)
    rows = []
    for group in molecule.conformations["AVR"].groups:
        if group.residue_type in residue_types:
            atom = group.atom
            rows.append({"chain": atom.chain_id, "residue_number": int(atom.res_num), "insertion_code": (atom.icode or "").strip(),
                         "residue_name": group.residue_type, "pka": float(group.pka_value), "model_pka": float(group.model_pka),
                         "shift": float(group.pka_value) - float(group.model_pka), "buried": float(group.buried)})
    return pd.DataFrame(rows, columns=["chain", "residue_number", "insertion_code", "residue_name", "pka", "model_pka", "shift", "buried"])


def add_hydrogens(pdb_path: str, output_path: Optional[str] = None, *, flip: bool = True, timeout: int = 600) -> str:
    """Run ``reduce`` on a PDB file and return the path of the hydrogenated file."""
    exe = shutil.which("reduce")
    if exe is None:
        raise FileNotFoundError("the `reduce` executable was not found on PATH")
    out = Path(output_path) if output_path else Path(tempfile.mkdtemp()) / (Path(pdb_path).stem + "_H.pdb")
    command = [exe, "-Quiet"] + (["-FLIP"] if flip else []) + [str(pdb_path)]
    with open(out, "w") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.DEVNULL, check=False, timeout=timeout)
    return str(out)


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    v1, v2 = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0))))


def virtual_histidine_hydrogens(ring_atoms: Dict[str, np.ndarray], existing_h: Sequence[np.ndarray]) -> List[tuple]:
    """Donor/hydrogen coordinate pairs for ring nitrogens of a histidine that lack a hydrogen.

    The hydrogen is placed 1.0 Å from the nitrogen along the line from the ring centre,
    which tests the doubly protonated form without changing the original model.
    """
    ring = [ring_atoms[n] for n in HIS_RING if n in ring_atoms]
    if len(ring) < 4:
        return []
    centre = np.mean(ring, axis=0)
    pairs = []
    for name in ("ND1", "NE2"):
        if name not in ring_atoms:
            continue
        n = ring_atoms[name]
        if any(np.linalg.norm(n - h) < 1.3 for h in existing_h):
            continue
        u = (n - centre) / np.linalg.norm(n - centre)
        pairs.append((n, n + u))
    return pairs


def directional_hbonds(donor_hydrogen_pairs: Sequence[tuple], acceptors: np.ndarray, *, h_acceptor_max: float = 2.5,
                       angle_min: float = 120.0, donor_acceptor_max: float = 3.5) -> bool:
    """Whether any donor hydrogen points at any acceptor.

    ``donor_hydrogen_pairs`` is a sequence of ``(donor_xyz, hydrogen_xyz)``; ``acceptors``
    is an ``(n, 3)`` array. True when some pair has H...A <= ``h_acceptor_max``, D-H...A >=
    ``angle_min`` degrees, and D...A <= ``donor_acceptor_max``.
    """
    acceptors = np.asarray(acceptors, float).reshape(-1, 3)
    for donor, hydrogen in donor_hydrogen_pairs:
        donor, hydrogen = np.asarray(donor, float), np.asarray(hydrogen, float)
        for acceptor in acceptors:
            if (np.linalg.norm(hydrogen - acceptor) <= h_acceptor_max and np.linalg.norm(donor - acceptor) <= donor_acceptor_max
                    and _angle(donor, hydrogen, acceptor) >= angle_min):
                return True
    return False


# Hydrogens placed from heavy atoms alone, for planar donors. Hydroxyl and Lys NZ hydrogens rotate and are not placed.
NH_LENGTH = 1.01
SULFUR_H_MAX = 2.9           # H...S limit of an H-bond-like contact to a sulfur
SULFUR_ANGLE_MIN = 120.0     # N-H...S angle limit
# (residue, atom) -> the two ring neighbors; the H lies on the exocyclic bisector (a His ring N is treated as if protonated).
RING_DONORS = {("TRP", "NE1"): ("CD1", "CE2"), ("HIS", "ND1"): ("CG", "CE1"), ("HIS", "NE2"): ("CD2", "CE1"), ("ARG", "NE"): ("CD", "CZ")}
# (residue, atom) -> (bonded C, reference atom defining the plane); two H at 120 degrees from the C-N bond.
AMIDE_DONORS = {("ASN", "ND2"): ("CG", "OD1"), ("GLN", "NE2"): ("CD", "OE1"), ("ARG", "NH1"): ("CZ", "NE"), ("ARG", "NH2"): ("CZ", "NE")}


def backbone_amide_hydrogen(n: np.ndarray, c_prev: np.ndarray, ca: np.ndarray) -> np.ndarray:
    """Amide H on the bisector of the C(i-1)-N and CA-N bonds, ``NH_LENGTH`` from N."""
    u = (n - c_prev) / np.linalg.norm(n - c_prev)
    v = (n - ca) / np.linalg.norm(n - ca)
    return n + NH_LENGTH * (u + v) / np.linalg.norm(u + v)


def ring_nh_hydrogen(n: np.ndarray, neighbor_a: np.ndarray, neighbor_b: np.ndarray) -> np.ndarray:
    """H on the exocyclic bisector of a ring N (Trp NE1, a His ring N as if protonated, Arg NE)."""
    u, v = neighbor_a - n, neighbor_b - n
    direction = -(u / np.linalg.norm(u) + v / np.linalg.norm(v))
    return n + NH_LENGTH * direction / np.linalg.norm(direction)


def amide_nh2_hydrogens(n: np.ndarray, c: np.ndarray, reference: np.ndarray) -> List[np.ndarray]:
    """The two H of an NH2 group, in the plane of ``reference``, 120 degrees from the C-N bond."""
    u = (n - c) / np.linalg.norm(n - c)
    w = reference - c
    w = w - np.dot(w, u) * u
    w /= np.linalg.norm(w)
    return [n + NH_LENGTH * (np.cos(np.radians(60)) * u + sign * np.sin(np.radians(60)) * w) for sign in (1.0, -1.0)]


def hydrogen_bond_geometry(h: np.ndarray, donor: np.ndarray, acceptor: np.ndarray) -> tuple:
    """(H...acceptor distance, donor-H...acceptor angle in degrees)."""
    return float(np.linalg.norm(acceptor - h)), _angle(donor, h, acceptor)


def carboxylate_distances(atoms: pd.DataFrame, acceptors: pd.DataFrame) -> pd.DataFrame:
    """Closest Asp/Glu carboxylate oxygen to any acceptor atom, per structure and position.

    ``atoms`` and ``acceptors`` follow the tables of :mod:`modules.conserved_interactions`.
    """
    rows = []
    for sid, frame in atoms.groupby("structure_id", sort=False):
        acc = acceptors[acceptors["structure_id"] == sid][["x", "y", "z"]].to_numpy(float)
        if len(acc) == 0:
            continue
        for position, res in frame.groupby("position", sort=False):
            name = res["resname"].iloc[0].upper()
            if name not in CARBOXYL_OXYGENS:
                continue
            ox = res[res["atom"].isin(CARBOXYL_OXYGENS[name])][["x", "y", "z"]].to_numpy(float)
            if len(ox):
                d = np.sqrt(((ox[:, None, :] - acc[None, :, :]) ** 2).sum(-1))
                rows.append({"structure_id": sid, "position": position, "resname": name, "oo_distance": float(d.min())})
    return pd.DataFrame(rows, columns=["structure_id", "position", "resname", "oo_distance"])
