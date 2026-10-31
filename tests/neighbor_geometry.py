"""Exact-geometry two-residue chains for the neighbor-contact tests.

Residue 1 and residue 2 have real backbone bond lengths. Residue 2 is turned about the peptide (C1-N2) bond so that one
of its atoms is at a chosen distance from a point near residue 1, and side-chain atoms are placed at an exact bond length
from their parent and an exact distance from a target atom, so a window can be probed 0.01 A either side.
"""

import numpy as np
import pandas as pd

RES1 = {"N": (0.0, 0.0, 0.0), "CA": (1.46, 0.0, 0.0), "C": (2.0, 1.4, 0.0), "O": (1.3, 2.4, 0.0), "CB": (2.0, -0.8, 1.2)}
RES2 = {"N": (3.2, 1.9, 0.2), "CA": (4.4, 1.2, 0.6), "C": (5.6, 2.1, 0.3), "O": (5.7, 3.3, 0.1), "CB": (4.4, 0.4, 1.9)}
COLUMNS = ["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"]


def rotate(points, origin, axis, angle):
    """Rodrigues rotation of one point (or an array of points) about the line through ``origin`` along ``axis``."""
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    shifted = np.asarray(points, float) - origin
    c, s = np.cos(angle), np.sin(angle)
    return origin + shifted * c + np.cross(axis, shifted) * s + np.outer(shifted @ axis, axis).reshape(shifted.shape) * (1 - c)


def turned_residue2(atom, reach, distance):
    """Residue 2 turned about the C1-N2 bond so that its ``atom`` is as close as possible to ``distance`` from the point ``reach``."""
    origin, axis = np.array(RES1["C"]), np.array(RES2["N"]) - np.array(RES1["C"])
    best = min(np.linspace(0, 2 * np.pi, 3600, endpoint=False),
               key=lambda a: abs(np.linalg.norm(rotate(RES2[atom], origin, axis, a) - np.asarray(reach, float)) - distance))
    return {name: rotate(xyz, origin, axis, best) for name, xyz in RES2.items()}


def point(parent, target, bond, distance):
    """A point ``bond`` from ``parent`` and ``distance`` from ``target`` (one of the circle of solutions)."""
    parent, target = np.asarray(parent, float), np.asarray(target, float)
    length = np.linalg.norm(target - parent)
    axis = (target - parent) / length
    cosine = (bond ** 2 + length ** 2 - distance ** 2) / (2 * bond * length)
    assert -1 <= cosine <= 1, f"no point {bond} from the parent and {distance} from the target (they are {length:.2f} apart)"
    side = np.cross(axis, [0.0, 0.0, 1.0])
    side /= np.linalg.norm(side)
    return parent + bond * (cosine * axis + np.sqrt(1 - cosine ** 2) * side)


def frame(resname1, resname2, side1, backbone2, side2):
    """The atoms table of the chain; ``side1`` and ``side2`` map side-chain atom names to coordinates."""
    rows = [("s", 1, 1, resname1, name, name[0], *map(float, xyz)) for name, xyz in {**RES1, **side1}.items()]
    rows += [("s", 2, 2, resname2, name, name[0], *map(float, xyz)) for name, xyz in {**backbone2, **side2}.items()]
    table = pd.DataFrame(rows, columns=COLUMNS)
    xyz = table[["x", "y", "z"]].to_numpy(float)
    distance = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    across = table["position"].to_numpy()[:, None] != table["position"].to_numpy()[None]
    close = np.argwhere(across & (distance < 1.9))
    assert {tuple(sorted((table.atom.iloc[a], table.atom.iloc[b]))) for a, b in close} <= {("C", "N")}, "an unintended cross-residue bond"
    return table


def hydroxyl_to_carbonyl(distance):
    """Ser 2 hydroxyl O (a donor whose hydrogen is not placed) ``distance`` from the carbonyl O of residue 1 (Ala)."""
    backbone = turned_residue2("CB", RES1["O"], 3.7)
    return frame("ALA", "SER", {}, backbone, {"OG": point(backbone["CB"], RES1["O"], 1.43, distance)})


def leucine_valine(distance):
    """Leu 1 CG ``distance`` from Val 2 CG1: two side-chain carbons beyond C-beta."""
    cg = np.array(RES1["CB"]) + 1.53 * np.array([0.0, -0.3, 1.0]) / np.linalg.norm([0.0, -0.3, 1.0])
    backbone = turned_residue2("CB", cg, 4.0)
    return frame("LEU", "VAL", {"CG": cg}, backbone, {"CG1": point(backbone["CB"], cg, 1.53, distance)})


def phenylalanine_serine(distance):
    """Phe 1 ring atom CG ``distance`` from Ser 2 OG: a ring atom against a side-chain O beyond C-beta."""
    cg = np.array(RES1["CB"]) + 1.51 * np.array([0.0, -0.3, 1.0]) / np.linalg.norm([0.0, -0.3, 1.0])
    backbone = turned_residue2("CB", cg, 4.0)
    return frame("PHE", "SER", {"CG": cg}, backbone, {"OG": point(backbone["CB"], cg, 1.43, distance)})
