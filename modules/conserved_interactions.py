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
  missing contacts;
* :func:`backbone_only_audit` -- flags backbone-only models, which cannot show a
  side-chain contact.

Input tables
------------
``atoms``: ``structure_id, position, seq_index, resname, atom, element, x, y, z`` (one
row per heavy atom; ``seq_index`` is the residue's order in its own chain, used to
exclude sequence neighbors). ``features``: ``structure_id, feature, atom, element, kind,
x, y, z`` with ``kind`` in ``{"metal", "oxygen"}``. ``waters``: ``structure_id, x, y, z``.
``structures``: ``structure_id, group, sequence_id``.

Geometry is judged from heavy atoms, so "polar" means a hydrogen-bond-capable N/O pair, not a hydrogen bond. Only
the opt-in flags (``polar_s_geom``, the neighbor ``polar`` class) place a hydrogen from the heavy atoms of the backbone or a planar group,
a model and not an observation (see :mod:`modules.protonation`).
"""

from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import inspect

import numpy as np
import pandas as pd

from . import coordination_table, protonation
from .coordination_api import METAL_ION_RESNAMES

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
SULFUR_FLAGS = ("polar_s", "polar_s_bb", "polar_s_geom")
# Opt-in sulfur-acceptor typing (``residue_contacts(sulfur_acceptor=True)``): acceptor atom -> largest donor...S heavy-atom
# distance. Windows from docs/scope_checks/03_sulfur_acceptor_survey.md. A disulfide SG is never an acceptor.
SULFUR_ACCEPTORS = {("CYS", "SG"): 3.7, ("MET", "SD"): 3.6}
NEIGHBOR_CHEMISTRY = ("polar", "ionic")                           # residue_contacts(neighbor_contacts=True) default
NEIGHBOR_CHEMISTRY_ALL = ("polar", "ionic", "hydrophobic", "aromatic")
NEIGHBOR_H_MAX = 2.5            # neighbor polar pairs: H...acceptor limit when the donor's hydrogen can be placed (N-H...A >= 120 degrees too)
NEIGHBOR_MIN_BONDS = 4          # residue_contacts(neighbor_contacts=True): bonded 1-2, 1-3 and 1-4 atom pairs are not contacts
BOND_DISTANCE = 1.9             # heavy atoms closer than this are bonded when counting bonds
POLAR_DISTANCE = 3.5            # residue_contacts: N/O pair for the polar flags
CHEMISTRY_DISTANCE = 4.0        # residue_contacts: hydrophobic, aromatic and ionic windows
SULFUR_SIDE_CHAIN_MAX = 3.6     # side-chain donors: the share pointing an N-H at the sulfur falls sharply beyond this
# Either His ring N may carry the H (which one depends on the tautomer, unknown from heavy atoms), so both are capable
# donors, except an N that a metal coordinates, which is deprotonated; ``metals`` supplies that exclusion.
SULFUR_DONORS = frozenset({("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("ASN", "ND2"), ("GLN", "NE2"), ("ARG", "NE"),
                           ("ARG", "NH1"), ("ARG", "NH2"), ("LYS", "NZ"), ("TRP", "NE1"), ("HIS", "ND1"), ("HIS", "NE2")})
DISULFIDE_DISTANCE = 2.3
FEATURE_CLASSES = ("coord", "polar", "ionic", "contact")
TIER_COLUMNS = coordination_table.TIER_COLUMNS
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
    map the residue number is the position. ``seq_index`` follows the residue numbering, so
    residues that are unobserved or skipped (hetero residues such as modified amino acids) still
    count as steps in the sequence: two residues are neighbors only when their numbers are
    consecutive (an insertion code or a number that does not increase counts as one step).
    Hydrogens are dropped, and an atom with alternate locations keeps its first conformer (``A`` when
    there is one).
    """
    from .io_utils import unpack_pdb_file

    structure, _ = unpack_pdb_file(str(path))
    model = list(structure)[model_index]
    sid = structure_id or Path(path).stem
    rows = []
    order = 0
    previous = None
    for residue in model[chain_id]:
        if residue.id[0] != " ":
            continue
        number = int(residue.id[1])
        order += 1 if previous is None or number <= previous else number - previous
        previous = number
        key = (number, residue.id[2].strip())
        position = key[0] if position_map is None else position_map.get(key)
        if position is None:
            continue
        # get_unpacked_list() yields each conformer of a disordered atom; iterating the residue
        # would give one DisorderedAtom whose altloc is that of the highest-occupancy conformer,
        # so the atom would vanish whenever B outweighs A.
        # Each atom keeps its first conformer (blank, else the lowest letter, which is A whenever there is one).
        chosen: Dict[str, object] = {}
        for atom in residue.get_unpacked_list():
            if atom.element.upper() in ("H", "D"):
                continue
            kept = chosen.get(atom.get_id())
            if kept is None or atom.get_altloc() < kept.get_altloc():
                chosen[atom.get_id()] = atom
        for atom in chosen.values():
            x, y, z = (float(v) for v in atom.coord)
            rows.append((sid, position, order, residue.resname.strip().upper(), atom.get_id(), atom.element.upper(), x, y, z))
    return pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])


def load_metals(path: str, *, structure_id: Optional[str] = None, model_index: int = 0) -> pd.DataFrame:
    """Metal ions of one model as the ``metals`` table (``structure_id, element, x, y, z``), all chains."""
    from .io_utils import unpack_pdb_file

    structure, _ = unpack_pdb_file(str(path))
    sid = structure_id or Path(path).stem
    rows = []
    for residue in list(structure)[model_index].get_residues():
        if residue.id[0] == " " or residue.id[0] == "W":
            continue
        for atom in residue:
            if atom.element.upper() in METAL_ION_RESNAMES:
                rows.append((sid, atom.element.upper(), *(float(v) for v in atom.coord)))
    return pd.DataFrame(rows, columns=["structure_id", "element", "x", "y", "z"])


def _metal_coordinated(frame: pd.DataFrame, xyz: np.ndarray, metals: pd.DataFrame) -> np.ndarray:
    """Atoms of ``frame`` that a metal coordinates: ``bond`` tier, or the global 2.8 A when the pair is not in the table."""
    coordinated = np.zeros(len(frame), bool)
    if metals.empty:
        return coordinated
    element = metals["element"].str.upper().to_numpy()
    i, j, d = _close_pairs(xyz, metals[["x", "y", "z"]].to_numpy(float), 2.8)
    res, name = frame["resname"].str.upper().to_numpy(), frame["atom"].to_numpy()
    for atom, metal, distance in zip(i, j, d):
        if coordination_table.coordination_tier(element[metal], res[atom], name[atom], float(distance)).tier in ("bond", "default"):
            coordinated[atom] = True
    return coordinated


def placed_hydrogens(frame: pd.DataFrame, xyz: np.ndarray) -> np.ndarray:
    """``(n, 2, 3)`` hydrogens placed on the planar donors of ``frame`` (NaN where there is none): see :mod:`modules.protonation`."""
    out = np.full((len(frame), 2, 3), np.nan)
    pos, seq = frame["position"].to_numpy(), frame["seq_index"].to_numpy()
    res, name = frame["resname"].str.upper().to_numpy(), frame["atom"].to_numpy()
    index = {(p, n): k for k, (p, n) in enumerate(zip(pos, name))}
    position_of = dict(zip(seq, pos))
    for k in range(len(frame)):
        key = (res[k], name[k])
        if name[k] == "N" and res[k] != "PRO":
            previous = position_of.get(seq[k] - 1)
            c, ca = index.get((previous, "C")), index.get((pos[k], "CA"))
            if c is not None and ca is not None and np.linalg.norm(xyz[c] - xyz[k]) < 1.6:      # a peptide bond joins them
                out[k, 0] = protonation.backbone_amide_hydrogen(xyz[k], xyz[c], xyz[ca])
        elif key in protonation.RING_DONORS:
            a, b = (index.get((pos[k], n)) for n in protonation.RING_DONORS[key])
            if a is not None and b is not None:
                out[k, 0] = protonation.ring_nh_hydrogen(xyz[k], xyz[a], xyz[b])
        elif key in protonation.AMIDE_DONORS:
            c, ref = (index.get((pos[k], n)) for n in protonation.AMIDE_DONORS[key])
            if c is not None and ref is not None:
                out[k] = protonation.amide_nh2_hydrogens(xyz[k], xyz[c], xyz[ref])
    return out


def _bond_steps(xyz: np.ndarray, members: np.ndarray, start: int) -> np.ndarray:
    """Bonds from atom ``start`` to each atom of ``members`` (indices into ``xyz``) within one residue; inf when unreachable."""
    local = xyz[members]
    adjacent = np.linalg.norm(local[:, None, :] - local[None, :, :], axis=-1) <= BOND_DISTANCE
    steps = np.full(len(members), np.inf)
    first = int(np.nonzero(members == start)[0][0])
    steps[first] = 0
    frontier = [first]
    while frontier:
        following = []
        for atom in frontier:
            for other in np.nonzero(adjacent[atom] & np.isinf(steps))[0]:
                steps[other] = steps[atom] + 1
                following.append(int(other))
        frontier = following
    return steps


def _neighbor_flags(elem, res, name, bb, a, b, d, cutoff, xyz, hydrogens, excluded) -> Dict[str, np.ndarray]:
    """The neighbor-contact chemistry of atom pairs ``(a, b)`` at distances ``d``, one boolean array per class.

    ``polar``: a hydrogen-bond pair with at least one side-chain atom: a donor (backbone N except Pro, or a side-chain atom
    of ``SULFUR_DONORS``) and an acceptor (any O, or a His ring N) within the polar window; when the donor's hydrogen can be placed
    it must point at the acceptor (H...A <= ``NEIGHBOR_H_MAX``, N-H...A >= 120 degrees), otherwise (hydroxyl, Lys NZ) the distance
    decides. A donor or acceptor that a metal coordinates (``excluded``) does not count. ``ionic``: a charged side-chain pair
    within the chemistry window. ``hydrophobic``: a C/S pair, both atoms side-chain atoms beyond C-beta, within the chemistry
    window. ``aromatic``: a ring atom against a side-chain atom beyond C-beta, within the chemistry window. Backbone and C-beta
    proximity between neighbors is not counted as hydrophobic or aromatic chemistry (it is set by the backbone).
    """
    beyond = ~bb & (name != "CB")
    ring = np.array([n in RING_ATOMS.get(r, ()) for n, r in zip(name, res)])
    positive = np.array([(r, n) in POSITIVE_ATOMS for n, r in zip(name, res)])
    negative = np.array([(r, n) in NEGATIVE_ATOMS for n, r in zip(name, res)])
    donor = (((name == "N") & (res != "PRO")) | np.array([(r, n) in SULFUR_DONORS for r, n in zip(res, name)])) & ~excluded
    acceptor = ((elem == "O") | ((res == "HIS") & np.isin(name, ["ND1", "NE2"]))) & ~excluded
    has_hydrogen = ~np.isnan(hydrogens[:, 0, 0])
    near_polar, near_chemistry = d <= min(POLAR_DISTANCE, cutoff), d <= min(CHEMISTRY_DISTANCE, cutoff)

    def hydrogen_bond(x, y):
        return donor[x] & acceptor[y] & near_polar & np.where(has_hydrogen[x], _points_at(hydrogens, xyz, x, y, acceptor, NEIGHBOR_H_MAX), True)
    return {
        "polar": (hydrogen_bond(a, b) | hydrogen_bond(b, a)) & (~bb[a] | ~bb[b]),
        "ionic": ((positive[a] & negative[b]) | (negative[a] & positive[b])) & near_chemistry,
        "hydrophobic": np.isin(elem[a], ["C", "S"]) & np.isin(elem[b], ["C", "S"]) & beyond[a] & beyond[b] & near_chemistry,
        "aromatic": ((ring[a] & beyond[b]) | (ring[b] & beyond[a])) & near_chemistry,
    }


def _neighbor_atom_pairs(xyz, pos, seq, name, bb, qualifies, i, j):
    """Atom pairs ``(i, j)`` from ``qualifies`` between residues one apart in sequence that are at least ``NEIGHBOR_MIN_BONDS`` bonds apart.

    The bond count runs through the peptide bond: steps from the first residue's atom to its C, one for C-N, and steps from the
    second residue's N to its atom. Pairs are dropped when the peptide bond is absent (a gap or a missing atom) or an atom is not
    connected to its backbone.
    """
    members: Dict[int, np.ndarray] = {}
    cache: Dict[tuple, np.ndarray] = {}

    def residue(s):
        if s not in members:
            members[s] = np.nonzero(seq == s)[0]
        return members[s]

    def atom(s, label):
        hit = [k for k in residue(s) if name[k] == label]
        return hit[0] if hit else None

    def steps(s, start):
        if (s, start) not in cache:
            cache[(s, start)] = _bond_steps(xyz, residue(s), start)
        return cache[(s, start)]

    kept_i, kept_j = [], []
    for a, b in zip(i[qualifies], j[qualifies]):
        low, high = (a, b) if seq[a] < seq[b] else (b, a)
        carbon, nitrogen = atom(seq[low], "C"), atom(seq[high], "N")
        if carbon is None or nitrogen is None or np.linalg.norm(xyz[carbon] - xyz[nitrogen]) >= 1.6:
            continue
        before, after = steps(seq[low], carbon), steps(seq[high], nitrogen)
        near, far = before[np.nonzero(residue(seq[low]) == low)[0][0]], after[np.nonzero(residue(seq[high]) == high)[0][0]]
        if np.isfinite(near) and np.isfinite(far) and near + 1 + far >= NEIGHBOR_MIN_BONDS:
            kept_i.append(a)
            kept_j.append(b)
    return np.array(kept_i, dtype=int), np.array(kept_j, dtype=int)


def _points_at(hydrogens: np.ndarray, xyz: np.ndarray, donor: np.ndarray, acceptor: np.ndarray, is_acceptor: np.ndarray,
               h_max: float = protonation.SULFUR_H_MAX) -> np.ndarray:
    """Per pair: a hydrogen on atom ``donor[i]`` is within ``h_max`` of ``acceptor[i]`` at an angle of at least ``SULFUR_ANGLE_MIN``."""
    result = np.zeros(len(donor), bool)
    live = np.nonzero(is_acceptor[acceptor] & ~np.isnan(hydrogens[donor, 0, 0]))[0]
    if len(live) == 0:
        return result
    h = hydrogens[donor[live]]                                               # (m, 2, 3)
    n, s = xyz[donor[live]][:, None, :], xyz[acceptor[live]][:, None, :]
    to_n, to_s = n - h, s - h
    h_s = np.linalg.norm(to_s, axis=-1)
    with np.errstate(invalid="ignore"):
        cosine = (to_n * to_s).sum(-1) / (np.linalg.norm(to_n, axis=-1) * h_s)
    ok = (h_s <= h_max) & (cosine <= np.cos(np.radians(protonation.SULFUR_ANGLE_MIN)))
    result[live] = ok.any(axis=1)
    return result


# ----------------------------------------------------------------------------- contacts
def residue_contacts(atoms: pd.DataFrame, *, cutoff: float = 4.0, min_distance: float = 1.9,
                     min_separation: int = 2, chemistry: bool = True, sulfur_acceptor: bool = False,
                     metals: Optional[pd.DataFrame] = None, neighbor_contacts: bool = False,
                     neighbor_chemistry: Sequence[str] = NEIGHBOR_CHEMISTRY) -> pd.DataFrame:
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

    ``sulfur_acceptor=True`` (opt-in; needs ``chemistry``) adds ``polar_s``, ``polar_s_bb`` and ``polar_s_geom``.
    ``polar_s``: a hydrogen-bond-capable donor (backbone N except Pro, or a side-chain N/OH listed in
    ``SULFUR_DONORS``) within ``SULFUR_ACCEPTORS`` of a Cys SG (3.7 Å) or Met SD (3.6 Å) in another residue,
    with side-chain donors capped at ``SULFUR_SIDE_CHAIN_MAX`` (3.6 Å). ``polar_s_bb`` keeps only backbone-N
    donors. ``polar_s_geom`` places the hydrogen from the heavy atoms (backbone amide, Trp NE1, Arg, Asn/Gln
    amide, and a His ring N as if protonated, either tautomer) and requires H...S <= 2.9 Å and N-H...S >= 120°;
    it does not test hydroxyl or Lys NZ donors, whose hydrogens rotate. A backbone amide needs the preceding
    residue in the table (seq_index - 1) to place its H. Evidence: the N-H...S survey in atomic-resolution
    structures (docs/scope_checks/03_sulfur_acceptor_survey.md). A disulfide SG is never an acceptor. A donor atom that a metal
    coordinates is deprotonated and is not a donor: ``metals`` (``structure_id, element, x, y, z``, see
    :func:`load_metals`) excludes every N/O within the ``bond`` tier of :mod:`modules.coordination_table`
    of a metal (a pair not in the table follows the global 2.8 Å); without ``metals`` nothing is
    excluded, so pass it for any structure with metal ions. Like ``polar``, a flag marks a
    capable pair, not a proven hydrogen bond, and ``polar``, ``hydrophobic`` and the other columns are
    unchanged. Sequence neighbors stay excluded by ``min_separation`` unless ``neighbor_contacts`` is set.

    ``neighbor_contacts=True`` (opt-in; needs ``chemistry`` and ``min_separation`` above 1) adds the contacts
    between residues one apart in sequence that the separation rule hides, judged per atom pair instead of
    per residue pair. A pair of atoms counts when it is at least ``NEIGHBOR_MIN_BONDS`` (4) bonds apart
    through the peptide bond (bonded 1-2, 1-3 and 1-4 atoms never count) and belongs to one of the classes in
    ``neighbor_chemistry`` (default ``polar`` and ``ionic``; ``NEIGHBOR_CHEMISTRY_ALL`` adds ``hydrophobic`` and
    ``aromatic``): ``polar`` is a hydrogen bond with at least one side-chain atom (a donor to an acceptor within
    3.5 Å, and where the donor's hydrogen can be placed it must point at the acceptor; a Pro N and an atom a metal
    coordinates never count), ``ionic`` a charged side-chain pair, ``hydrophobic`` a C/S pair with both atoms beyond
    C-beta, ``aromatic`` a ring atom against a side-chain atom beyond C-beta (see :func:`_neighbor_flags`; a ring beside the
    neighbor's backbone is not counted, it would mark every aromatic residue).
    A residue pair with such an atom pair gets a row with ``min_distance`` and ``closest`` of its closest qualifying
    atom pair, ``n_contacts`` the number of qualifying pairs, the flags of the selected classes, the others and the
    sulfur flags False, and the extra column ``neighbor`` True (False on every other row). Evidence:
    docs/scope_checks/05_sequence_neighbors_survey.md and 06_neighbor_chemistry_survey.md.
    """
    if sulfur_acceptor and not chemistry:
        raise ValueError("sulfur_acceptor needs chemistry=True")
    if neighbor_contacts and (not chemistry or min_separation <= 1):
        raise ValueError("neighbor_contacts needs chemistry=True and min_separation of at least 2")
    neighbor_chemistry = tuple(neighbor_chemistry)
    if neighbor_contacts and (not neighbor_chemistry or set(neighbor_chemistry) - set(NEIGHBOR_CHEMISTRY_ALL)):
        raise ValueError(f"neighbor_chemistry must be a non-empty subset of {NEIGHBOR_CHEMISTRY_ALL}, not {neighbor_chemistry!r}")
    columns = ["structure_id", "position_a", "position_b", "resname_a", "resname_b", "min_distance", "n_contacts"]
    if chemistry:
        columns += ["closest", *CONTACT_FLAGS, *(SULFUR_FLAGS if sulfur_acceptor else ()), *(("neighbor",) if neighbor_contacts else ())]
    frames: List[pd.DataFrame] = []
    for sid, frame in _per_structure(_heavy(atoms)):
        frame = frame.reset_index(drop=True)
        xyz = frame[["x", "y", "z"]].to_numpy(float)
        i, j, d = _close_pairs(xyz, xyz, cutoff)
        pos = frame["position"].to_numpy()
        seq = frame["seq_index"].to_numpy()
        name = frame["atom"].to_numpy()
        elem = frame["element"].str.upper().to_numpy()
        res = frame["resname"].str.upper().to_numpy()
        bb = np.array([n in BACKBONE for n in name])
        keep = (i < j) & (pos[i] != pos[j]) & (np.abs(seq[i] - seq[j]) >= min_separation)
        neighbor = np.zeros(int(keep.sum()), bool)
        if neighbor_contacts:
            coordinated = (_metal_coordinated(frame, xyz, metals[metals["structure_id"] == sid]) if metals is not None
                           else np.zeros(len(frame), bool))
            hydrogens = placed_hydrogens(frame, xyz)
            classes = _neighbor_flags(elem, res, name, bb, i, j, d, cutoff, xyz, hydrogens, coordinated)
            wanted = np.zeros(len(i), bool)
            for chosen in neighbor_chemistry:
                wanted |= classes[chosen]
            candidates = (i < j) & (pos[i] != pos[j]) & (np.abs(seq[i] - seq[j]) == 1) & wanted
            ni, nj = _neighbor_atom_pairs(xyz, pos, seq, name, bb, candidates, i, j)
            if len(ni):
                nd = np.linalg.norm(xyz[ni] - xyz[nj], axis=-1)
                i, j, d = np.concatenate([i[keep], ni]), np.concatenate([j[keep], nj]), np.concatenate([d[keep], nd])
                neighbor = np.concatenate([neighbor, np.ones(len(ni), bool)])
            else:
                i, j, d = i[keep], j[keep], d[keep]
        else:
            i, j, d = i[keep], j[keep], d[keep]
        if len(i) == 0:
            continue
        swap = pos[i] > pos[j]
        a, b = np.where(swap, j, i), np.where(swap, i, j)
        pair = pd.DataFrame({"position_a": pos[a], "position_b": pos[b], "d": d, "closest": np.where(bb[a], "bb", "sc").astype(object) + "_" + np.where(bb[b], "bb", "sc").astype(object)})
        if chemistry:
            polar = np.isin(elem[a], ["N", "O"]) & np.isin(elem[b], ["N", "O"]) & (d <= min(POLAR_DISTANCE, cutoff))
            ring = np.array([n in RING_ATOMS.get(r, ()) for n, r in zip(name, res)])
            posi = np.array([(r, n) in POSITIVE_ATOMS for n, r in zip(name, res)])
            nega = np.array([(r, n) in NEGATIVE_ATOMS for n, r in zip(name, res)])
            within = d <= min(CHEMISTRY_DISTANCE, cutoff)
            pair["polar"] = polar
            pair["polar_sc"] = polar & (~bb[a] | ~bb[b])
            pair["polar_bb_bb"] = polar & bb[a] & bb[b]
            pair["hydrophobic"] = np.isin(elem[a], ["C", "S"]) & np.isin(elem[b], ["C", "S"]) & ~bb[a] & ~bb[b] & within
            pair["aromatic"] = (ring[a] | ring[b]) & within
            pair["ionic"] = ((posi[a] & nega[b]) | (nega[a] & posi[b])) & within
            if sulfur_acceptor:
                acceptor_limit = np.array([SULFUR_ACCEPTORS.get((r, n), 0.0) for r, n in zip(res, name)])
                sg = np.nonzero((res == "CYS") & (name == "SG"))[0]
                disulfide = np.zeros(len(frame), bool)
                if len(sg) > 1:
                    si, sj, _ = _close_pairs(xyz[sg], xyz[sg], DISULFIDE_DISTANCE)
                    disulfide[sg[si[si != sj]]] = True
                acceptor_limit[disulfide] = 0.0
                backbone_donor = (name == "N") & (res != "PRO")
                side_donor = np.array([(r, n) in SULFUR_DONORS for r, n in zip(res, name)])
                coordinated = np.zeros(len(frame), bool)
                if metals is not None:
                    coordinated = _metal_coordinated(frame, xyz, metals[metals["structure_id"] == sid])
                    backbone_donor &= ~coordinated
                    side_donor &= ~coordinated
                side_limit = np.minimum(acceptor_limit, SULFUR_SIDE_CHAIN_MAX)
                pair["polar_s_bb"] = ((backbone_donor[a] & (d <= acceptor_limit[b])) | (backbone_donor[b] & (d <= acceptor_limit[a]))) & (d <= cutoff)
                pair["polar_s"] = pair["polar_s_bb"] | (side_donor[a] & (d <= side_limit[b])) | (side_donor[b] & (d <= side_limit[a]))
                hydrogens = placed_hydrogens(frame, xyz)
                hydrogens[coordinated] = np.nan
                pair["polar_s_geom"] = (_points_at(hydrogens, xyz, a, b, acceptor_limit > 0) | _points_at(hydrogens, xyz, b, a, acceptor_limit > 0)) & (d <= cutoff)
        if neighbor_contacts:
            pair["nbr"] = neighbor
            classes = _neighbor_flags(elem, res, name, bb, a, b, d, cutoff, xyz, hydrogens, coordinated)
            for flag in ("polar", "ionic", "hydrophobic", "aromatic"):
                own = classes[flag] if flag in neighbor_chemistry else np.zeros(len(d), bool)
                pair[flag] = np.where(neighbor, own, pair[flag])        # a neighbor row follows the neighbor rules and the chosen classes
            pair["polar_sc"] = np.where(neighbor, pair["polar"], pair["polar_sc"])
            pair["polar_bb_bb"] &= ~neighbor
            for flag in SULFUR_FLAGS if sulfur_acceptor else ():
                pair[flag] &= ~neighbor
        pair["ra"], pair["rb"] = res[a], res[b]
        keys = ["position_a", "position_b"]
        flags = (list(CONTACT_FLAGS) + (list(SULFUR_FLAGS) if sulfur_acceptor else [])) if chemistry else []
        closest = pair.sort_values("d", kind="stable").drop_duplicates(keys)[keys + ["d", "closest", "ra", "rb"]]
        extra = {"neighbor": ("nbr", "any")} if neighbor_contacts else {}
        agg = pair.groupby(keys, sort=False).agg(n_contacts=("d", "size"), **{f: (f, "any") for f in flags}, **extra).reset_index()
        out = closest.merge(agg, on=keys)
        out = out[out["d"] >= min_distance].rename(columns={"d": "min_distance", "ra": "resname_a", "rb": "resname_b"})
        out.insert(0, "structure_id", sid)
        frames.append(out)
    used = dict(cutoff=cutoff, min_distance=min_distance, min_separation=min_separation, sulfur_acceptor=sulfur_acceptor,
                neighbor_contacts=neighbor_contacts, neighbor_chemistry=neighbor_chemistry)
    if not frames:
        return _with_rules(pd.DataFrame(columns=columns), "residue_contacts", **used)
    return _with_rules(pd.concat(frames, ignore_index=True)[columns], "residue_contacts", **used)


# ----------------------------------------------------------------------------- cofactor features
def cofactor_contacts(atoms: pd.DataFrame, features: pd.DataFrame, *, coordination_distance: float = 2.8,
                      polar_distance: float = 3.5, ionic_distance: float = 4.0, contact_distance: float = 4.0,
                      record_distance: float = 8.0, coordination_model: Optional[str] = None) -> pd.DataFrame:
    """Contacts between residues and cofactor features, one row per structure x feature x position.

    Classes: ``coord`` (a ``metal`` feature atom to a residue N, O, or S within
    ``coordination_distance``), ``polar`` (an ``oxygen`` feature atom to a residue N or O
    within ``polar_distance``), ``ionic`` (an ``oxygen`` feature atom to Lys NZ or Arg
    NE/NH1/NH2 within ``ionic_distance``), ``contact`` (any heavy-atom pair within
    ``contact_distance``). A feature groups all its atoms; residues farther than
    ``record_distance`` from every atom of the feature are not recorded (absent).

    ``coordination_model="by_pair"`` (opt-in, experimental) adds ``coord_tier`` (``bond``, ``weak``,
    ``default`` or ``none``), ``coord_donor_class`` and ``coord_rule`` from
    :mod:`modules.coordination_table`: the cutoffs depend on the metal element and the donor
    atom, and a pair not in the table follows ``coordination_distance`` and is marked ``default``.
    The ``coord`` column is unchanged, so default output and the manuscript reproduction are identical.
    """
    if coordination_model not in (None, "by_pair"):
        raise ValueError(f"coordination_model must be None or 'by_pair', not {coordination_model!r}")
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
            metal_element = fgrp["element"].str.upper().to_numpy()
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
                if coordination_model:
                    donor = is_metal & np.isin(elem[ii], ["N", "O", "S", "SE"])
                    call = coordination_table.strongest(
                        (float(dist), coordination_table.coordination_tier(metal_element[b], res[a], name[a], float(dist),
                                                                           default_max=coordination_distance))
                        for a, b, dist in zip(ii[donor], jj[donor], dd[donor]))
                    row.update(zip(TIER_COLUMNS, call))
                rows.append(row)
    columns = ["structure_id", "feature", "position", "resname", "min_distance", *FEATURE_CLASSES]
    return _with_rules(pd.DataFrame(rows, columns=columns + (list(TIER_COLUMNS) if coordination_model else [])), "cofactor_contacts",
                       coordination_distance=coordination_distance, polar_distance=polar_distance, ionic_distance=ionic_distance,
                       contact_distance=contact_distance, record_distance=record_distance, coordination_model=coordination_model)


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
    return _with_rules(pd.DataFrame(rows, columns=["structure_id", "feature", "position", "bridged", "direct", "indirect"]), "water_bridges",
                       bridge=bridge, direct=direct)


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
    return _with_rules(pd.DataFrame(out, columns=["structure_id", "position", "steps"]), "relay_steps", link=link, max_steps=max_steps)


# ----------------------------------------------------------------------------- the vote
def sequence_votes(observations: pd.DataFrame, structures: pd.DataFrame, keys: Sequence[str], flag: str, *,
                   condition: Optional[str] = None, structure_fraction: float = 0.5) -> pd.DataFrame:
    """The vote of every sequence for every key: ``keys``, ``group``, ``sequence_id`` and a boolean ``vote``.

    The first half of :func:`sequence_vote`: within a (key, sequence) the vote is positive when at least
    ``structure_fraction`` of the sequence's structures (those present in ``observations`` for that key, and where
    ``condition`` is true when given) are positive for ``flag``. Exposed so that the stability of the second half (the support
    rule over sequences) can be studied without recomputing any geometry.
    """
    meta = structures[["structure_id", "group", "sequence_id"]].drop_duplicates("structure_id")
    data = observations.merge(meta, on="structure_id", how="inner")
    return _per_sequence_votes(data, list(keys), flag, condition, structure_fraction)


def _per_sequence_votes(data: pd.DataFrame, keys: Sequence[str], flag: str, condition: Optional[str], structure_fraction: float) -> pd.DataFrame:
    sub = data if condition is None else data[data[condition].astype(bool)]
    share = sub.assign(_flag=sub[flag].astype(float)).groupby(list(keys) + ["group", "sequence_id"], sort=False)["_flag"].mean()
    return (share >= structure_fraction).rename("vote").reset_index()


def vote_coverage(votes: pd.DataFrame, groups: Sequence[str], *, support: float = 0.90, min_sequences: int = 5) -> pd.DataFrame:
    """How deep the vote behind each key is, from the table of :func:`sequence_votes`.

    One row per key (the columns of ``votes`` other than ``group``, ``sequence_id`` and ``vote``) with, for each group ``g``:
    ``spare_negatives_<g>`` (how many further negative sequences the key could take and still pass: ``floor(n (1 - support))``
    minus its negatives, so negative when the group fails on support) and ``spare_voters_<g>`` (voting sequences above
    ``min_sequences``; negative when the group is too small to count). Also ``conserved_in`` (the groups, joined by ``+``,
    whose support rule the key passes), ``failing_groups`` (the others), and ``held_back_by``: when exactly one group fails and
    leaving out any single one of its negative sequences would make it pass, those sequences' ids joined by ``;``
    (a group that fails because it has too few voting sequences is never held back by a sequence). ``one_negative_away`` is the
    same without the minimum: exactly one group fails, and leaving out any single one of its negative sequences would bring its
    positive share up to ``support`` (the group may still be below the minimum number of voting sequences, as a group of exactly
    ``min_sequences`` always is once a sequence is left out). A key is conserved (:func:`sequence_vote`) exactly when ``failing_groups`` is empty.
    """
    groups = list(groups)
    keys = [c for c in votes.columns if c not in ("group", "sequence_id", "vote")]
    grouped = votes.groupby(keys + ["group"], sort=False)["vote"]
    n = grouped.count().unstack("group").reindex(columns=groups).fillna(0).astype(int)
    pos = grouped.sum().unstack("group").reindex(columns=groups).fillna(0).astype(int)
    out = pd.DataFrame(index=n.index)
    passing, rescued, by_support_only = {}, {}, {}
    for g in groups:
        negatives = n[g] - pos[g]
        out[f"spare_negatives_{g}"] = np.floor(n[g] * (1 - support) + 1e-9).astype(int) - negatives
        out[f"spare_voters_{g}"] = n[g] - min_sequences
        passing[g] = (n[g] >= min_sequences) & (pos[g] / n[g].where(n[g] > 0) >= support)
        by_support = (negatives >= 1) & (pos[g] / (n[g] - 1).where(n[g] > 1) >= support)
        rescued[g] = by_support & (n[g] - 1 >= min_sequences)
        by_support_only[g] = by_support
    passes = pd.DataFrame(passing)
    out["conserved_in"] = ["+".join(g for g in groups if row[g]) for row in passes.to_dict("records")]
    out["failing_groups"] = ["+".join(g for g in groups if not row[g]) for row in passes.to_dict("records")]
    out["held_back_by"] = ""
    out["one_negative_away"] = ""
    negative_ids = votes[~votes["vote"].astype(bool)].groupby(keys + ["group"], sort=False)["sequence_id"].apply(lambda s: ";".join(sorted(s)))
    single = (~passes).sum(axis=1) == 1
    for g in groups:
        for column, rescue in (("held_back_by", rescued[g]), ("one_negative_away", by_support_only[g])):
            for key in out.index[single & ~passes[g] & rescue]:
                out.loc[key, column] = negative_ids.get((*key, g) if isinstance(key, tuple) else (key, g), "")
    return out.reset_index()


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
        per_seq = _per_sequence_votes(data, keys, flag, condition, structure_fraction)
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


def contact_observations(atoms: pd.DataFrame, *, cutoff: float = 4.0, min_separation: int = 2, sulfur_acceptor: bool = False,
                         metals: Optional[pd.DataFrame] = None, neighbor_contacts: bool = False,
                         neighbor_chemistry: Sequence[str] = NEIGHBOR_CHEMISTRY) -> pd.DataFrame:
    """One row per (structure, candidate pair): ``contact`` plus the chemistry of the pair where it is in contact.

    The candidates are the position pairs in contact anywhere in the set, evaluated in every structure that maps both
    positions (a structure with both positions present and no contact is a negative). This is the table that
    :func:`conserved_contacts` votes on.
    """
    contacts = residue_contacts(atoms, cutoff=cutoff, min_separation=min_separation, sulfur_acceptor=sulfur_acceptor, metals=metals,
                                neighbor_contacts=neighbor_contacts, neighbor_chemistry=neighbor_chemistry)
    present = _heavy(atoms).groupby("structure_id")["position"].apply(set)
    candidates = list(contacts[["position_a", "position_b"]].drop_duplicates().itertuples(index=False, name=None))
    in_contact = set(zip(contacts["structure_id"], contacts["position_a"], contacts["position_b"]))
    rows = [(sid, a, b, (sid, a, b) in in_contact) for sid, have in present.items() for a, b in candidates if a in have and b in have]
    obs = pd.DataFrame(rows, columns=["structure_id", "position_a", "position_b", "contact"])
    obs = obs.merge(contacts, on=["structure_id", "position_a", "position_b"], how="left")
    obs["contact"] = obs["contact"].astype(bool)
    for c in CONTACT_CLASSES:
        obs[c] = obs["closest"] == c
    return obs


def conserved_contacts(atoms: pd.DataFrame, structures: pd.DataFrame, *, groups: Optional[Sequence[str]] = None,
                       cutoff: float = 4.0, min_separation: int = 2, support: float = 0.90, min_sequences: int = 5,
                       structure_fraction: float = 0.5, sulfur_acceptor: bool = False,
                       metals: Optional[pd.DataFrame] = None, neighbor_contacts: bool = False,
                       neighbor_chemistry: Sequence[str] = NEIGHBOR_CHEMISTRY, coverage: bool = False) -> pd.DataFrame:
    """Contacts passing the sequence vote, with the conserved chemistry of each.

    Every pair of positions mapped in a structure is a candidate; the contact vote counts
    a structure as negative when both positions are present but not in contact. Chemistry
    flags are voted over contact-positive structures only (``closest`` classes and the
    flags in :data:`CONTACT_FLAGS`). ``sulfur_acceptor`` and ``metals`` are passed to
    :func:`residue_contacts` (opt-in), which adds the ``polar_s*`` flags to the voted chemistry;
    ``neighbor_contacts`` adds the polar and ionic contacts between sequence neighbors, voted like any other
    pair (a structure is negative when both positions are present and the pair has no such atom pair);
    ``neighbor_chemistry`` selects its classes. ``coverage=True`` (opt-in) adds the depth of the vote behind every contact from
    :func:`vote_coverage`: ``conserved_in``, ``failing_groups``, ``held_back_by``, ``one_negative_away`` and the spare negatives and voters per
    group; the default columns and the conserved set do not change.
    """
    obs = contact_observations(atoms, cutoff=cutoff, min_separation=min_separation, sulfur_acceptor=sulfur_acceptor, metals=metals,
                               neighbor_contacts=neighbor_contacts, neighbor_chemistry=neighbor_chemistry)
    keys = ["position_a", "position_b"]
    base = sequence_vote(obs, structures, keys, ["contact"], groups=groups, support=support,
                         min_sequences=min_sequences, structure_fraction=structure_fraction)
    chem = sequence_vote(obs, structures, keys, list(CONTACT_CLASSES) + list(CONTACT_FLAGS) + (list(SULFUR_FLAGS) if sulfur_acceptor else []),
                         condition="contact", groups=groups,
                         support=support, min_sequences=min_sequences, structure_fraction=structure_fraction)
    result = base.merge(chem, on=keys, how="left")
    if coverage:
        used_groups = list(groups) if groups is not None else sorted(structures["group"].unique())
        votes = sequence_votes(obs, structures, keys, "contact", structure_fraction=structure_fraction)
        result = result.merge(vote_coverage(votes, used_groups, support=support, min_sequences=min_sequences), on=keys, how="left")
    return _with_rules(result, "conserved_contacts", coverage=coverage, cutoff=cutoff, min_separation=min_separation,
                       sulfur_acceptor=sulfur_acceptor, neighbor_contacts=neighbor_contacts, neighbor_chemistry=neighbor_chemistry, support=support, min_sequences=min_sequences, structure_fraction=structure_fraction)


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


def backbone_only_audit(atoms: pd.DataFrame, *, max_heavy_atoms: int = 5, min_fraction: float = 0.5) -> pd.DataFrame:
    """Flag structures that are backbone-only models (side chains absent beyond CB).

    A structure is flagged when more than ``min_fraction`` of its mapped residues other than
    Gly and Ala have at most ``max_heavy_atoms`` distinct heavy atoms. Such a model cannot show a
    side-chain contact, so every side-chain contact reads as absent for its sequence. The
    fraction is 1.0 when a structure maps no residue other than Gly or Ala. Returns
    ``structure_id, n_residues, fraction_truncated, backbone_only``.
    """
    heavy = _heavy(atoms)
    heavy = heavy[~heavy["resname"].str.upper().isin(["GLY", "ALA"])]
    counts = heavy.groupby(["structure_id", "position"], sort=False)["atom"].nunique()
    truncated = (counts <= max_heavy_atoms).groupby(level="structure_id", sort=False)
    table = pd.DataFrame({"n_residues": truncated.size(), "fraction_truncated": truncated.mean()})
    missing = pd.Index(atoms["structure_id"].unique()).difference(table.index)
    if len(missing):
        table = pd.concat([table, pd.DataFrame({"n_residues": 0, "fraction_truncated": 1.0}, index=missing)])
    table["backbone_only"] = table["fraction_truncated"] > min_fraction
    table.index.name = "structure_id"
    return table.reset_index().sort_values("structure_id", ignore_index=True)


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


# ----------------------------------------------------------------------------- provenance
RULE_TABLES = ("residue_contacts", "cofactor_contacts", "water_bridges", "relay_steps", "conserved_contacts")


def class_rules(table: Optional[str] = None, **parameters) -> pd.DataFrame:
    """The rule behind every reported column, with the cutoffs actually used.

    One row per ``(table, name)`` with ``rule`` (a sentence), ``cutoff_A`` (the distance that decides it, when
    there is one) and ``status`` (``default`` or ``opt-in``). ``table`` is one of :data:`RULE_TABLES`, or None for
    all. ``parameters`` are the arguments of the function that produced the table (an omitted one is that
    function's default, read from its signature), so the text and the numbers cannot drift from the code. Every
    result of this module carries its rules in ``frame.attrs["rules"]`` (see :func:`rules_of`).
    """
    if table is not None and table not in RULE_TABLES:
        raise ValueError(f"table must be one of {RULE_TABLES} or None, not {table!r}")
    rows: List[dict] = []

    def add(owner, name, rule, cutoff=np.nan, status="default"):
        rows.append({"table": owner, "name": name, "rule": rule, "cutoff_A": cutoff, "status": status})

    def given(function, name):
        value = parameters.get(name)
        return inspect.signature(function).parameters[name].default if value is None else value

    if table in (None, "residue_contacts", "conserved_contacts"):
        cutoff, min_distance = given(residue_contacts, "cutoff"), given(residue_contacts, "min_distance")
        min_separation = given(residue_contacts, "min_separation")
        polar, chem = min(POLAR_DISTANCE, cutoff), min(CHEMISTRY_DISTANCE, cutoff)
        owner = "residue_contacts"
        add(owner, "contact", f"closest heavy-atom pair of two residues between {min_distance:g} and {cutoff:g} A, residues at least "
            f"{min_separation} apart in sequence; one row per residue pair", cutoff)
        add(owner, "closest", f"backbone/side-chain class (bb_bb, bb_sc, sc_bb, sc_sc) of the single closest atom pair; backbone atoms are {', '.join(sorted(BACKBONE))}")
        add(owner, "polar", f"any N/O atom pair within {polar:g} A (a capable pair judged from heavy atoms alone, and sulfur is never counted)", polar)
        add(owner, "polar_sc", f"a polar pair with at least one side-chain atom, within {polar:g} A", polar)
        add(owner, "polar_bb_bb", f"a polar pair with both atoms in the backbone, within {polar:g} A", polar)
        add(owner, "hydrophobic", f"a side-chain C/S pair within {chem:g} A", chem)
        add(owner, "aromatic", f"a ring atom of Phe, Tyr, Trp or His of either residue within {chem:g} A of any heavy atom of the other", chem)
        add(owner, "ionic", f"Lys NZ or Arg NE/NH1/NH2 to Asp OD1/OD2 or Glu OE1/OE2 within {chem:g} A", chem)
        windows = ", ".join(f"{r.title()} {a} {limit:g} A" for (r, a), limit in SULFUR_ACCEPTORS.items())
        opt = "opt-in"
        add(owner, "polar_s", f"a backbone N (not Pro) within {windows} of a sulfur in another residue, or a listed side-chain N/OH "
            f"(Ser, Thr, Tyr, Asn, Gln, Arg, Lys, Trp, His) within {SULFUR_SIDE_CHAIN_MAX:g} A; a disulfide SG (another SG within "
            f"{DISULFIDE_DISTANCE:g} A) never accepts; a donor a metal coordinates is excluded when `metals` is given; a capable pair, not a proven hydrogen bond",
            SULFUR_SIDE_CHAIN_MAX, opt)
        add(owner, "polar_s_bb", f"polar_s restricted to backbone N donors, within {windows}", max(SULFUR_ACCEPTORS.values()), opt)
        add(owner, "polar_s_geom", f"the hydrogen placed from the heavy atoms (backbone amide, Trp NE1, Arg, Asn/Gln amide, a His ring N as if protonated) "
            f"points at the sulfur: H...S <= {protonation.SULFUR_H_MAX:g} A and N-H...S >= {protonation.SULFUR_ANGLE_MIN:g} degrees; "
            f"hydroxyl and Lys NZ are not tested", protonation.SULFUR_H_MAX, opt)
        selected = tuple(parameters.get("neighbor_chemistry") or NEIGHBOR_CHEMISTRY)
        add(owner, "neighbor", f"opt-in (neighbor_contacts): residues one apart in sequence, judged per atom pair; an atom pair counts when it is at least "
            f"{NEIGHBOR_MIN_BONDS} bonds apart through the peptide bond and belongs to a selected class (default {', '.join(NEIGHBOR_CHEMISTRY)}; "
            f"all: {', '.join(NEIGHBOR_CHEMISTRY_ALL)}); a row carries only the selected classes, the sulfur flags are False, and min_distance and closest "
            f"are those of the closest qualifying pair", np.nan, opt)
        donors = "backbone N except Pro, or Ser OG, Thr OG1, Tyr OH, Asn ND2, Gln NE2, Arg N, Lys NZ, Trp NE1, a His ring N"
        add(owner, "neighbor_polar", f"a hydrogen-bond candidate with at least one side-chain atom: a donor ({donors}) to an acceptor (any O or a His ring N) within "
            f"{polar:g} A; when the donor's hydrogen can be placed it must point at the acceptor (H...A <= {NEIGHBOR_H_MAX:g} A, N-H...A >= "
            f"{protonation.SULFUR_ANGLE_MIN:g} degrees), otherwise (hydroxyl, Lys NZ) the distance decides; a metal-coordinated atom does not count", polar, opt)
        add(owner, "neighbor_ionic", f"Lys NZ or Arg NE/NH1/NH2 to Asp OD1/OD2 or Glu OE1/OE2 within {chem:g} A", chem, opt)
        add(owner, "neighbor_hydrophobic", f"a C/S pair, both atoms side-chain atoms beyond C-beta, within {chem:g} A (backbone and C-beta proximity is conformation)", chem, opt)
        add(owner, "neighbor_aromatic", f"a ring atom of Phe, Tyr, Trp or His against a side-chain atom beyond C-beta within {chem:g} A (a ring beside the backbone is not counted)", chem, opt)
        for row in rows:
            if row["name"].startswith("neighbor_") and parameters.get("neighbor_contacts", False) and row["name"][len("neighbor_"):] not in selected:
                row["rule"] += " [not selected: neighbor_chemistry excludes it]"
        if not parameters.get("neighbor_contacts", False) and table in ("residue_contacts", "conserved_contacts"):
            for row in rows:
                if row["name"].startswith("neighbor"):
                    row["rule"] += " [not computed in this table: neighbor_contacts was off]"
        if not parameters.get("sulfur_acceptor", False) and table in ("residue_contacts", "conserved_contacts"):
            for row in rows:
                if row["status"] == opt and row["name"].startswith("polar_s"):
                    row["rule"] += " [not computed in this table: sulfur_acceptor was off]"
    if table in (None, "cofactor_contacts"):
        coord = given(cofactor_contacts, "coordination_distance")
        polar, ionic = given(cofactor_contacts, "polar_distance"), given(cofactor_contacts, "ionic_distance")
        contact, record = given(cofactor_contacts, "contact_distance"), given(cofactor_contacts, "record_distance")
        owner = "cofactor_contacts"
        add(owner, "coord", f"a metal feature atom to a residue N, O or S within {coord:g} A (one global cutoff)", coord)
        add(owner, "polar", f"an oxygen feature atom to a residue N or O within {polar:g} A", polar)
        add(owner, "ionic", f"an oxygen feature atom to Lys NZ or Arg NE/NH1/NH2 within {ionic:g} A", ionic)
        add(owner, "contact", f"any feature atom to any residue heavy atom within {contact:g} A", contact)
        add(owner, "(rows)", f"residues farther than {record:g} A from every atom of the feature are not recorded", record)
        if table is None or parameters.get("coordination_model"):
            opt = "opt-in"
            add(owner, "coord_tier", f"by metal and donor class (coordination_model='by_pair'): bond, weak, default (pair not in the table: "
                f"follows {coord:g} A) or none; the row used is in coord_rule", status=opt)
            for entry in coordination_table.RULES:
                weak = f", weak <= {entry.weak_max:g}" if entry.weak_max is not None else ""
                add(owner, f"coord_tier: {entry.metal.title()} {entry.donor_class}", f"bond <= {entry.bond_max:g}{weak} A; {entry.source}", entry.bond_max, opt)
    if table in (None, "water_bridges"):
        bridge, direct = given(water_bridges, "bridge"), given(water_bridges, "direct")
        add("water_bridges", "bridged", f"a modeled water oxygen within {bridge:g} A of a feature oxygen and within {bridge:g} A of a residue N/O", bridge)
        add("water_bridges", "direct", f"a residue N/O within {direct:g} A of the feature oxygen", direct)
        add("water_bridges", "indirect", "bridged and not direct")
    if table in (None, "relay_steps"):
        link, steps = given(relay_steps, "link"), given(relay_steps, "max_steps")
        add("relay_steps", "steps", f"fewest links from the source atoms to a residue through residue N/O atoms and waters; two nodes link within {link:g} A; "
            f"a residue reached only by its backbone N or O is a terminal; at most {steps} links", link)
    if table in (None, "conserved_contacts"):
        fraction, support = given(sequence_vote, "structure_fraction"), given(sequence_vote, "support")
        minimum = given(sequence_vote, "min_sequences")
        coverage_notes = [
            ("conserved_in", f"the groups whose support rule the contact passes (at least {support:.0%} of at least {minimum} voting sequences positive), joined by '+'; "
                             "the contact is conserved when all groups pass"),
            ("failing_groups", "the groups it does not pass"),
            ("spare_negatives_<group>", f"further negative sequences the contact could take in that group and still pass: floor(n x {1 - support:.0%}) minus its negatives; "
                                        "negative when the group fails on support"),
            ("spare_voters_<group>", f"voting sequences in that group above the minimum of {minimum}; negative when the group is too small to count"),
            ("held_back_by", "when exactly one group fails and leaving out any single one of its negative sequences would make it pass, those sequence ids; "
                             "a group that is too small is never held back by a sequence"),
            ("one_negative_away", "as held_back_by but ignoring the minimum number of voting sequences: leaving out one negative sequence would bring the "
                                  "failing group's positive share up to the support (a group at the minimum is always below it afterwards)"),
        ]
        for name, text in coverage_notes:
            note = "" if parameters.get("coverage", False) else " [not computed in this table: coverage was off]"
            add("conserved_contacts", name, text + note, np.nan, "opt-in")
        minimum = given(sequence_vote, "min_sequences")
        add("conserved_contacts", "conserved_<flag>", f"one vote per sequence: positive when at least {fraction:.0%} of the sequence's structures show it; "
            f"conserved when every group has at least {minimum} voting sequences and at least {support:.0%} of them positive; chemistry flags are voted "
            f"over contact-positive structures only", support)
    return pd.DataFrame(rows, columns=["table", "name", "rule", "cutoff_A", "status"])


def _with_rules(frame: pd.DataFrame, table: str, **parameters) -> pd.DataFrame:
    frame.attrs["rules"] = class_rules(table, **parameters).to_dict("records")
    return frame


def rules_of(frame: pd.DataFrame) -> pd.DataFrame:
    """The rules a result of this module was produced with (empty if the frame did not come from here)."""
    return pd.DataFrame(frame.attrs.get("rules", []), columns=["table", "name", "rule", "cutoff_A", "status"])


def format_rules(rules: pd.DataFrame) -> str:
    """The rules as a Markdown table, for a log or a report."""
    lines = ["| table | name | rule | cutoff (A) | status |", "|---|---|---|---|---|"]
    for r in rules.itertuples():
        cutoff = "" if pd.isna(r.cutoff_A) else f"{r.cutoff_A:g}"
        lines.append(f"| {r.table} | {r.name} | {r.rule} | {cutoff} | {r.status} |")
    return "\n".join(lines) + "\n"
