"""Opt-in, chemically typed atom-level interaction evidence.

This layer wraps PDBe Arpeggio (CREDO interaction rules) and deliberately
does not alter CNI's distance-defined shell contacts.  It writes independent
evidence and residue-pair tables so that the two analyses can be compared.
"""
from __future__ import annotations

import csv
import json
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

import numpy as np

ARPEGGIO_VERSION = "1.4.4"
ARPEGGIO_COMMIT = "258855b8ba13447f2776b232ca32884d637c6a9c"
TYPE_MAP = {
    "hbond": "hydrogen_bond",
    "weak_hbond": "weak_hydrogen_bond",
    "ionic": "ionic",
    "hydrophobic": "hydrophobic",
    "aromatic": "aromatic_atom_contact",
    "metal_complex": "metal_coordination",
    "covalent": "covalent",
    "CATIONPI": "cation_pi",
    "DONORPI": "donor_pi",
    "METSULPHURPI": "sulfur_pi",
}
ATOM_TYPES = frozenset(TYPE_MAP.values()) | frozenset({"disulfide", "aromatic_stacking"})
WATER_NAMES = frozenset({"HOH", "WAT", "DOD", "H2O"})

EVIDENCE_COLUMNS = [
    "structure_id", "analysis_layer", "interaction_type", "criterion",
    "src_resname", "src_chain", "src_resnum", "src_insertion_code",
    "src_atom", "src_element", "src_altloc", "dst_resname", "dst_chain",
    "dst_resnum", "dst_insertion_code", "dst_atom", "dst_element",
    "dst_altloc", "distance_A", "donor_h_acceptor_angle_deg",
    "plane_angle_deg", "centroid_vector_angle_deg", "angle_note", "criterion_source",
]
PAIR_COLUMNS = [
    "structure_id", "analysis_layer", "src_resname", "src_chain", "src_resnum",
    "src_insertion_code", "dst_resname", "dst_chain", "dst_resnum",
    "dst_insertion_code", "baseline_contact", "interaction_types", "outcome",
]


def angle_degrees(a: Sequence[float], vertex: Sequence[float], c: Sequence[float]) -> float:
    """Return the A-vertex-C angle in degrees, with a stable zero check."""
    va = np.asarray(a, dtype=float) - np.asarray(vertex, dtype=float)
    vc = np.asarray(c, dtype=float) - np.asarray(vertex, dtype=float)
    denom = float(np.linalg.norm(va) * np.linalg.norm(vc))
    if denom == 0:
        raise ValueError("angle points must be distinct")
    return math.degrees(math.acos(float(np.clip(np.dot(va, vc) / denom, -1.0, 1.0))))


def geometry_rule_accepts(
    interaction_type: str, *, distance_A: float, chemistry_supported: bool,
    donor_h_acceptor_angle_deg: Optional[float] = None,
    substituent_angle_deg: Optional[float] = None,
    plane_angle_deg: Optional[float] = None,
    centroid_vector_angle_deg: Optional[float] = None,
    covalent_bond: bool = False,
    cysteine_sg_pair: bool = False,
) -> bool:
    """Small QA oracle for the pinned geometric limits (not a replacement classifier)."""
    kind = interaction_type
    if not chemistry_supported:
        return False
    if kind == "hydrogen_bond":
        return distance_A <= 3.9 and donor_h_acceptor_angle_deg is not None and donor_h_acceptor_angle_deg >= 90
    if kind == "weak_hydrogen_bond":
        return (distance_A <= 3.6 and donor_h_acceptor_angle_deg is not None and
                donor_h_acceptor_angle_deg >= 130 and substituent_angle_deg is not None and
                30 <= substituent_angle_deg <= 150)
    if kind == "ionic": return distance_A <= 4.0
    if kind == "hydrophobic": return distance_A <= 4.5
    if kind == "aromatic_atom_contact": return distance_A <= 4.0
    if kind == "aromatic_stacking": return distance_A <= 6.0 and plane_angle_deg is not None and centroid_vector_angle_deg is not None and plane_angle_deg <= 90 and centroid_vector_angle_deg <= 90
    if kind in {"cation_pi", "donor_pi"}:
        return distance_A <= 4.5 and centroid_vector_angle_deg is not None and centroid_vector_angle_deg <= 30
    if kind == "sulfur_pi": return distance_A <= 6.0
    if kind == "metal_coordination": return distance_A <= 2.8
    if kind == "disulfide": return covalent_bond and cysteine_sg_pair
    if kind == "covalent": return covalent_bond
    return False


def interaction_type_eligibility(residue_a: Mapping[str, Any], residue_b: Mapping[str, Any]) -> dict[str, bool]:
    """Estimate which type rules are assessable from observed endpoint atoms.

    Inputs need `resname` and either `atom_names` or `atoms`. This is a
    denominator helper only; Arpeggio remains responsible for positive labels.
    It does not count a missing atom as a negative interaction observation.
    """
    def inventory(residue: Mapping[str, Any]) -> tuple[str, set[str]]:
        name = _clean(residue.get("resname")).upper()
        atoms = residue.get("atom_names", residue.get("atoms", []))
        if isinstance(atoms, Mapping):
            atoms = atoms.keys()
        names = {str(atom.get("name", "") if isinstance(atom, Mapping) else atom).strip().upper() for atom in atoms}
        names.discard("")
        return name, names

    name_a, atoms_a = inventory(residue_a)
    name_b, atoms_b = inventory(residue_b)
    donors = {"N", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "OG", "OG1", "OH", "SG"}
    acceptors = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2", "SD"}
    positive_atoms = {"NZ", "NE", "NH1", "NH2", "ND1", "NE2"}
    negative_atoms = {"OD1", "OD2", "OE1", "OE2"}
    aromatic_atoms = {
        "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
        "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
        "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
        "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    }
    def donor_atoms(name: str, atoms: set[str]) -> set[str]:
        valid = atoms & donors
        if name == "PRO":
            valid.discard("N")
        return valid
    def acceptor_atoms(atoms: set[str]) -> set[str]:
        return atoms & acceptors
    donor_a, donor_b = donor_atoms(name_a, atoms_a), donor_atoms(name_b, atoms_b)
    acceptor_a, acceptor_b = acceptor_atoms(atoms_a), acceptor_atoms(atoms_b)
    ring_a = aromatic_atoms.get(name_a, set()) <= atoms_a and name_a in aromatic_atoms
    ring_b = aromatic_atoms.get(name_b, set()) <= atoms_b and name_b in aromatic_atoms
    positive_a = bool(atoms_a & positive_atoms) and name_a in {"LYS", "ARG", "HIS"}
    positive_b = bool(atoms_b & positive_atoms) and name_b in {"LYS", "ARG", "HIS"}
    negative_a = bool(atoms_a & negative_atoms) and name_a in {"ASP", "GLU"}
    negative_b = bool(atoms_b & negative_atoms) and name_b in {"ASP", "GLU"}
    carbon_a = atoms_a - {"C", "CA", "O", "N", "OXT"}
    carbon_b = atoms_b - {"C", "CA", "O", "N", "OXT"}
    carbonyl_sidechain = {"ASP": {"CG"}, "ASN": {"CG"}, "GLU": {"CD"}, "GLN": {"CD"}}
    hydrophobic_a = {name for name in carbon_a if name.startswith("C") and name not in carbonyl_sidechain.get(name_a, set())}
    hydrophobic_b = {name for name in carbon_b if name.startswith("C") and name not in carbonyl_sidechain.get(name_b, set())}
    return {
        "hydrogen_bond": bool((donor_a and acceptor_b) or (donor_b and acceptor_a)),
        "weak_hydrogen_bond": bool((carbon_a and acceptor_b) or (carbon_b and acceptor_a)),
        "ionic": bool((positive_a and negative_b) or (positive_b and negative_a)),
        "hydrophobic": bool(hydrophobic_a and hydrophobic_b),
        "aromatic_atom_contact": ring_a and ring_b,
        "aromatic_stacking": ring_a and ring_b,
        "cation_pi": bool((positive_a and ring_b) or (positive_b and ring_a)),
        "donor_pi": bool((donor_a and ring_b) or (donor_b and ring_a)),
        "sulfur_pi": bool((name_a == "MET" and "SD" in atoms_a and ring_b) or (name_b == "MET" and "SD" in atoms_b and ring_a)),
        "disulfide": name_a == "CYS" and name_b == "CYS" and "SG" in atoms_a and "SG" in atoms_b,
        "covalent": name_a == "CYS" and name_b == "CYS" and "SG" in atoms_a and "SG" in atoms_b,
        "metal_coordination": False,
    }


def _residue_inventory(path: Path) -> dict[tuple[str, str, str, str], dict[str, Any]]:
    from Bio.PDB import MMCIFParser, PDBParser
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    structure = parser.get_structure("typed_inventory", str(path))
    result = {}
    for chain in next(structure.get_models()):
        for residue in chain:
            heavy = [atom for atom in residue.get_unpacked_list() if str(atom.element or "").upper() not in {"H", "D"}]
            key = (str(chain.id), str(residue.id[1]), str(residue.id[2]).strip(), residue.resname.upper())
            if heavy:
                result[key] = {"resname": residue.resname.upper(), "atom_names": {atom.name.strip().upper() for atom in heavy}}
    return result


def _complete_eligibility(edge: dict[str, Any], inventory: Mapping[tuple[str, str, str, str], Mapping[str, Any]]) -> None:
    def get_endpoint(prefix: str):
        key = tuple(_clean(edge.get(f"{prefix}_{field}")) for field in ("chain", "resnum", "insertion_code", "resname"))
        return inventory.get(key)
    a, b = get_endpoint("src"), get_endpoint("dst")
    try:
        supplied = json.loads(edge.get("eligible_by_type") or "{}") if isinstance(edge.get("eligible_by_type"), str) else edge.get("eligible_by_type") or {}
    except json.JSONDecodeError:
        supplied = {}
    computed = interaction_type_eligibility(a, b) if a and b else {kind: False for kind in ATOM_TYPES}
    computed.update({str(key): bool(value) for key, value in supplied.items()})
    edge["eligible_by_type"] = computed
    if edge.get("eligible") in (None, ""):
        edge["eligible"] = bool(a and b)


def _clean(value: Any) -> str:
    if value is None or str(value) in {".", "?", "None", "nan"}:
        return ""
    return str(value).strip()


def _endpoint(atom: Mapping[str, Any], element_lookup: Mapping[tuple, Any]) -> dict[str, str]:
    comp = _clean(atom.get("label_comp_id")).upper()
    chain = _clean(atom.get("auth_asym_id"))
    number = _clean(atom.get("auth_seq_id"))
    insertion = _clean(atom.get("pdbx_PDB_ins_code"))
    name = _clean(atom.get("auth_atom_id"))
    atom_names = [part.strip() for part in name.split(",") if part.strip()]
    atom_info = [element_lookup.get((chain, number, insertion, comp, atom_name), {}) for atom_name in atom_names]
    normalized = []
    for value in atom_info:
        if isinstance(value, str):
            value = {"element": value, "altloc": ""}
        normalized.append(value)
    elements = [value.get("element", "") for value in normalized if value.get("element")]
    altlocs = sorted({value.get("altloc", "") for value in normalized if value.get("altloc")})
    return {
        "resname": comp, "chain": chain, "resnum": number,
        "insertion_code": insertion, "atom": name,
        "element": ",".join(elements),
        "altloc": "" if not altlocs else altlocs[0] if len(altlocs) == 1 else "ambiguous:" + "|".join(altlocs),
    }


def _layer(a: Mapping[str, Any], b: Mapping[str, Any]) -> str:
    names = {str(a.get("label_comp_type", "")).lower(), str(b.get("label_comp_type", "")).lower()}
    residues = {str(a.get("label_comp_id", "")).upper(), str(b.get("label_comp_id", "")).upper()}
    if names & {"water", "solvent", "w"} or residues & WATER_NAMES:
        return "water_mediated_or_water_endpoint"
    if names <= {"p", "polypeptide(l)", "polypeptide(d)"}:
        return "direct_protein_protein"
    if not names or "non-polymer" in names or "branched" in names or "l" in names:
        return "protein_cofactor_or_ligand"
    return "protein_cofactor_or_ligand"


def _criteria(contact: str, a: Mapping[str, Any], b: Mapping[str, Any]) -> str:
    if contact == "hydrogen_bond":
        return "Arpeggio 1.4.4 CREDO hbond typing; donor-acceptor <=3.9 A and donor-H-acceptor >=90 degrees"
    if contact == "weak_hydrogen_bond":
        return "Arpeggio 1.4.4 CREDO weak-hbond geometry (distance and angular criteria in pinned config)"
    if contact == "ionic":
        return "Arpeggio 1.4.4 oppositely ionisable atom types, distance <=4.0 A"
    if contact == "hydrophobic":
        return "Arpeggio 1.4.4 hydrophobe atom types, distance <=4.5 A"
    if contact == "aromatic_atom_contact":
        return "Arpeggio 1.4.4 aromatic atom types and atom-atom distance <=4.0 A"
    if contact == "aromatic_stacking":
        return "Arpeggio 1.4.4 aromatic ring centroids <=6.0 A, with plane geometry classified by dihedral and centroid-vector angles"
    if contact == "metal_coordination":
        return "Arpeggio 1.4.4 metal and acceptor atom types, distance <=2.8 A"
    if contact == "covalent":
        if _clean(a.get("label_comp_id")).upper() == _clean(b.get("label_comp_id")).upper() == "CYS":
            return "Arpeggio covalent bond; cysteine SG-SG evidence reported as disulfide candidate"
        return "Arpeggio 1.4.4 Open Babel perceived covalent bond"
    if contact == "cation_pi":
        return "Arpeggio CREDO cation-ring centroid distance <=4.5 A and ring-normal theta <=30 degrees"
    if contact == "donor_pi":
        return "Arpeggio CREDO donor-ring geometry; centroid distance <=4.5 A and ring-normal theta <=30 degrees"
    if contact == "sulfur_pi":
        return "Arpeggio CREDO methionine sulfur-aromatic centroid distance <=6.0 A"
    return "Arpeggio 1.4.4 typed geometry"


def normalize_arpeggio_contacts(
    contacts: Iterable[Mapping[str, Any]], *, structure_id: str = "",
    element_lookup: Optional[Mapping[tuple, Any]] = None,
) -> list[dict[str, Any]]:
    """Normalize Arpeggio JSON to one row per atom-pair type, preserving layers."""
    elements = element_lookup or {}
    rows = []
    for item in contacts:
        a, b = item.get("bgn", {}), item.get("end", {})
        if not isinstance(a, Mapping) or not isinstance(b, Mapping):
            continue
        layer = _layer(a, b)
        if item.get("type") == "atom-atom":
            raw_types = item.get("contact", [])
            for raw in raw_types:
                interaction = TYPE_MAP.get(str(raw).lower())
                if not interaction:
                    continue
                # A protein-protein disulfide is a covalent subtype, never a generic proximity label.
                if interaction == "covalent" and str(a.get("auth_atom_id", "")).upper() == str(b.get("auth_atom_id", "")).upper() == "SG" and str(a.get("label_comp_id", "")).upper() == str(b.get("label_comp_id", "")).upper() == "CYS":
                    interaction = "disulfide"
                ea, eb = _endpoint(a, elements), _endpoint(b, elements)
                rows.append({
                    "structure_id": structure_id, "analysis_layer": layer,
                    "interaction_type": interaction,
                    "criterion": _criteria("covalent" if interaction == "disulfide" else interaction, a, b),
                    **{f"src_{k}": v for k, v in ea.items()},
                    **{f"dst_{k}": v for k, v in eb.items()},
                    "distance_A": item.get("distance"),
                    "donor_h_acceptor_angle_deg": None, "plane_angle_deg": None,
                    "centroid_vector_angle_deg": None,
                    "angle_note": "angle not applicable" if interaction not in {"hydrogen_bond", "weak_hydrogen_bond", "cation_pi", "donor_pi"} else "not yet enriched",
                    "criterion_source": f"PDBe Arpeggio {ARPEGGIO_VERSION} ({ARPEGGIO_COMMIT})",
                })
        elif item.get("type") in {"atom-plane", "plane-plane", "group-plane", "group-group"}:
            labels = item.get("contact", [])
            if isinstance(labels, str):
                labels = [labels]
            for label in labels:
                interaction = TYPE_MAP.get(str(label).upper())
                if interaction:
                    ea, eb = _endpoint(a, elements), _endpoint(b, elements)
                    rows.append({
                        "structure_id": structure_id, "analysis_layer": layer,
                        "interaction_type": interaction,
                        "criterion": _criteria(interaction, a, b),
                        **{f"src_{k}": v for k, v in ea.items()},
                        **{f"dst_{k}": v for k, v in eb.items()},
                        "distance_A": item.get("distance"),
                        "donor_h_acceptor_angle_deg": None, "plane_angle_deg": None,
                        "centroid_vector_angle_deg": None,
                        "angle_note": "angle not applicable" if interaction not in {"hydrogen_bond", "weak_hydrogen_bond", "cation_pi", "donor_pi"} else "not yet enriched",
                        "criterion_source": f"PDBe Arpeggio {ARPEGGIO_VERSION} ({ARPEGGIO_COMMIT})",
                    })
                    continue
                # Arpeggio encodes aromatic stacking as face/edge geometry labels (FT/ET/..).
                if not str(label).upper() in {"FF", "OF", "EE", "FT", "OT", "ET", "FE", "OE", "EF"}:
                    continue
                ea, eb = _endpoint(a, elements), _endpoint(b, elements)
                rows.append({
                    "structure_id": structure_id, "analysis_layer": layer,
                    "interaction_type": "aromatic_stacking",
                    "criterion": f"Arpeggio 1.4.4 ring-plane geometry category {label}; centroid and angular criteria in pinned config",
                    **{f"src_{k}": v for k, v in ea.items()},
                    **{f"dst_{k}": v for k, v in eb.items()},
                    "distance_A": item.get("distance"),
                    "donor_h_acceptor_angle_deg": None, "plane_angle_deg": None,
                    "centroid_vector_angle_deg": None,
                    "angle_note": "plane angle not provided by Arpeggio JSON export",
                    "criterion_source": f"PDBe Arpeggio {ARPEGGIO_VERSION} ({ARPEGGIO_COMMIT})",
                })
    return rows


def _element_lookup(path: Path) -> dict[tuple, dict[str, str]]:
    """Read atom elements from Bio.PDB without changing/choosing altlocs."""
    from Bio.PDB import MMCIFParser, PDBParser
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    structure = parser.get_structure("typed", str(path))
    variants: dict[tuple, set[tuple[str, str]]] = defaultdict(set)
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue.get_unpacked_list():
                    key = (str(chain.id), str(residue.id[1]), str(residue.id[2]).strip(),
                           residue.resname.upper(), atom.name.strip())
                    variants[key].add((str(atom.element or "").upper(), str(atom.get_altloc() or "").strip()))
    result = {}
    for key, atoms in variants.items():
        elements = {item[0] for item in atoms}
        altlocs = sorted({item[1] for item in atoms if item[1]})
        result[key] = {
            "element": sorted(elements)[0] if elements else "",
            "altloc": "" if not altlocs else altlocs[0] if len(altlocs) == 1 else "ambiguous:" + "|".join(altlocs),
        }
    return result


def _model_count(path: Path) -> int:
    from Bio.PDB import MMCIFParser, PDBParser
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    structure = parser.get_structure("typed_model_count", str(path))
    return len(list(structure.get_models()))


def _coordinate_map(path: Path) -> dict[tuple[str, str, str, str, str], list[tuple[str, np.ndarray]]]:
    from Bio.PDB import MMCIFParser, PDBParser
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    structure = parser.get_structure("typed_coords", str(path))
    result: dict[tuple[str, str, str, str, str], list[tuple[str, np.ndarray]]] = defaultdict(list)
    model = next(structure.get_models())
    for chain in model:
        for residue in chain:
            for atom in residue.get_unpacked_list():
                key = (chain.id, str(residue.id[1]), str(residue.id[2]).strip(), residue.resname.upper(), atom.name.strip())
                result[key].append((str(atom.element or "").upper(), np.asarray(atom.coord, dtype=float)))
    return result


def enrich_hbond_angles(evidence: Sequence[dict[str, Any]], coordinate_path: os.PathLike[str] | str) -> None:
    """Measure donor-H-acceptor angles from prepared hydrogen coordinates."""
    path = Path(coordinate_path)
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    cif = MMCIF2Dict(str(path))
    elements = cif.get("_atom_site.type_symbol", [])
    xs, ys, zs = (cif.get(f"_atom_site.Cartn_{axis}", []) for axis in "xyz")
    hydrogen_coords = [np.asarray([float(xs[i]), float(ys[i]), float(zs[i])])
                       for i, element in enumerate(elements) if element.upper() in {"H", "D"}]
    hydrogen_array = np.asarray(hydrogen_coords, dtype=float).reshape((-1, 3))
    coords = None
    for row in evidence:
        kind = row.get("interaction_type")
        if kind in {"cation_pi", "donor_pi", "aromatic_stacking"}:
            heavy_map = row.get("_heavy_coordinate_map")
            if heavy_map is not None:
                def endpoint_coords(prefix: str) -> list[np.ndarray]:
                    key = tuple(_clean(row.get(f"{prefix}_{field}")) for field in ("chain", "resnum", "insertion_code", "resname"))
                    names = [part.strip() for part in _clean(row.get(f"{prefix}_atom")).split(",") if part.strip()]
                    return [xyz for name in names for _, xyz in heavy_map.get((*key, name), [])]
                src_coords, dst_coords = endpoint_coords("src"), endpoint_coords("dst")
                rings = [src_coords, dst_coords]
                if kind == "aromatic_stacking" and len(src_coords) >= 3 and len(dst_coords) >= 3:
                    normal_a, normal_b = _plane_normal(src_coords), _plane_normal(dst_coords)
                    row["plane_angle_deg"] = round(math.degrees(math.acos(float(np.clip(abs(np.dot(normal_a, normal_b)), 0, 1)))), 2)
                    centroid_vector = np.mean(dst_coords, axis=0) - np.mean(src_coords, axis=0)
                    vector_norm = float(np.linalg.norm(centroid_vector))
                    if vector_norm:
                        row["centroid_vector_angle_deg"] = round(math.degrees(math.acos(float(np.clip(abs(np.dot(normal_a, centroid_vector / vector_norm)), 0, 1)))), 2)
                    row["angle_note"] = "ring-normal dihedral and centroid-vector angle measured from original coordinates; Arpeggio JSON reports geometry category only"
                elif kind in {"cation_pi", "donor_pi"}:
                    ring, atom_coords = (src_coords, dst_coords) if len(src_coords) >= 3 else (dst_coords, src_coords)
                    if len(ring) >= 3 and atom_coords:
                        center = np.mean(ring, axis=0)
                        normal = _plane_normal(ring)
                        vector = atom_coords[0] - center
                        vector_norm = float(np.linalg.norm(vector))
                        if vector_norm:
                            row["centroid_vector_angle_deg"] = round(math.degrees(math.acos(float(np.clip(abs(np.dot(normal, vector / vector_norm)), 0, 1)))), 2)
                            row["angle_note"] = "ring-normal theta measured from original coordinates"
            continue
        if kind not in {"hydrogen_bond", "weak_hydrogen_bond"}:
            continue
        if row.get("analysis_layer") == "water_mediated_or_water_endpoint":
            row["angle_note"] = "Arpeggio water hydrogen-bond rule does not require explicit donor-H angle"
            continue
        def endpoint(prefix: str) -> tuple[str, ...]:
            return tuple(_clean(row.get(f"{prefix}_{field}")) for field in ("chain", "resnum", "insertion_code", "resname", "atom"))
        left, right = endpoint("src"), endpoint("dst")
        # Open Babel's hydrogenated mmCIF output drops auth chain/residue identifiers.
        # The donor/acceptor heavy atoms are resolved from the original structure
        # in run_arpeggio and supplied through a temporary attached coordinate map.
        heavy_map = row.get("_heavy_coordinate_map")
        if heavy_map is None:
            row["angle_note"] = "unavailable: Arpeggio hydrogenated export lacks endpoint identifiers"
            continue
        angles = []
        for donor_key, acceptor_key in ((left, right), (right, left)):
            donor_atoms, acceptor_atoms = heavy_map.get(donor_key, []), heavy_map.get(acceptor_key, [])
            if not donor_atoms or not acceptor_atoms:
                continue
            for donor_element, donor_coord in donor_atoms:
                if donor_element not in {"N", "O", "S"}:
                    continue
                near_hydrogens = hydrogen_array[np.linalg.norm(hydrogen_array - donor_coord[None, :], axis=1) <= 1.35]
                for h_coord in near_hydrogens:
                    for _, acceptor_coord in acceptor_atoms:
                        angles.append(angle_degrees(donor_coord, h_coord, acceptor_coord))
        if angles:
            row["donor_h_acceptor_angle_deg"] = round(max(angles), 2)
            row["angle_note"] = "measured from Arpeggio prepared coordinate file"
        else:
            row["angle_note"] = "unavailable: donor H absent or endpoint-to-hydrogen mapping ambiguous"


def _plane_normal(coords: Sequence[np.ndarray]) -> np.ndarray:
    centered = np.asarray(coords, dtype=float) - np.mean(coords, axis=0)
    _, _, vectors = np.linalg.svd(centered, full_matrices=False)
    normal = vectors[-1]
    return normal / np.linalg.norm(normal)


def run_arpeggio(
    structure_path: os.PathLike[str] | str, output_dir: os.PathLike[str] | str,
    *, executable: str = "pdbe-arpeggio", ph: float = 7.4,
    vdw_compensation: float = 0.1, interacting_cutoff_A: float = 5.0,
    write_hydrogenated: bool = True,
    mapped_edges_csv: Optional[os.PathLike[str] | str] = None,
    selections: Optional[Sequence[str]] = None,
) -> dict[str, Any]:
    """Run the optional Arpeggio executable and emit separate evidence tables.

    Requires the pinned ``pdbe-arpeggio==1.4.4`` optional dependency. CNI's
    legacy analysis and files are not invoked or modified by this function.
    """
    structure = Path(structure_path).resolve()
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    model_count = _model_count(structure)
    if model_count != 1:
        raise ValueError(f"Arpeggio JSON omits model IDs; analyze one model per invocation (found {model_count})")
    args = [executable, str(structure), "-m", "-ph", str(ph), "-co", str(vdw_compensation),
            "-i", str(interacting_cutoff_A), "-o", str(out)]
    if selections:
        args.extend(["-s", *selections])
    if write_hydrogenated:
        args.append("-wh")
    completed = subprocess.run(args, check=False, capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(f"Arpeggio failed ({completed.returncode}): {completed.stderr[-2000:]}")
    json_path = out / f"{structure.name.split('.')[0]}.json"
    if not json_path.exists():
        raise FileNotFoundError(f"Arpeggio did not produce expected JSON: {json_path}")
    data = json.loads(json_path.read_text())
    evidence = normalize_arpeggio_contacts(data, structure_id=structure.stem,
                                           element_lookup=_element_lookup(structure))
    stem = structure.name.split('.')[0]
    candidates = [out / f"{stem}_hydrogenated{structure.suffix}", out / f"{stem}_hydrogenated.mmcif", out / f"{stem}_hydrogenated.pdb"]
    hydrogenated = next((candidate for candidate in candidates if candidate.exists()), None)
    if write_hydrogenated and hydrogenated is not None:
        heavy_map = _coordinate_map(structure)
        for row in evidence:
            row["_heavy_coordinate_map"] = heavy_map
        enrich_hbond_angles(evidence, hydrogenated)
        for row in evidence:
            row.pop("_heavy_coordinate_map", None)
    else:
        for row in evidence:
            if row.get("angle_note") in {"not yet enriched", "plane angle not provided by Arpeggio JSON export"}:
                if row.get("interaction_type") in {"hydrogen_bond", "weak_hydrogen_bond", "cation_pi", "donor_pi", "aromatic_stacking"}:
                    row["angle_note"] = "unavailable: prepared coordinate output is disabled or missing"
    # Atom-atom type rows and ring geometry rows are independent; retain all types.
    mapped_edges = []
    if mapped_edges_csv is not None:
        with Path(mapped_edges_csv).open(newline="") as handle:
            mapped_edges = list(csv.DictReader(handle))
        inventory = _residue_inventory(structure)
        for edge in mapped_edges:
            if not edge.get("structure_id") or str(edge.get("structure_id")).lower() == structure.stem.lower():
                _complete_eligibility(edge, inventory)
    _write_csv(out / "atom_interaction_evidence.csv", evidence)
    layer_files = {
        "direct_protein_protein": "direct_protein_protein_interactions.csv",
        "protein_cofactor_or_ligand": "protein_cofactor_interactions.csv",
        "water_mediated_or_water_endpoint": "water_endpoint_interactions.csv",
    }
    for layer, filename in layer_files.items():
        _write_csv(out / filename, [row for row in evidence if row.get("analysis_layer") == layer])
    pairs = summarize_residue_pairs(
        evidence,
        baseline_pairs=[row for row in mapped_edges if str(row.get("baseline_contact", "")).lower() in {"true", "1", "yes"}],
    )
    _write_csv(out / "residue_pair_interactions.csv", pairs)
    structure_summary = summarize_structure_interactions(evidence)
    _write_csv(out / "per_structure_interaction_summary.csv", structure_summary)
    mapped_outputs = []
    if mapped_edges_csv is not None:
        observations = summarize_mapped_edges(evidence, mapped_edges, structure_id=structure.stem)
        _write_csv(out / "mapped_edge_type_observations.csv", observations)
        _write_csv(out / "sequence_balanced_edge_type_summary.csv", sequence_balanced_type_summary(observations))
        mapped_outputs = ["mapped_edge_type_observations.csv", "sequence_balanced_edge_type_summary.csv"]
    manifest = {
        "analysis": "opt-in chemical interaction typing; separate from CNI distance-defined shell contacts",
        "input_structure": str(structure),
        "software": _arpeggio_runtime_versions(executable),
        "settings": {"pH": ph, "vdw_compensation_A": vdw_compensation,
                     "interacting_cutoff_A": interacting_cutoff_A,
                     "include_sequence_adjacent": False, "write_hydrogenated": write_hydrogenated,
                     "use_ambiguities": False, "selections": list(selections or []),
                     "mapped_edges_csv": str(Path(mapped_edges_csv).resolve()) if mapped_edges_csv else None},
        "criteria": {"hydrogen_bond_DA_max_A": 3.9, "hydrogen_bond_DHA_min_deg": 90,
                     "weak_hydrogen_bond_DA_max_A": 3.6, "weak_hydrogen_bond_DHA_min_deg": 130,
                     "weak_hydrogen_bond_CX_angle_deg": [30, 150], "ionic_max_A": 4.0,
                     "hydrophobic_max_A": 4.5, "aromatic_atom_max_A": 4.0,
                     "aromatic_ring_centroid_max_A": 6.0, "aromatic_ring_dihedral_max_deg": 90,
                     "aromatic_ring_centroid_vector_angle_max_deg": 90,
                     "cation_pi_centroid_max_A": 4.5,
                     "cation_pi_theta_max_deg": 30, "methionine_sulfur_pi_max_A": 6.0,
                     "metal_coordination_max_A": 2.8},
        "hydrogen_assumption": "Arpeggio/Open Babel adds hydrogens at configured pH only if no H/D exists; any input H/D suppresses addition for the whole structure.",
        "alternate_locations": "Arpeggio evaluates all observed alternate-location conformers; missing atoms are not imputed.",
        "model_count": model_count,
        "histidine_protonation": "No histidine tautomer is forced by this adapter; Open Babel uses the requested pH only when it adds hydrogens. Existing H/D coordinates are trusted as supplied.",
        "outputs": [json_path.name, "atom_interaction_evidence.csv", *layer_files.values(),
                    "residue_pair_interactions.csv", "per_structure_interaction_summary.csv", *mapped_outputs],
        "hydrogenated_coordinate_file": hydrogenated.name if hydrogenated is not None else None,
        "atom_evidence_rows": len(evidence),
        "stdout": completed.stdout[-1000:], "stderr": completed.stderr[-1000:],
    }
    (out / "interaction_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def run_arpeggio_batch(
    structure_paths: Iterable[os.PathLike[str] | str], output_dir: os.PathLike[str] | str,
    *, mapped_edges_csv: Optional[os.PathLike[str] | str] = None,
    executable: str = "pdbe-arpeggio", ph: float = 7.4,
    vdw_compensation: float = 0.1, interacting_cutoff_A: float = 5.0,
    selections_by_structure: Optional[Mapping[str, Sequence[str]]] = None,
) -> dict[str, Any]:
    """Run structures independently and combine mapped edge sequence votes."""
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    manifests = []
    observations = []
    for path in structure_paths:
        item = Path(path)
        manifest = run_arpeggio(
            item, out / item.stem, executable=executable, ph=ph,
            vdw_compensation=vdw_compensation, interacting_cutoff_A=interacting_cutoff_A,
            mapped_edges_csv=mapped_edges_csv,
            selections=(selections_by_structure or {}).get(item.stem),
        )
        manifests.append({"structure_id": item.stem, "manifest": manifest})
        mapped_path = out / item.stem / "mapped_edge_type_observations.csv"
        if mapped_path.exists():
            with mapped_path.open(newline="") as handle:
                for row in csv.DictReader(handle):
                    try:
                        row["eligible_by_type"] = json.loads(row.get("eligible_by_type") or "{}")
                    except json.JSONDecodeError:
                        row["eligible_by_type"] = {}
                    row["eligible"] = str(row.get("eligible", "")).lower() in {"true", "1", "yes"}
                    observations.append(row)
    if observations:
        _write_csv(out / "sequence_balanced_edge_type_summary.csv", sequence_balanced_type_summary(observations))
    batch_manifest = {
        "analysis": "typed interactions; per-structure Arpeggio results and sequence-balanced mapped-edge summary",
        "software": _arpeggio_runtime_versions(executable),
        "settings": {"pH": ph, "vdw_compensation_A": vdw_compensation,
                     "interacting_cutoff_A": interacting_cutoff_A,
                     "mapped_edges_csv": str(Path(mapped_edges_csv).resolve()) if mapped_edges_csv else None},
        "structures": manifests,
        "mapped_edge_observations": len(observations),
        "outputs": ["<structure_id>/atom_interaction_evidence.csv", "<structure_id>/residue_pair_interactions.csv"] + (["sequence_balanced_edge_type_summary.csv"] if observations else []),
    }
    (out / "typed_interaction_batch_manifest.json").write_text(json.dumps(batch_manifest, indent=2) + "\n")
    return batch_manifest


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    columns = EVIDENCE_COLUMNS if path.name == "atom_interaction_evidence.csv" or path.name.endswith("_interactions.csv") else PAIR_COLUMNS if path.name == "residue_pair_interactions.csv" else None
    if path.name == "per_structure_interaction_summary.csv":
        columns = ["structure_id", "analysis_layer", "interaction_type", "atom_evidence_rows", "residue_pair_count"]
    with path.open("w", newline="") as handle:
        if columns is None:
            columns = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, list, tuple)) else value
                             for key, value in row.items()})


def _arpeggio_runtime_versions(executable: str) -> dict[str, str]:
    executable_path = shutil.which(executable) or executable
    python = str(Path(executable_path).resolve().parent / "python")
    probe = (
        "import arpeggio, gemmi; from openbabel import openbabel as ob; "
        "import sys; print(arpeggio.__version__+'|'+gemmi.__version__+'|'+ob.OBReleaseVersion()+'|'+sys.version.split()[0])"
    )
    try:
        result = subprocess.run([python, "-c", probe], capture_output=True, text=True, check=True, timeout=10)
        arpeggio, gemmi, openbabel, runtime_python = result.stdout.strip().split("|", 3)
    except (OSError, subprocess.SubprocessError, ValueError):
        arpeggio, gemmi, openbabel = ARPEGGIO_VERSION, "unknown", "unknown"
        try:
            version = subprocess.run([python, "--version"], capture_output=True, text=True, check=True, timeout=10)
            runtime_python = (version.stdout or version.stderr).strip().removeprefix("Python ")
        except (OSError, subprocess.SubprocessError):
            runtime_python = "unknown"
    return {"pdbe_arpeggio": arpeggio, "arpeggio_commit": ARPEGGIO_COMMIT,
            "gemmi": gemmi, "openbabel": openbabel, "python": runtime_python}


def summarize_structure_interactions(evidence: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    counts: dict[tuple[str, str], dict[str, Any]] = {}
    pairs: dict[tuple[str, str], set[tuple[str, ...]]] = defaultdict(set)
    for row in evidence:
        key = (str(row.get("analysis_layer", "")), str(row.get("interaction_type", "")))
        record = counts.setdefault(key, {
            "structure_id": str(row.get("structure_id", "")), "analysis_layer": key[0],
            "interaction_type": key[1], "atom_evidence_rows": 0,
        })
        record["atom_evidence_rows"] += 1
        pairs[key].add(_pair_key(row))
    for key, record in counts.items():
        record["residue_pair_count"] = len(pairs[key])
    return [counts[key] for key in sorted(counts)]


def _pair_key(row: Mapping[str, Any]) -> tuple[str, ...]:
    a = tuple(_clean(row.get(f"src_{x}")) for x in ("resname", "chain", "resnum", "insertion_code"))
    b = tuple(_clean(row.get(f"dst_{x}")) for x in ("resname", "chain", "resnum", "insertion_code"))
    return tuple(sorted((a, b)))


def summarize_residue_pairs(
    evidence: Iterable[Mapping[str, Any]], *, baseline_pairs: Iterable[Mapping[str, Any]] = (),
) -> list[dict[str, Any]]:
    """Collapse atom-level evidence by residue pair; optionally include baseline-only pairs."""
    grouped: dict[tuple, dict[str, Any]] = {}
    for row in evidence:
        key = _pair_key(row)
        current = grouped.setdefault(key, {"row": row, "types": set(), "layers": set()})
        current["types"].add(str(row["interaction_type"]))
        current["layers"].add(str(row["analysis_layer"]))
    for item in baseline_pairs:
        key = _pair_key(item)
        baseline_item = dict(item)
        baseline_item.setdefault("baseline_contact", True)
        grouped.setdefault(key, {"row": baseline_item, "types": set(), "layers": {"direct_protein_protein"}})
    result = []
    for key, value in sorted(grouped.items()):
        row = value["row"]
        types = sorted(value["types"])
        a, b = key
        layer = sorted(value["layers"])[0] if value["layers"] else "direct_protein_protein"
        baseline_value = row.get("baseline_contact")
        baseline = None if baseline_value is None or str(baseline_value).strip() == "" else str(baseline_value).lower() in {"true", "1", "yes"}
        result.append({
            "structure_id": _clean(row.get("structure_id")), "analysis_layer": layer,
            "src_resname": a[0], "src_chain": a[1], "src_resnum": a[2], "src_insertion_code": a[3],
            "dst_resname": b[0], "dst_chain": b[1], "dst_resnum": b[2], "dst_insertion_code": b[3],
            "baseline_contact": baseline,
            "interaction_types": ";".join(types),
            "outcome": ";".join(types) if types else "proximity_only_unclassified" if baseline is True else "no_typed_interaction",
        })
    return result


def sequence_balanced_type_summary(observations: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Summarize per-type support with explicit assessable denominators.

    Each row describes one structure × residue edge, with ``sequence_hash``,
    ``group``, ``eligible`` and ``interaction_types`` (iterable or semicolon text).
    Replicate structures cast a sequence-level positive vote at >= half of
    eligible replicates, matching the manuscript's existing rule.
    """
    records = []
    for raw in observations:
        row = dict(raw)
        value = row.get("interaction_types", [])
        if isinstance(value, str):
            types = {part for part in value.split(";") if part}
        else:
            types = set(value or [])
        row["types"] = types
        records.append(row)
    groups = sorted({str(row.get("group", "")) for row in records})
    edges = sorted({str(row.get("edge_id", row.get("edge", "all_edges"))) for row in records})
    types = sorted(ATOM_TYPES | {t for row in records for t in row["types"]})
    result = []
    for edge in edges:
      for group in groups:
        group_rows = [r for r in records if str(r.get("group", "")) == group and str(r.get("edge_id", r.get("edge", "all_edges"))) == edge]
        sequences = sorted({str(r.get("sequence_hash", "")) for r in group_rows if r.get("sequence_hash")})
        for kind in types:
            eligible_sequences = 0
            positive_sequences = 0
            eligible_structures = 0
            positive_structures = 0
            for seq in sequences:
                reps = [r for r in group_rows if str(r.get("sequence_hash")) == seq
                        and _type_eligible(r, kind)]
                if not reps:
                    continue
                eligible_sequences += 1
                eligible_structures += len(reps)
                positives = sum(kind in r["types"] for r in reps)
                positive_structures += positives
                if positives >= math.ceil(len(reps) / 2):
                    positive_sequences += 1
            result.append({
                "edge_id": edge, "group": group, "interaction_type": kind,
                "eligible_sequences": eligible_sequences,
                "positive_sequences": positive_sequences,
                "sequence_support_fraction": positive_sequences / eligible_sequences if eligible_sequences else None,
                "eligible_structures": eligible_structures, "positive_structures": positive_structures,
                "replicate_rule": "sequence positive when >= half eligible structures are positive",
            })
    return result


def summarize_mapped_edges(
    evidence: Iterable[Mapping[str, Any]], mapped_edges: Iterable[Mapping[str, Any]], *, structure_id: str,
) -> list[dict[str, Any]]:
    """Join a user's mapped edge ledger to atom evidence for this structure.

    The CSV requires `edge_id`, `group`, `sequence_hash`, `src_*`, `dst_*`,
    `eligible`, and may carry JSON in `eligible_by_type`. Endpoint columns use
    auth chain/number/insertion and residue names.
    """
    evidence_by_pair: dict[tuple[str, ...], set[str]] = defaultdict(set)
    for row in evidence:
        evidence_by_pair[_pair_key(row)].add(str(row.get("interaction_type", "")))
    result = []
    for raw in mapped_edges:
        row = dict(raw)
        if row.get("structure_id") and str(row["structure_id"]).lower() != structure_id.lower():
            continue
        try:
            eligible_by_type = json.loads(row.get("eligible_by_type") or "{}")
        except (TypeError, json.JSONDecodeError):
            eligible_by_type = {}
        key = _pair_key(row)
        kinds = sorted(evidence_by_pair.get(key, set()))
        baseline = str(row.get("baseline_contact", "")).lower() in {"true", "1", "yes"}
        result.append({
            **row, "structure_id": structure_id, "eligible_by_type": eligible_by_type,
            "interaction_types": ";".join(kinds),
            "outcome": ";".join(kinds) if kinds else "proximity_only_unclassified" if baseline else "no_baseline_contact_or_typed_evidence",
        })
    return result


def _type_eligible(row: Mapping[str, Any], interaction_type: str) -> bool:
    """Use type-specific assessability when provided; otherwise generic eligibility."""
    by_type = row.get("eligible_by_type", {})
    if isinstance(by_type, Mapping) and interaction_type in by_type:
        value = by_type[interaction_type]
    else:
        value = row.get("eligible")
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)

