"""Adapters from legacy SSCNA CSV reports to the canonical table contract.

Legacy reports are intentionally left untouched.  This module makes their
information consumable by the current comparison layer and records where the
legacy semantics are weaker: a legacy ``cofactor->pcs`` link is a geometric
primary-shell contact, not proof of validated direct coordination.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, Union

import pandas as pd

from .coordination_api import (
    ATOM_COLUMNS,
    CONTACT_COLUMNS,
    LINK_COLUMNS,
    RESIDUE_COLUMNS,
)
from .motif_registry import motif_for_atom


PathLike = Union[str, Path]


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if value != value:
            return ""
    except Exception:
        pass
    return str(value).strip()


def _float_or_none(value: Any) -> Optional[float]:
    try:
        if value in (None, ""):
            return None
        number = float(value)
        return number if number == number else None
    except (TypeError, ValueError):
        return None


def _structure_id_from_path(path: Path) -> str:
    name = path.name
    for suffix in ("_Coord_Breakdown_atoms.csv", "_Coord_Breakdown.csv", "_Coord_Links.csv"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for suffix in (".mmcif", ".cif", ".pdb"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def _read_csv_rows(path: PathLike) -> List[Dict[str, str]]:
    with Path(path).expanduser().open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _read_breakdown_sections(path: PathLike) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    """Read the summary and optional atom sections from legacy breakdown CSV."""
    breakdown_path = Path(path).expanduser()
    with breakdown_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))

    atom_header_index = next(
        (
            index
            for index, row in enumerate(rows)
            if row[:5] == ["Residue Category", "Residue Name", "Residue Number", "Chain", "Atom"]
        ),
        None,
    )
    def _serialize(csv_rows: Iterable[List[str]]) -> str:
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerows(row for row in csv_rows if row)
        return output.getvalue()

    summary_rows = rows[:atom_header_index] if atom_header_index is not None else rows
    summary_text = io.StringIO(_serialize(summary_rows))
    summaries = [dict(row) for row in csv.DictReader(summary_text)] if summary_text.getvalue() else []

    atoms: List[Dict[str, str]] = []
    if atom_header_index is not None:
        atom_text = io.StringIO(_serialize(rows[atom_header_index:]))
        atoms = [dict(row) for row in csv.DictReader(atom_text)] if atom_text.getvalue() else []
    return summaries, atoms


def _legacy_motif(residue: str, atom: str, moiety: Any = "") -> str:
    explicit = _text(moiety).lower()
    if explicit and explicit not in {"unknown_moiety", "unknown_motif"}:
        return explicit
    return motif_for_atom(residue, atom)


def _atom_key(category: Any, residue: Any, number: Any, chain: Any, atom: Any) -> Tuple[str, str, str, str, str]:
    return (
        _text(category).upper(),
        _text(residue).upper(),
        _text(number),
        _text(chain),
        _text(atom).upper(),
    )


def _residue_key(category: Any, residue: Any, number: Any, chain: Any) -> Tuple[str, str, str, str]:
    return (
        _text(category).upper(),
        _text(residue).upper(),
        _text(number),
        _text(chain),
    )


def _empty_analysis_tables() -> Dict[str, pd.DataFrame]:
    return {
        "atoms": pd.DataFrame(columns=ATOM_COLUMNS),
        "residues": pd.DataFrame(columns=RESIDUE_COLUMNS),
        "links": pd.DataFrame(columns=LINK_COLUMNS),
        "contacts": pd.DataFrame(columns=CONTACT_COLUMNS),
    }


def legacy_csvs_to_analysis(
    coord_breakdown_path: PathLike,
    coord_links_path: Optional[PathLike] = None,
    *,
    coord_breakdown_atoms_path: Optional[PathLike] = None,
    structure_id: Optional[str] = None,
    site_id: str = "all",
) -> Dict[str, pd.DataFrame]:
    """Map legacy CSV artifacts to the current analysis-table contract.

    ``coord_breakdown_path`` may contain both legacy sections separated by a
    blank line. If the atom section is absent, provide the separate
    ``coord_breakdown_atoms_path`` emitted by the cohesive legacy pipeline.
    Coordinates and atom-level identity are preserved when available; legacy
    reports do not contain model IDs, insertion codes, hetero flags, or
    element columns, so those fields are left empty rather than inferred.
    """
    breakdown_path = Path(coord_breakdown_path).expanduser()
    summaries, embedded_atoms = _read_breakdown_sections(breakdown_path)
    if coord_breakdown_atoms_path is not None:
        atom_rows = _read_csv_rows(coord_breakdown_atoms_path)
    else:
        atom_rows = embedded_atoms
    link_rows = _read_csv_rows(coord_links_path) if coord_links_path is not None else []

    resolved_structure_id = structure_id or _structure_id_from_path(breakdown_path)
    tables = _empty_analysis_tables()

    atom_records: List[Dict[str, Any]] = []
    motif_by_atom: Dict[Tuple[str, str, str, str, str], str] = {}
    residue_records: Dict[Tuple[str, str, str, str], Dict[str, Any]] = {}
    for row in atom_rows:
        category = _text(row.get("Category", row.get("Residue Category")))
        residue = _text(row.get("Residue", row.get("Residue Name"))).upper()
        number = _text(row.get("Residue Number"))
        chain = _text(row.get("Chain"))
        atom = _text(row.get("Atom", row.get("atom_name"))).upper()
        motif = _legacy_motif(residue, atom, row.get("Moiety", row.get("motif")))
        atom_key = _atom_key(category, residue, number, chain, atom)
        motif_by_atom[atom_key] = motif
        interactor = _text(row.get("InteractorFlag", row.get("interactor_flag"))).lower()
        category_upper = category.upper()
        role = "active_site_component"
        if interactor == "interactor" and category_upper == "PCS":
            role = "primary_coordinator"
        elif interactor == "interactor" and category_upper in {"SCS", "TCS", "SHELL3"}:
            role = "network_context"

        record = {
            "structure_id": resolved_structure_id,
            "site_id": site_id,
            "shell": "Cofactor" if category_upper == "COFACTOR" else category_upper,
            "residue_name": residue,
            "residue_number": number,
            "chain": chain,
            "insertion_code": "",
            "hetero_flag": "",
            "model_id": "",
            "atom_name": atom,
            "element": "",
            "motif": motif,
            "coordination_role": role,
            "x": _float_or_none(row.get("x")),
            "y": _float_or_none(row.get("y")),
            "z": _float_or_none(row.get("z")),
        }
        atom_records.append(record)

        residue_key = _residue_key(category, residue, number, chain)
        residue_record = residue_records.setdefault(
            residue_key,
            {
                "shell": record["shell"],
                "residue_name": residue,
                "residue_number": number,
                "chain": chain,
                "atoms": set(),
                "motifs": set(),
                "roles": set(),
            },
        )
        residue_record["atoms"].add(atom)
        residue_record["motifs"].add(motif)
        residue_record["roles"].add(role)

    cofactor_records = [record for record in atom_records if record["shell"] == "Cofactor"]
    cofactor_names = ";".join(dict.fromkeys(record["residue_name"] for record in cofactor_records))
    cofactor_numbers = ";".join(dict.fromkeys(record["residue_number"] for record in cofactor_records))
    cofactor_chains = ";".join(dict.fromkeys(record["chain"] for record in cofactor_records))

    distance_by_residue: Dict[Tuple[str, str, str], float] = {}
    canonical_links: List[Dict[str, Any]] = []
    canonical_contacts: List[Dict[str, Any]] = []
    for row in link_rows:
        link_type = _text(row.get("link_type")).lower()
        source_category = "Cofactor" if link_type == "cofactor->pcs" else "PCS" if link_type == "pcs->scs" else "SCS"
        destination_category = "PCS" if link_type == "cofactor->pcs" else "SCS" if link_type == "pcs->scs" else "TCS"
        source_residue = _text(row.get("src_resname")).upper()
        source_number = _text(row.get("src_resnum"))
        source_chain = _text(row.get("src_chain"))
        source_atom = _text(row.get("src_atom")).upper()
        destination_residue = _text(row.get("dst_resname")).upper()
        destination_number = _text(row.get("dst_resnum"))
        destination_chain = _text(row.get("dst_chain"))
        destination_atom = _text(row.get("dst_atom")).upper()
        distance = _float_or_none(row.get("distance_A"))
        source_motif = motif_by_atom.get(
            _atom_key(source_category, source_residue, source_number, source_chain, source_atom),
            _legacy_motif(source_residue, source_atom, row.get("src_moiety")),
        )
        destination_motif = motif_by_atom.get(
            _atom_key(destination_category, destination_residue, destination_number, destination_chain, destination_atom),
            _legacy_motif(destination_residue, destination_atom, row.get("dst_moiety")),
        )
        canonical_links.append(
            {
                "structure_id": resolved_structure_id,
                "link_type": link_type,
                "src_resname": source_residue,
                "src_resnum": source_number,
                "src_chain": source_chain,
                "src_atom": source_atom,
                "src_moiety": source_motif,
                "dst_resname": destination_residue,
                "dst_resnum": destination_number,
                "dst_chain": destination_chain,
                "dst_atom": destination_atom,
                "dst_moiety": destination_motif,
                "distance_A": distance,
                "src_insertion_code": "",
                "src_hetero_flag": "",
                "dst_insertion_code": "",
                "dst_hetero_flag": "",
                "src_element": "",
                "dst_element": "",
                "direct_coordination": False,
            }
        )
        if distance is not None:
            key = (destination_residue, destination_number, destination_chain)
            previous = distance_by_residue.get(key)
            if previous is None or distance < previous:
                distance_by_residue[key] = distance

        if link_type == "cofactor->pcs":
            contact_role = "legacy_primary_contact"
        elif link_type == "pcs->scs":
            contact_role = "legacy_secondary_contact"
        else:
            contact_role = "legacy_tertiary_contact"
        canonical_contacts.append(
            {
                "structure_id": resolved_structure_id,
                "site_id": site_id,
                "src_shell": source_category,
                "dst_shell": destination_category,
                "link_type": link_type,
                "contact_role": contact_role,
                "direct_coordination": False,
                "src_resname": source_residue,
                "src_resnum": source_number,
                "src_chain": source_chain,
                "src_insertion_code": "",
                "src_hetero_flag": "",
                "src_model_id": "",
                "src_atom": source_atom,
                "src_element": "",
                "src_motif": source_motif,
                "dst_resname": destination_residue,
                "dst_resnum": destination_number,
                "dst_chain": destination_chain,
                "dst_insertion_code": "",
                "dst_hetero_flag": "",
                "dst_model_id": "",
                "dst_atom": destination_atom,
                "dst_element": "",
                "dst_motif": destination_motif,
                "distance_A": distance,
            }
        )

    residue_rows: List[Dict[str, Any]] = []
    for residue_key, record in residue_records.items():
        category, residue, number, chain = residue_key
        residue_rows.append(
            {
                "structure_id": resolved_structure_id,
                "site_id": site_id,
                "cofactor_residue_name": cofactor_names,
                "cofactor_residue_number": cofactor_numbers,
                "cofactor_chain": cofactor_chains,
                "cofactor_insertion_code": "",
                "shell": record["shell"],
                "residue_name": residue,
                "residue_number": number,
                "chain": chain,
                "insertion_code": "",
                "hetero_flag": "",
                "atoms_involved": ",".join(sorted(record["atoms"])),
                "minimum_distance_A": None
                if category == "COFACTOR"
                else distance_by_residue.get((residue, number, chain)),
            }
        )

    tables["atoms"] = pd.DataFrame(atom_records, columns=ATOM_COLUMNS)
    tables["links"] = pd.DataFrame(canonical_links, columns=LINK_COLUMNS)
    tables["contacts"] = pd.DataFrame(canonical_contacts, columns=CONTACT_COLUMNS)
    tables["residues"] = pd.DataFrame(residue_rows, columns=RESIDUE_COLUMNS)
    return tables


def legacy_csvs_to_signature(
    coord_breakdown_path: PathLike,
    coord_links_path: Optional[PathLike] = None,
    *,
    coord_breakdown_atoms_path: Optional[PathLike] = None,
    structure_id: Optional[str] = None,
    site_id: str = "all",
    residue_map: Optional[Mapping[Any, Any]] = None,
) -> Dict[str, Any]:
    """Build a canonical comparison signature from legacy CSV artifacts."""
    from .network_comparison import build_network_signature

    tables = legacy_csvs_to_analysis(
        coord_breakdown_path,
        coord_links_path,
        coord_breakdown_atoms_path=coord_breakdown_atoms_path,
        structure_id=structure_id,
        site_id=site_id,
    )
    signature = build_network_signature(
        tables,
        structure_id=structure_id or _structure_id_from_path(Path(coord_breakdown_path)),
        site_id=site_id,
        residue_map=residue_map,
    )
    signature["source"] = {
        "format": "legacy_csv",
        "coord_breakdown": str(Path(coord_breakdown_path).expanduser()),
        "coord_links": str(Path(coord_links_path).expanduser()) if coord_links_path else None,
    }
    return signature


__all__ = ["legacy_csvs_to_analysis", "legacy_csvs_to_signature"]
