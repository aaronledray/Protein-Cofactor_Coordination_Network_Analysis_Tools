"""Evaluation plots and the cohesive interactive coordination-network viewer."""

import csv
import html
import json
import logging
import re
from collections import OrderedDict
from typing import Any, Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.colors import sample_colorscale

from .display import show_matplotlib
from .motif_registry import motif_for_atom
from .structure_processing import generate_residue_bonds

logger = logging.getLogger(__name__)

__all__ = ["plot_evaluation_results", "plot_interactive_cohesive_network"]


def plot_evaluation_results(
    df_results: pd.DataFrame,
    mode: str = "CA_only",
    output_csv: Optional[str] = None,
    highlight_labels: Optional[List[str]] = None,
) -> None:
    """
    Plot a horizontal bar chart based on the evaluation mode and (optionally) save
    a CSV including matched residue numbers.

    Parameters
    ----------
    df_results : pandas.DataFrame
        Must contain:
          - "Query_File"
          - "Matches" (if mode == "CA_only")
          - "Score"   (if mode == "CA_CB_vectors")
        (Any extra columns like "Matched_Residues" are preserved for CSV output.)
    mode : {"CA_only", "CA_CB_vectors"}
        Chooses what to show on the x-axis.
    output_csv : str or None
        If provided, write `df_results` to this CSV path.
    highlight_labels : list of substrings or None
        Any y-tick label containing one of these substrings gets a yellow background.
    """
    if highlight_labels is None:
        highlight_labels = []

    # Select the series to plot
    if mode == "CA_only":
        if "Matches" not in df_results.columns:
            raise ValueError("df_results must have a 'Matches' column for mode='CA_only'.")
        values = df_results["Matches"]
        xlabel = "Number of Matches (CA)"
        title = "CA-only Matches for Each Query File"
    elif mode == "CA_CB_vectors":
        if "Score" not in df_results.columns:
            raise ValueError("df_results must have a 'Score' column for mode='CA_CB_vectors'.")
        values = df_results["Score"]
        xlabel = "CA→CB Vector Score"
        title = "CA→CB Vector Scores for Each Query File"
    else:
        raise ValueError("Unknown mode. Use 'CA_only' or 'CA_CB_vectors'.")

    labels = df_results["Query_File"]

    # Plot
    plt.figure(figsize=(10, 6))
    bars = plt.barh(labels, values, color="grey")
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel("Query File", fontsize=12)
    plt.title(title, fontsize=14)
    plt.tight_layout()

    # Highlight requested labels
    ax = plt.gca()
    for label in ax.get_yticklabels():
        text = label.get_text()
        if any(substr in text for substr in highlight_labels):
            label.set_bbox({"facecolor": "yellow", "edgecolor": "none", "pad": 2})

    show_matplotlib(plt)

    # Optional CSV
    if output_csv:
        df_results.to_csv(output_csv, index=False)
        logger.info("Results written to %s", output_csv)



















# UPDATED FOR COORD DOTTED LINES:







def plot_interactive_cohesive_network(
    structure,
    cofactor_atoms: Sequence[Dict],
    pcs_atoms: Sequence[Dict],
    scs_atoms: Sequence[Dict],
    bond_lookup_table: Dict[str, List[List[str]]],
    *,
    contacts: Optional[Any] = None,
    focused_atoms_by_shell: Optional[Dict[str, Sequence[Dict]]] = None,
    pdb_name: str = "structure.pdb",
    cofactor_resname: str = "cofactor",
    atom_type_colors: Optional[Dict[str, str]] = None,
    output_filename: str = "coordination_network.html",
    first_model_only: bool = False,
    compact_html: bool = False,
    conservation_by_residue: Optional[Dict[tuple, float]] = None,
    conservation_label: str = "Profile conservation",
    external_residue_network_csv: Optional[str] = None,
) -> None:
    """Write a layered, standalone interactive coordination-network viewer.

    This viewer combines full-protein context, focused shell atoms, motif
    coloring, and exact motif-contact overlays. The legacy viewer function
    above remains available for compatibility, but is not used by the main
    SSCNA CLI output path.
    """
    if atom_type_colors is None:
        atom_type_colors = {
            "C": "#444444", "N": "#2166ac", "O": "#d73027", "S": "#fdae61",
            "FE": "#e66101", "MN": "#762a83", "CA": "#1b7837", "CU": "#b8860b",
            "ZN": "#5e3c99", "MO": "#018571",
        }

    backbone_names = {"N", "H", "CA", "HA", "C", "O", "OXT"}

    def _value(atom: Dict, *names: str, default: Any = "") -> Any:
        for name in names:
            if name in atom and atom[name] is not None:
                return atom[name]
        return default

    def _motif_of(atom: Dict) -> str:
        residue_name = str(_value(atom, "residue", "residue_name", "src_resname", "dst_resname")).upper()
        atom_name = str(_value(atom, "name", "atom_name", "src_atom", "dst_atom")).upper()
        label = motif_for_atom(residue_name, atom_name)
        if label == "unknown_motif" and atom_name in backbone_names:
            label = "backbone"
        return label

    def _normalize_atom(atom: Dict, shell: str) -> Dict:
        coordinates = atom.get("coordinates")
        if coordinates is None:
            coordinates = [atom.get("x"), atom.get("y"), atom.get("z")]
        return {
            "name": _value(atom, "name", "atom_name"),
            "element": _value(atom, "element"),
            "residue": _value(atom, "residue", "residue_name"),
            "residue_number": _value(atom, "residue_number"),
            "chain": _value(atom, "chain"),
            "insertion_code": _value(atom, "insertion_code"),
            "hetero_flag": _value(atom, "hetero_flag"),
            "model_id": _value(atom, "model_id"),
            "motif": _value(atom, "motif", default=_motif_of(atom)),
            "coordination_role": _value(
                atom,
                "coordination_role",
                default="active_site_component" if shell != "Protein" else "full_protein",
            ),
            "shell": shell,
            "coordinates": np.asarray(coordinates, dtype=float),
        }

    def _structure_atoms() -> List[Dict]:
        models = list(structure)
        if first_model_only:
            models = models[:1]
        atoms: List[Dict] = []
        for model in models:
            for chain in model:
                for residue in chain:
                    residue_id = residue.get_id()
                    for atom in residue:
                        atoms.append(_normalize_atom({
                            "name": atom.get_name(),
                            "element": getattr(atom, "element", ""),
                            "residue": residue.get_resname(),
                            "residue_number": residue_id[1],
                            "chain": chain.id,
                            "insertion_code": str(residue_id[2] or "").strip(),
                            "hetero_flag": str(residue_id[0] or "").strip(),
                            "model_id": model.id,
                            "coordinates": np.asarray(atom.coord, dtype=float),
                        }, "Protein"))
        return atoms

    if focused_atoms_by_shell is None:
        focused_atoms_by_shell = {
            "Cofactor": cofactor_atoms,
            "PCS": pcs_atoms,
            "SCS": scs_atoms,
        }
    shell_groups = {
        str(shell): [_normalize_atom(atom, str(shell)) for atom in atoms]
        for shell, atoms in focused_atoms_by_shell.items()
    }
    focused_atoms = [atom for atoms in shell_groups.values() for atom in atoms]
    full_atoms = _structure_atoms()

    def _residue_key(atom: Dict, include_model: bool = True) -> tuple:
        return (
            str(atom.get("model_id", "")) if include_model else "",
            str(atom.get("residue", "")), str(atom.get("residue_number", "")),
            str(atom.get("chain", "")), str(atom.get("insertion_code", "")),
            str(atom.get("hetero_flag", "")),
        )

    def _conservation_value(atom: Dict) -> float:
        if conservation_by_residue is None:
            return 0.0
        key = _residue_key(atom)
        value = conservation_by_residue.get(key)
        if value is None:
            value = conservation_by_residue.get(_residue_key(atom, include_model=False), 0.0)
        try:
            return max(0.0, min(1.0, float(value)))
        except (TypeError, ValueError):
            return 0.0

    focused_residue_keys = {_residue_key(atom) for atom in focused_atoms}
    motif_atoms = [
        atom for atom in full_atoms
        if _residue_key(atom) in focused_residue_keys
        and str(atom["motif"]) not in {"unknown_motif", "backbone"}
    ]
    backbone_atoms = [
        atom for atom in full_atoms
        if _residue_key(atom) in focused_residue_keys
        and str(atom["motif"]) == "backbone"
    ]

    def _atom_key(atom: Dict, include_model: bool = True) -> tuple:
        return (
            str(atom.get("model_id", "")) if include_model else "",
            str(atom.get("residue", "")), str(atom.get("residue_number", "")),
            str(atom.get("chain", "")), str(atom.get("insertion_code", "")),
            str(atom.get("hetero_flag", "")), str(atom.get("name", "")),
        )

    coordinates_by_key = {_atom_key(atom): atom["coordinates"] for atom in full_atoms}
    coordinates_by_base_key = {
        _atom_key(atom, include_model=False): atom["coordinates"] for atom in full_atoms
    }

    external_graph_groups: OrderedDict = OrderedDict()
    external_graph_unmapped: List[str] = []
    external_graph_input_count = 0
    if external_residue_network_csv:
        aliases = {
            "source_chain": ("source_chain", "source_chain_id", "src_chain"),
            "source_number": ("source_residue_number", "source_resnum", "src_resnum", "source_number"),
            "target_chain": ("target_chain", "target_chain_id", "dst_chain"),
            "target_number": ("target_residue_number", "target_resnum", "dst_resnum", "target_number"),
            "source_icode": ("source_insertion_code", "source_icode", "src_insertion_code"),
            "target_icode": ("target_insertion_code", "target_icode", "dst_insertion_code"),
        }

        def _residue_number(value: str) -> str:
            # CSVs saved by spreadsheets commonly serialize integer residue
            # numbers as e.g. "42.0". Normalize those to the PDB integer ID.
            try:
                number = float(value)
                if number.is_integer():
                    return str(int(number))
            except (TypeError, ValueError):
                pass
            return str(value).strip()

        structure_residues: Dict[tuple, Dict[str, Any]] = {}
        representative_models = list(structure)[:1]
        for model in representative_models:
            for chain in model:
                for residue in chain:
                    residue_id = residue.get_id()
                    ca = residue.child_dict.get("CA")
                    key = (str(chain.id).strip(), str(residue_id[1]), str(residue_id[2] or "").strip())
                    if ca is not None:
                        structure_residues[key] = {
                            "coordinates": np.asarray(ca.coord, dtype=float),
                            "resname": str(residue.get_resname()).strip(),
                            "model_id": str(model.id),
                        }

        with open(external_residue_network_csv, newline="", encoding="utf-8-sig") as csv_file:
            reader = csv.DictReader(csv_file)
            fieldnames = list(reader.fieldnames or [])
            normalized_fields = {str(field).strip().lower(): field for field in fieldnames}
            required = ("source_chain", "source_number", "target_chain", "target_number")
            missing = [
                name for name in required
                if not any(alias in normalized_fields for alias in aliases[name])
            ]
            if missing:
                raise ValueError(
                    "External residue-network CSV requires source_chain, "
                    "source_residue_number, target_chain, and target_residue_number columns"
                )

            def _row_value(row: Dict[str, Any], logical_name: str) -> str:
                actual = next(
                    (normalized_fields[alias] for alias in aliases[logical_name] if alias in normalized_fields),
                    "",
                )
                return str(row.get(actual, "") or "").strip()

            for row_index, row in enumerate(reader, start=2):
                external_graph_input_count += 1
                source_chain = _row_value(row, "source_chain")
                source_number = _residue_number(_row_value(row, "source_number"))
                target_chain = _row_value(row, "target_chain")
                target_number = _residue_number(_row_value(row, "target_number"))
                source_icode = _row_value(row, "source_icode")
                target_icode = _row_value(row, "target_icode")
                source_label = f"{source_chain}:{source_number}{source_icode}"
                target_label = f"{target_chain}:{target_number}{target_icode}"
                source = structure_residues.get((source_chain, source_number, source_icode))
                target = structure_residues.get((target_chain, target_number, target_icode))
                if source is None or target is None:
                    missing_ends = []
                    if source is None:
                        missing_ends.append(f"source {source_label}")
                    if target is None:
                        missing_ends.append(f"target {target_label}")
                    external_graph_unmapped.append(
                        f"CSV row {row_index}: {source_label} → {target_label} "
                        f"({', '.join(missing_ends)} has no matching Cα residue)"
                    )
                    continue

                group_column = next(
                    (normalized_fields[name] for name in ("edge_group", "group", "category") if name in normalized_fields),
                    "",
                )
                color_column = next(
                    (normalized_fields[name] for name in ("color", "edge_color") if name in normalized_fields),
                    "",
                )
                group = str(row.get(group_column, "") or "").strip() or "Imported edges"
                requested_color = str(row.get(color_column, "") or "").strip()
                safe_color = requested_color if re.fullmatch(r"#[0-9a-fA-F]{3,8}|[A-Za-z]+", requested_color) else ""
                color = safe_color or ""
                key = (group, color)
                if key not in external_graph_groups:
                    palette = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#6A3D9A")
                    external_graph_groups[key] = {
                        "group": group,
                        "color": color or palette[len(external_graph_groups) % len(palette)],
                        "x": [], "y": [], "z": [], "text": [], "edge_count": 0,
                    }

                metadata_fields = []
                endpoint_and_style_columns = {
                    normalized_fields[alias]
                    for names in aliases.values()
                    for alias in names
                    if alias in normalized_fields
                }
                endpoint_and_style_columns.update(
                    normalized_fields[name]
                    for name in ("edge_group", "group", "category", "color", "edge_color")
                    if name in normalized_fields
                )
                for column, value in row.items():
                    value = str(value or "").strip()
                    if column not in endpoint_and_style_columns and value:
                        metadata_fields.append(f"{html.escape(str(column))}: {html.escape(value)}")
                hover = (
                    "<b>External residue-level graph</b><br>"
                    "Imported residue-pair edge; not a CNI coordination contact<br>"
                    f"Source: {html.escape(source['resname'])} {html.escape(source_label)}<br>"
                    f"Target: {html.escape(target['resname'])} {html.escape(target_label)}<br>"
                    f"Edge group: {html.escape(group)}"
                )
                if requested_color:
                    hover += f"<br>Edge color: {html.escape(requested_color)}"
                if metadata_fields:
                    hover += "<br>" + "<br>".join(metadata_fields)
                trace = external_graph_groups[key]
                trace["x"].extend([float(source["coordinates"][0]), float(target["coordinates"][0]), None])
                trace["y"].extend([float(source["coordinates"][1]), float(target["coordinates"][1]), None])
                trace["z"].extend([float(source["coordinates"][2]), float(target["coordinates"][2]), None])
                trace["text"].extend([hover, hover, None])
                trace["edge_count"] += 1

    external_graph_traces = list(external_graph_groups.values())

    if contacts is None:
        contact_rows: List[Dict] = []
    elif hasattr(contacts, "to_dict"):
        contact_rows = contacts.to_dict("records")
    else:
        contact_rows = list(contacts)

    def _contact_atom(row: Dict, prefix: str) -> Dict:
        return {
            "model_id": row.get(f"{prefix}_model_id", ""),
            "residue": row.get(f"{prefix}_resname", ""),
            "residue_number": row.get(f"{prefix}_resnum", ""),
            "chain": row.get(f"{prefix}_chain", ""),
            "insertion_code": row.get(f"{prefix}_insertion_code", ""),
            "hetero_flag": row.get(f"{prefix}_hetero_flag", ""),
            "name": row.get(f"{prefix}_atom", ""),
        }

    # Preserve compatibility with callers that pass atom records without the
    # additive coordination_role column by deriving primary endpoints from
    # the contact rows.
    primary_contact_keys = set()
    for row in contact_rows:
        direct_value = row.get("direct_coordination", False)
        is_direct = direct_value is True or str(direct_value).strip().lower() in {"true", "1", "yes"}
        is_primary = is_direct or str(row.get("contact_role", "")).strip().lower() == "primary_motif_contact"
        if is_primary:
            primary_contact_keys.update(
                _atom_key(_contact_atom(row, prefix)) for prefix in ("src", "dst")
            )
    for atom in focused_atoms:
        if (
            atom.get("coordination_role") == "active_site_component"
            and _atom_key(atom) in primary_contact_keys
        ):
            atom["coordination_role"] = "primary_coordinator"

    def _shell_number(label: Any) -> Optional[int]:
        normalized = str(label or "").strip().upper()
        named = {"COFACTOR": 0, "PCS": 1, "SCS": 2, "TCS": 3}
        if normalized in named:
            return named[normalized]
        if normalized.startswith("SHELL"):
            try:
                return int(normalized[5:])
            except ValueError:
                return None
        return None

    contact_segments: Dict[str, List[tuple]] = {
        "primary": [], "secondary": [], "tertiary": [],
    }
    contact_candidates: Dict[str, List[tuple]] = {
        "primary": [], "secondary": [], "tertiary": [],
    }
    rooted_atom_keys = set()

    def _contact_atom_keys(atom: Dict) -> set:
        # Keep both model-aware and model-agnostic keys so viewer callers that
        # omit model identity still get the same rooted-network behavior.
        return {
            _atom_key(atom),
            _atom_key(atom, include_model=False),
        }

    for row in contact_rows:
        source = _contact_atom(row, "src")
        destination = _contact_atom(row, "dst")
        source_coordinates = coordinates_by_key.get(_atom_key(source))
        destination_coordinates = coordinates_by_key.get(_atom_key(destination))
        if source_coordinates is None:
            source_coordinates = coordinates_by_base_key.get(_atom_key(source, include_model=False))
        if destination_coordinates is None:
            destination_coordinates = coordinates_by_base_key.get(_atom_key(destination, include_model=False))
        if source_coordinates is None or destination_coordinates is None:
            continue
        direct_value = row.get("direct_coordination", False)
        is_direct = direct_value is True or str(direct_value).strip().lower() in {"true", "1", "yes"}
        contact_role = str(row.get("contact_role", "")).strip().lower()
        is_inferred_hbond = contact_role == "primary_motif_contact"
        link_type = str(row.get("link_type", "")).strip().lower()
        dst_shell_number = _shell_number(row.get("dst_shell"))
        if is_direct or is_inferred_hbond:
            category = "primary"
        elif link_type == "pcs->scs" or dst_shell_number == 2:
            category = "secondary"
        elif link_type == "scs->tcs" or (dst_shell_number is not None and dst_shell_number >= 3):
            category = "tertiary"
        else:
            # A non-direct cofactor->PCS pair is only a distance-qualified
            # proximity pair (for example, O2 near porphyrin atoms), not a
            # coordination edge. Do not draw it as a network contact.
            continue
        if is_inferred_hbond:
            destination_motif = str(row.get("dst_motif", "")).strip().lower()
            destination_element = str(row.get("dst_element", "")).strip().upper()
            if destination_motif in {"water", "hydroxyl", "phenol", "thiol"} and destination_element in {"O", "S"}:
                hbond_title = (
                    "Inferred O–H···O hydrogen bond"
                    if destination_motif == "water" and destination_element == "O"
                    else "Inferred O–H···X hydrogen bond"
                )
                label = (
                    f"{hbond_title}<br>"
                    f"Acceptor: {row.get('src_resname', '?')} {row.get('src_resnum', '?')}:{row.get('src_atom', '?')}"
                    f" [{row.get('src_motif', '?')}]<br>"
                    f"Donor atom: {row.get('dst_resname', '?')} {row.get('dst_resnum', '?')}:{row.get('dst_atom', '?')}"
                    f" [{row.get('dst_motif', '?')}]<br>"
                    f"Heavy-atom distance: {float(row.get('distance_A', 0.0)):.3f} Å"
                )
            else:
                label = (
                    "Primary heme-motif contact<br>"
                    f"{row.get('src_resname', '?')} {row.get('src_resnum', '?')}:{row.get('src_atom', '?')}"
                    f" [{row.get('src_motif', '?')}] → "
                    f"{row.get('dst_resname', '?')} {row.get('dst_resnum', '?')}:{row.get('dst_atom', '?')}"
                    f" [{row.get('dst_motif', '?')}]<br>"
                    f"Distance: {float(row.get('distance_A', 0.0)):.3f} Å"
                )
        else:
            label = (
                f"{row.get('link_type', 'contact')}<br>"
                f"{row.get('src_resname', '?')} {row.get('src_resnum', '?')}:{row.get('src_atom', '?')}"
                f" [{row.get('src_motif', '?')}] → "
                f"{row.get('dst_resname', '?')} {row.get('dst_resnum', '?')}:{row.get('dst_atom', '?')}"
                f" [{row.get('dst_motif', '?')}]<br>"
                f"Distance: {float(row.get('distance_A', 0.0)):.3f} Å"
            )
        source_keys = _contact_atom_keys(source)
        destination_keys = _contact_atom_keys(destination)
        contact_candidates[category].append(
            (
                source_keys,
                destination_keys,
                source_coordinates,
                destination_coordinates,
                label,
            )
        )
        if category == "primary":
            rooted_atom_keys.update(source_keys)
            rooted_atom_keys.update(destination_keys)

    # Only expose deeper-shell overlays that descend from an actual Primary
    # contact. Shell membership remains available in the atom layers, but an
    # unrooted PCS->SCS or SCS->TCS geometric branch is not drawn as part of
    # the coordination network.
    for source_keys, destination_keys, source_coordinates, destination_coordinates, label in contact_candidates["primary"]:
        contact_segments["primary"].append((source_coordinates, destination_coordinates, label))
        rooted_atom_keys.update(source_keys)
        rooted_atom_keys.update(destination_keys)

    reachable_atom_keys = set(rooted_atom_keys)
    network_atom_keys = set(rooted_atom_keys)
    for category in ("secondary", "tertiary"):
        for source_keys, destination_keys, source_coordinates, destination_coordinates, label in contact_candidates[category]:
            if not source_keys.intersection(reachable_atom_keys):
                continue
            contact_segments[category].append((source_coordinates, destination_coordinates, label))
            reachable_atom_keys.update(destination_keys)
            network_atom_keys.update(source_keys)
            network_atom_keys.update(destination_keys)

    # A focused shell contains both actual network atoms and nearby active-site
    # context. Keep those as separate marker layers so "Coordination Network
    # only" does not silently show unconnected geometric context atoms.
    network_atoms = [
        atom for atom in focused_atoms
        if str(atom.get("shell", "")) == "Cofactor"
        or _contact_atom_keys(atom).intersection(network_atom_keys)
    ]
    network_atom_ids = {id(atom) for atom in network_atoms}
    context_atoms = [atom for atom in focused_atoms if id(atom) not in network_atom_ids]
    cofactor_residue_keys = {
        _residue_key(atom) for atom in focused_atoms
        if str(atom.get("shell", "")) == "Cofactor"
    }
    network_residue_keys = {_residue_key(atom) for atom in network_atoms}
    network_motif_residue_keys = {
        _residue_key(atom)
        for atom in network_atoms
        if str(atom.get("motif", "")) not in {"unknown_motif", "backbone"}
    }
    network_motif_atoms = [
        atom for atom in motif_atoms
        if _residue_key(atom) in network_motif_residue_keys
        and _residue_key(atom) not in cofactor_residue_keys
    ]
    network_backbone_atoms = [
        atom for atom in full_atoms
        if _residue_key(atom) in network_residue_keys
        and str(atom.get("motif", "")) == "backbone"
    ]

    def _line_arrays(segments: Sequence[tuple]) -> tuple:
        xs: List[Any] = []
        ys: List[Any] = []
        zs: List[Any] = []
        hover: List[Any] = []
        for first, second, label in segments:
            xs.extend([float(first[0]), float(second[0]), None])
            ys.extend([float(first[1]), float(second[1]), None])
            zs.extend([float(first[2]), float(second[2]), None])
            hover.extend([label, label, None])
        return xs, ys, zs, hover

    def _bond_arrays(bonds: Sequence[tuple]) -> tuple:
        xs: List[Any] = []
        ys: List[Any] = []
        zs: List[Any] = []
        for first, second in bonds:
            xs.extend([float(first[0]), float(second[0]), None])
            ys.extend([float(first[1]), float(second[1]), None])
            zs.extend([float(first[2]), float(second[2]), None])
        return xs, ys, zs

    coordinator_atoms = [
        atom for atom in focused_atoms
        if str(atom.get("coordination_role", "")) == "primary_coordinator"
    ]
    cofactor_full_atoms = [
        atom for atom in full_atoms
        if _residue_key(atom) in cofactor_residue_keys
    ]
    non_cofactor_focused_atoms = [
        atom for atom in focused_atoms
        if _residue_key(atom) not in cofactor_residue_keys
    ]
    non_cofactor_full_atoms = [
        atom for atom in full_atoms
        if _residue_key(atom) not in cofactor_residue_keys
    ]
    cofactor_bonds = generate_residue_bonds(cofactor_full_atoms, bond_lookup_table)
    focused_bonds = generate_residue_bonds(non_cofactor_focused_atoms, bond_lookup_table)
    network_residue_atoms = [
        atom for atom in full_atoms
        if _residue_key(atom) in network_residue_keys
        and _residue_key(atom) not in cofactor_residue_keys
    ]
    network_residue_bonds = generate_residue_bonds(network_residue_atoms, bond_lookup_table)
    motif_bonds = generate_residue_bonds(network_motif_atoms, bond_lookup_table)
    protein_atom_names_by_residue: Dict[tuple, set] = {}
    for atom in full_atoms:
        protein_atom_names_by_residue.setdefault(_residue_key(atom), set()).add(
            str(atom.get("name", "")).strip().upper()
        )
    protein_residue_keys = {
        residue_key for residue_key, atom_names in protein_atom_names_by_residue.items()
        if {"N", "CA", "C"}.issubset(atom_names)
    }
    protein_atoms = [
        atom for atom in full_atoms
        if _residue_key(atom) in protein_residue_keys
    ]
    # The Backbone bond mode is the whole-protein structural context: include
    # each protein residue's backbone and its R-group bonds, while excluding
    # cofactors, waters, and other non-protein residues.
    protein_backbone_bonds = generate_residue_bonds(protein_atoms, bond_lookup_table)
    full_bonds = generate_residue_bonds(non_cofactor_full_atoms, bond_lookup_table)
    network_coordinates = np.asarray([atom["coordinates"] for atom in network_atoms], dtype=float)
    context_coordinates = np.asarray([atom["coordinates"] for atom in context_atoms], dtype=float)
    coordinator_coordinates = np.asarray([atom["coordinates"] for atom in coordinator_atoms], dtype=float)
    network_motif_coordinates = np.asarray([atom["coordinates"] for atom in network_motif_atoms], dtype=float)
    network_backbone_coordinates = np.asarray([atom["coordinates"] for atom in network_backbone_atoms], dtype=float)
    motif_coordinates = np.asarray([atom["coordinates"] for atom in motif_atoms], dtype=float)
    backbone_coordinates = np.asarray([atom["coordinates"] for atom in backbone_atoms], dtype=float)
    full_coordinates = np.asarray([atom["coordinates"] for atom in full_atoms], dtype=float)
    shell_colors = {"Cofactor": "#111111", "PCS": "#2166ac", "SCS": "#c51b7d", "TCS": "#1b7837"}
    fallback_shell_colors = ["#762a83", "#e08214", "#008837", "#7f3b08"]
    unique_motifs = sorted({str(atom["motif"]) for atom in focused_atoms})
    motif_palette = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d"]
    motif_colors = {motif: motif_palette[index % len(motif_palette)] for index, motif in enumerate(unique_motifs)}
    network_shell_colors = [
        shell_colors.get(str(atom["shell"]), fallback_shell_colors[index % len(fallback_shell_colors)])
        for index, atom in enumerate(network_atoms)
    ]
    context_shell_colors = [
        shell_colors.get(str(atom["shell"]), fallback_shell_colors[index % len(fallback_shell_colors)])
        for index, atom in enumerate(context_atoms)
    ]
    coordinator_shell_colors = [
        shell_colors.get(str(atom["shell"]), fallback_shell_colors[index % len(fallback_shell_colors)])
        for index, atom in enumerate(coordinator_atoms)
    ]
    network_element_colors = [atom_type_colors.get(str(atom["element"]).upper(), "#777777") for atom in network_atoms]
    network_motif_colors = [motif_colors[str(atom["motif"])] for atom in network_atoms]
    context_element_colors = [atom_type_colors.get(str(atom["element"]).upper(), "#777777") for atom in context_atoms]
    context_motif_colors = [motif_colors[str(atom["motif"])] for atom in context_atoms]
    network_motif_layer_colors = [motif_colors.get(str(atom["motif"]), "#777777") for atom in network_motif_atoms]
    network_backbone_layer_colors = ["#8c8c8c"] * len(network_backbone_atoms)
    network_motif_shell_colors = ["#777777"] * len(network_motif_atoms)
    network_backbone_shell_colors = ["#8c8c8c"] * len(network_backbone_atoms)
    network_motif_element_colors = [atom_type_colors.get(str(atom["element"]).upper(), "#777777") for atom in network_motif_atoms]
    network_backbone_element_colors = [atom_type_colors.get(str(atom["element"]).upper(), "#777777") for atom in network_backbone_atoms]
    coordinator_element_colors = [atom_type_colors.get(str(atom["element"]).upper(), "#777777") for atom in coordinator_atoms]
    coordinator_motif_colors = [motif_colors[str(atom["motif"])] for atom in coordinator_atoms]
    motif_layer_colors = [motif_colors.get(str(atom["motif"]), "#777777") for atom in motif_atoms]
    backbone_layer_colors = ["#8c8c8c"] * len(backbone_atoms)
    focused_motif_shell_colors = ["#777777"] * len(motif_atoms)
    focused_backbone_shell_colors = ["#8c8c8c"] * len(backbone_atoms)
    motif_element_colors = [atom_type_colors.get(str(atom["element"]).upper(), "#777777") for atom in motif_atoms]
    backbone_element_colors = [atom_type_colors.get(str(atom["element"]).upper(), "#777777") for atom in backbone_atoms]
    network_backbone_motif_colors = [motif_colors.get(str(atom["motif"]), "#777777") for atom in network_backbone_atoms]
    focused_backbone_motif_colors = [motif_colors.get(str(atom["motif"]), "#777777") for atom in backbone_atoms]

    def _conservation_colors(atoms: Sequence[Dict]) -> List[str]:
        return [sample_colorscale("YlGnBu", [_conservation_value(atom)])[0] for atom in atoms]

    network_conservation_colors = _conservation_colors(network_atoms)
    coordinator_conservation_colors = _conservation_colors(coordinator_atoms)
    context_conservation_colors = _conservation_colors(context_atoms)
    network_motif_conservation_colors = _conservation_colors(network_motif_atoms)
    network_backbone_conservation_colors = _conservation_colors(network_backbone_atoms)
    focused_motif_conservation_colors = _conservation_colors(motif_atoms)
    focused_backbone_conservation_colors = _conservation_colors(backbone_atoms)

    def _shell_hover(atoms: Sequence[Dict]) -> List[str]:
        return [
        f"Shell: {atom['shell']}<br>Residue: {atom['residue']} {atom['residue_number']} {atom['chain']}<br>"
        f"Atom: {atom['name']} ({atom['element']})<br>Motif: {atom['motif']}<br>"
        f"Role: {atom['coordination_role']}<br>Model: {atom['model_id']}"
        + (f"<br>{conservation_label}: {_conservation_value(atom):.0%}" if conservation_by_residue is not None else "")
            for atom in atoms
        ]
    network_hover = _shell_hover(network_atoms)
    context_hover = _shell_hover(context_atoms)
    coordinator_hover = [
        f"Actual primary coordinator<br>Shell: {atom['shell']}<br>"
        f"Residue: {atom['residue']} {atom['residue_number']} {atom['chain']}<br>"
        f"Atom: {atom['name']} ({atom['element']})<br>Motif: {atom['motif']}<br>Model: {atom['model_id']}"
        for atom in coordinator_atoms
    ]
    full_hover = [
        f"Residue: {atom['residue']} {atom['residue_number']} {atom['chain']}<br>"
        f"Atom: {atom['name']} ({atom['element']})<br>Motif: {atom['motif']}<br>Model: {atom['model_id']}"
        for atom in full_atoms
    ]
    motif_hover = [
        f"Motif layer<br>Residue: {atom['residue']} {atom['residue_number']} {atom['chain']}<br>"
        f"Atom: {atom['name']} ({atom['element']})<br>Motif: {atom['motif']}<br>Model: {atom['model_id']}"
        for atom in motif_atoms
    ]
    network_motif_hover = [
        f"Network motif layer<br>Residue: {atom['residue']} {atom['residue_number']} {atom['chain']}<br>"
        f"Atom: {atom['name']} ({atom['element']})<br>Motif: {atom['motif']}<br>Model: {atom['model_id']}"
        for atom in network_motif_atoms
    ]
    backbone_hover = [
        f"Backbone layer<br>Residue: {atom['residue']} {atom['residue_number']} {atom['chain']}<br>"
        f"Atom: {atom['name']} ({atom['element']})<br>Model: {atom['model_id']}"
        for atom in backbone_atoms
    ]
    network_backbone_hover = [
        f"Network backbone layer<br>Residue: {atom['residue']} {atom['residue_number']} {atom['chain']}<br>"
        f"Atom: {atom['name']} ({atom['element']})<br>Model: {atom['model_id']}"
        for atom in network_backbone_atoms
    ]
    network_x, network_y, network_z = zip(*network_coordinates) if network_coordinates.size else ([], [], [])
    context_x, context_y, context_z = zip(*context_coordinates) if context_coordinates.size else ([], [], [])
    coordinator_x, coordinator_y, coordinator_z = zip(*coordinator_coordinates) if coordinator_coordinates.size else ([], [], [])
    network_motif_x, network_motif_y, network_motif_z = zip(*network_motif_coordinates) if network_motif_coordinates.size else ([], [], [])
    network_backbone_x, network_backbone_y, network_backbone_z = zip(*network_backbone_coordinates) if network_backbone_coordinates.size else ([], [], [])
    motif_x, motif_y, motif_z = zip(*motif_coordinates) if motif_coordinates.size else ([], [], [])
    backbone_x, backbone_y, backbone_z = zip(*backbone_coordinates) if backbone_coordinates.size else ([], [], [])
    full_x, full_y, full_z = zip(*full_coordinates) if full_coordinates.size else ([], [], [])
    focused_bond_x, focused_bond_y, focused_bond_z = _bond_arrays(focused_bonds)
    cofactor_bond_x, cofactor_bond_y, cofactor_bond_z = _bond_arrays(cofactor_bonds)
    network_residue_bond_x, network_residue_bond_y, network_residue_bond_z = _bond_arrays(network_residue_bonds)
    motif_bond_x, motif_bond_y, motif_bond_z = _bond_arrays(motif_bonds)
    backbone_bond_x, backbone_bond_y, backbone_bond_z = _bond_arrays(protein_backbone_bonds)
    full_bond_x, full_bond_y, full_bond_z = _bond_arrays(full_bonds)
    primary_x, primary_y, primary_z, primary_hover = _line_arrays(contact_segments["primary"])
    secondary_x, secondary_y, secondary_z, secondary_hover = _line_arrays(contact_segments["secondary"])
    tertiary_x, tertiary_y, tertiary_z, tertiary_hover = _line_arrays(contact_segments["tertiary"])

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=full_x, y=full_y, z=full_z, mode="markers",
        marker=dict(size=2, color="#bdbdbd", opacity=0.22),
        text=full_hover, hoverinfo="text", name="Full protein", visible=False, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=network_x, y=network_y, z=network_z, mode="markers",
        marker=dict(size=5, color=network_shell_colors, opacity=0.35),
        text=network_hover, hoverinfo="text", name="Coordination Network atoms", visible=True, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=coordinator_x, y=coordinator_y, z=coordinator_z, mode="markers",
        marker=dict(size=9, color=coordinator_shell_colors, opacity=1.0, line=dict(color="#111111", width=1)),
        text=coordinator_hover, hoverinfo="text", name="Actual coordinators", visible=True, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=context_x, y=context_y, z=context_z, mode="markers",
        marker=dict(size=5, color=context_shell_colors, opacity=0.35),
        text=context_hover, hoverinfo="text", name="Active-site context", visible=False, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=network_motif_x, y=network_motif_y, z=network_motif_z, mode="markers",
        marker=dict(size=5, color=network_motif_layer_colors, opacity=0.95),
        text=network_motif_hover, hoverinfo="text", name="Network motif atoms", visible=False, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=network_backbone_x, y=network_backbone_y, z=network_backbone_z, mode="markers",
        marker=dict(size=4, color=network_backbone_layer_colors, opacity=0.9),
        text=network_backbone_hover, hoverinfo="text", name="Network backbone atoms", visible=False, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=motif_x, y=motif_y, z=motif_z, mode="markers",
        marker=dict(size=5, color=motif_layer_colors, opacity=0.95),
        text=motif_hover, hoverinfo="text", name="Focused motif atoms", visible=False, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=backbone_x, y=backbone_y, z=backbone_z, mode="markers",
        marker=dict(size=4, color=backbone_layer_colors, opacity=0.9),
        text=backbone_hover, hoverinfo="text", name="Focused backbone atoms", visible=False, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=cofactor_bond_x, y=cofactor_bond_y, z=cofactor_bond_z, mode="lines",
        line=dict(color="#222222", width=4), hoverinfo="skip", name="Cofactor bonds", visible=True, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=focused_bond_x, y=focused_bond_y, z=focused_bond_z, mode="lines",
        line=dict(color="#222222", width=4), hoverinfo="skip", name="Focused bonds", visible=False, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=network_residue_bond_x, y=network_residue_bond_y, z=network_residue_bond_z, mode="lines",
        line=dict(color="#222222", width=4), hoverinfo="skip", name="Network residue bonds", visible=True, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=full_bond_x, y=full_bond_y, z=full_bond_z, mode="lines",
        line=dict(color="#999999", width=1), hoverinfo="skip", name="Full-protein bonds", visible=False, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=motif_bond_x, y=motif_bond_y, z=motif_bond_z, mode="lines",
        line=dict(color="#555555", width=3), hoverinfo="skip", name="Motif bonds", visible=False, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=backbone_bond_x, y=backbone_bond_y, z=backbone_bond_z, mode="lines",
        line=dict(color="#aaaaaa", width=2), hoverinfo="skip", name="Backbone bonds (full protein + R-groups)", visible=False, showlegend=False,
    ))
    fig.add_trace(go.Scatter3d(
        x=primary_x, y=primary_y, z=primary_z, mode="lines",
        # Keep all contact levels gray and dotted, with a deliberately clear
        # visual hierarchy.  A 2 px tertiary trace becomes effectively
        # invisible in a large Plotly 3D scene, especially over the grid and
        # structural bonds, even when its data are present.
        line=dict(color="#6f7378", width=7, dash="dot"), opacity=1.0,
        text=primary_hover, hoverinfo="text",
        name="Primary coordination", visible=True,
    ))
    fig.add_trace(go.Scatter3d(
        x=secondary_x, y=secondary_y, z=secondary_z, mode="lines",
        line=dict(color="#6f7378", width=5, dash="dot"), opacity=1.0,
        text=secondary_hover, hoverinfo="text",
        name="Secondary coordination", visible=True,
    ))
    fig.add_trace(go.Scatter3d(
        x=tertiary_x, y=tertiary_y, z=tertiary_z, mode="lines",
        line=dict(color="#6f7378", width=3, dash="dot"), opacity=1.0,
        text=tertiary_hover, hoverinfo="text",
        name="Tertiary coordination", visible=True,
    ))
    for external_trace in external_graph_traces:
        fig.add_trace(go.Scatter3d(
            x=external_trace["x"], y=external_trace["y"], z=external_trace["z"],
            mode="lines",
            line=dict(color=external_trace["color"], width=6, dash="solid"),
            opacity=0.95,
            text=external_trace["text"], hoverinfo="text",
            name=f"External residue graph — {external_trace['group']} ({external_trace['edge_count']} edges)",
            visible=True, showlegend=False,
        ))

    full_index, network_index, coordinator_index, context_index = 0, 1, 2, 3
    network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index = 4, 5, 6, 7
    cofactor_bond_index, focused_bond_index, residue_bond_index, full_bond_index, motif_bond_index, backbone_bond_index = 8, 9, 10, 11, 12, 13
    primary_index, secondary_index, tertiary_index = 14, 15, 16
    all_trace_indices = list(range(17))
    external_graph_indices = list(range(17, 17 + len(external_graph_traces)))
    network_preset_visibility = [False, True, True, False, False, False, False, False, True, False, True, False, False, False, True, True, True] + [False] * len(external_graph_indices)
    network_motifs_preset_visibility = [False, True, True, False, True, False, False, False, True, False, False, False, True, False, True, True, True] + [False] * len(external_graph_indices)
    structural_context_preset_visibility = [False, True, True, True, False, False, True, True, True, False, False, False, False, True, True, True, True] + [False] * len(external_graph_indices)
    full_protein_preset_visibility = [True, True, True, False, False, False, False, False, True, False, False, True, False, False, True, True, True] + [False] * len(external_graph_indices)
    external_graph_preset_visibility = [True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False] + [True] * len(external_graph_indices)
    view_preset_buttons = [
        dict(label="Network", method="update", args=[{"visible": network_preset_visibility}, {f"annotations[{10 if external_residue_network_csv else 9}].text": "View: Network"}]),
        dict(label="Network + motifs", method="update", args=[{"visible": network_motifs_preset_visibility}, {f"annotations[{10 if external_residue_network_csv else 9}].text": "View: Network + motifs"}]),
        dict(label="Structural context", method="update", args=[{"visible": structural_context_preset_visibility}, {f"annotations[{10 if external_residue_network_csv else 9}].text": "View: Structural context"}]),
        dict(label="Full protein context", method="update", args=[{"visible": full_protein_preset_visibility}, {f"annotations[{10 if external_residue_network_csv else 9}].text": "View: Full protein context"}]),
    ]
    if external_residue_network_csv:
        view_preset_buttons.append(
            dict(label="External residue graph", method="update", args=[{"visible": external_graph_preset_visibility}, {f"annotations[10].text": "View: External residue graph (residue level)"}])
        )
    color_by_buttons = [
        dict(label="Coordination shells", method="restyle", args=[{"marker.color": [network_shell_colors, coordinator_shell_colors, context_shell_colors, network_motif_shell_colors, network_backbone_shell_colors, focused_motif_shell_colors, focused_backbone_shell_colors]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
        dict(label="Elements", method="restyle", args=[{"marker.color": [network_element_colors, coordinator_element_colors, context_element_colors, network_motif_element_colors, network_backbone_element_colors, motif_element_colors, backbone_element_colors]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
        dict(label="Motifs", method="restyle", args=[{"marker.color": [network_motif_colors, coordinator_motif_colors, context_motif_colors, network_motif_layer_colors, network_backbone_motif_colors, motif_layer_colors, focused_backbone_motif_colors]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
    ]
    if conservation_by_residue is not None:
        color_by_buttons.append(
            dict(
                label=conservation_label,
                method="restyle",
                args=[{"marker.color": [network_conservation_colors, coordinator_conservation_colors, context_conservation_colors, network_motif_conservation_colors, network_backbone_conservation_colors, focused_motif_conservation_colors, focused_backbone_conservation_colors]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]],
            )
        )
    external_layer_menu = []
    if external_residue_network_csv:
        external_layer_menu.append(
            dict(type="dropdown", direction="down", x=0.01, y=0.40, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                dict(label="On", method="restyle", args=[{"visible": [True] * len(external_graph_indices)}, external_graph_indices]),
                dict(label="Off", method="restyle", args=[{"visible": [False] * len(external_graph_indices)}, external_graph_indices]),
            ])
        )
    contact_menu_y = 0.29 if external_residue_network_csv else 0.40
    grid_menu_y = 0.20 if external_residue_network_csv else 0.29
    background_menu_y = 0.11 if external_residue_network_csv else 0.18
    labels_menu_y = 0.02 if external_residue_network_csv else 0.07
    fig.update_layout(
        title={"text": f"Cofactor coordination network: {cofactor_resname} in {pdb_name}", "x": 0.60, "font": {"size": 22}},
        scene=dict(
            xaxis=dict(title="X (Å)", showgrid=True, gridcolor="#d9e1ec", showbackground=True, backgroundcolor="#e5ecf6"),
            yaxis=dict(title="Y (Å)", showgrid=True, gridcolor="#d9e1ec", showbackground=True, backgroundcolor="#e5ecf6"),
            zaxis=dict(title="Z (Å)", showgrid=True, gridcolor="#d9e1ec", showbackground=True, backgroundcolor="#e5ecf6"),
            bgcolor="#e5ecf6",
            aspectmode="data",
            domain=dict(x=[0.20, 1.0], y=[0.0, 1.0]),
        ),
        margin=dict(l=0, r=0, t=70, b=0),
        legend=dict(x=0.78, y=0.98, xanchor="left", yanchor="top", bgcolor="rgba(255,255,255,0.85)", font=dict(size=10)),
        updatemenus=[
            dict(type="dropdown", direction="down", x=0.01, y=0.95, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                *view_preset_buttons,
            ]),
            dict(type="dropdown", direction="down", x=0.01, y=0.84, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                dict(label="On", method="restyle", args=[{"visible": [True]}, [full_index]]),
                dict(label="Off", method="restyle", args=[{"visible": [False]}, [full_index]]),
            ]),
            dict(type="dropdown", direction="down", x=0.01, y=0.73, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                dict(label="Network atoms only", method="restyle", args=[{"visible": [True, True, False, False, False, False, False]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
                dict(label="Network atoms + motifs", method="restyle", args=[{"visible": [True, True, False, True, False, False, False]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
                dict(label="Network atoms + backbone", method="restyle", args=[{"visible": [True, True, False, False, True, False, False]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
                dict(label="Active-site context only", method="restyle", args=[{"visible": [False, False, True, False, False, False, False]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
                dict(label="Primary coordinators only", method="restyle", args=[{"visible": [False, True, False, False, False, False, False]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
                dict(label="Network motifs only", method="restyle", args=[{"visible": [False, False, False, True, False, False, False]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
                dict(label="Network backbone only", method="restyle", args=[{"visible": [False, False, False, False, True, False, False]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
                dict(label="Focused motifs only", method="restyle", args=[{"visible": [False, False, False, False, False, True, False]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
                dict(label="Focused backbone only", method="restyle", args=[{"visible": [False, False, False, False, False, False, True]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
                dict(label="All focused atom layers", method="restyle", args=[{"visible": [True, True, True, False, False, True, True]}, [network_index, coordinator_index, context_index, network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index]]),
            ]),
            dict(type="dropdown", direction="down", x=0.01, y=0.62, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                *color_by_buttons,
            ]),
            dict(type="dropdown", direction="down", x=0.01, y=0.51, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                dict(label="Focused", method="restyle", args=[{"visible": [True, False, False, False, False]}, [focused_bond_index, residue_bond_index, full_bond_index, motif_bond_index, backbone_bond_index]]),
                dict(label="Residues", method="restyle", args=[{"visible": [False, True, False, False, False]}, [focused_bond_index, residue_bond_index, full_bond_index, motif_bond_index, backbone_bond_index]]),
                dict(label="Motif", method="restyle", args=[{"visible": [False, False, False, True, False]}, [focused_bond_index, residue_bond_index, full_bond_index, motif_bond_index, backbone_bond_index]]),
                dict(label="Backbone", method="restyle", args=[{"visible": [False, False, False, False, True]}, [focused_bond_index, residue_bond_index, full_bond_index, motif_bond_index, backbone_bond_index]]),
                dict(label="Full protein", method="restyle", args=[{"visible": [False, False, True, False, False]}, [focused_bond_index, residue_bond_index, full_bond_index, motif_bond_index, backbone_bond_index]]),
                dict(label="All", method="restyle", args=[{"visible": [True, True, True, True, True]}, [focused_bond_index, residue_bond_index, full_bond_index, motif_bond_index, backbone_bond_index]]),
                dict(label="None", method="restyle", args=[{"visible": [False, False, False, False, False]}, [focused_bond_index, residue_bond_index, full_bond_index, motif_bond_index, backbone_bond_index]]),
            ]),
            *external_layer_menu,
            dict(type="dropdown", direction="down", x=0.01, y=contact_menu_y, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                dict(label="All", method="restyle", args=[{"visible": [True, True, True]}, [primary_index, secondary_index, tertiary_index]]),
                dict(label="Primary", method="restyle", args=[{"visible": [True, False, False]}, [primary_index, secondary_index, tertiary_index]]),
                dict(label="Secondary", method="restyle", args=[{"visible": [False, True, False]}, [primary_index, secondary_index, tertiary_index]]),
                dict(label="Tertiary", method="restyle", args=[{"visible": [False, False, True]}, [primary_index, secondary_index, tertiary_index]]),
                dict(label="None", method="restyle", args=[{"visible": [False, False, False]}, [primary_index, secondary_index, tertiary_index]]),
            ]),
            dict(type="dropdown", direction="down", x=0.01, y=grid_menu_y, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                dict(label="On", method="relayout", args=[{
                    "scene.xaxis.showgrid": True,
                    "scene.yaxis.showgrid": True,
                    "scene.zaxis.showgrid": True,
                }]),
                dict(label="Off", method="relayout", args=[{
                    "scene.xaxis.showgrid": False,
                    "scene.yaxis.showgrid": False,
                    "scene.zaxis.showgrid": False,
                }]),
            ]),
            dict(type="dropdown", direction="up", x=0.01, y=background_menu_y, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                dict(label="On", method="relayout", args=[{
                    "scene.bgcolor": "#e5ecf6",
                    "scene.xaxis.showbackground": True,
                    "scene.yaxis.showbackground": True,
                    "scene.zaxis.showbackground": True,
                    "scene.xaxis.backgroundcolor": "#e5ecf6",
                    "scene.yaxis.backgroundcolor": "#e5ecf6",
                    "scene.zaxis.backgroundcolor": "#e5ecf6",
                }]),
                dict(label="Off", method="relayout", args=[{
                    "scene.bgcolor": "rgba(0,0,0,0)",
                    "scene.xaxis.showbackground": False,
                    "scene.yaxis.showbackground": False,
                    "scene.zaxis.showbackground": False,
                }]),
            ]),
            dict(type="dropdown", direction="up", x=0.01, y=labels_menu_y, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                dict(label="On", method="relayout", args=[{
                    "scene.xaxis.showticklabels": True,
                    "scene.yaxis.showticklabels": True,
                    "scene.zaxis.showticklabels": True,
                    "scene.xaxis.title.text": "X (Å)",
                    "scene.yaxis.title.text": "Y (Å)",
                    "scene.zaxis.title.text": "Z (Å)",
                }]),
                dict(label="Off", method="relayout", args=[{
                    "scene.xaxis.showticklabels": False,
                    "scene.yaxis.showticklabels": False,
                    "scene.zaxis.showticklabels": False,
                    "scene.xaxis.title.text": "",
                    "scene.yaxis.title.text": "",
                    "scene.zaxis.title.text": "",
                }]),
            ]),
        ],
        annotations=[
            dict(text="View preset", x=0.01, y=0.995, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Display full protein structure", x=0.01, y=0.885, xref="paper", yref="paper", showarrow=False, font=dict(size=10, color="#52627a"), xanchor="left"),
            dict(text="Atom layers", x=0.01, y=0.775, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Color by", x=0.01, y=0.665, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Bonds", x=0.01, y=0.555, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="External graph (residue level)" if external_residue_network_csv else "Contacts", x=0.01, y=0.445, xref="paper", yref="paper", showarrow=False, font=dict(size=10 if external_residue_network_csv else 11, color="#52627a"), xanchor="left"),
            dict(text="Contacts" if external_residue_network_csv else "Grid", x=0.01, y=0.335 if external_residue_network_csv else 0.335, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Grid" if external_residue_network_csv else "Background", x=0.01, y=0.245 if external_residue_network_csv else 0.225, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Background" if external_residue_network_csv else "Labels", x=0.01, y=0.155 if external_residue_network_csv else 0.115, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            *([dict(text="Labels", x=0.01, y=0.065, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left")] if external_residue_network_csv else []),
            dict(text="View: Network", x=0.21, y=0.995, xref="paper", yref="paper", showarrow=False, bgcolor="rgba(255,255,255,0.8)", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=10, color="#52627a"), xanchor="left"),
        ],
    )
    viewer_config = {
        "displayModeBar": True,
        "modeBarButtonsToAdd": [
            "orbitRotation",
            "tableRotation",
            "pan3d",
            "zoom3d",
            "resetCameraDefault3d",
            "resetCameraLastSave3d",
            "hoverClosest3d",
            "toImage",
        ],
        "displaylogo": False,
        "responsive": True,
    }
    post_script = """
const coordinationPlot = document.getElementById('{plot_id}');
if (coordinationPlot) {
  const atomDetails = document.createElement('div');
  atomDetails.id = 'coordination-atom-details';
  atomDetails.textContent = 'Click an atom to inspect its details.';
  atomDetails.style.cssText = [
    'position: fixed', 'right: 14px', 'bottom: 14px', 'z-index: 20',
    'max-width: 290px', 'padding: 10px 12px', 'white-space: pre-line',
    'font: 12px Arial, sans-serif', 'line-height: 1.35',
    'color: #26364d', 'background: rgba(255,255,255,0.92)',
    'border: 1px solid #b8c4d6', 'border-radius: 4px',
    'box-shadow: 0 1px 4px rgba(0,0,0,0.12)', 'pointer-events: none'
  ].join(';');
  document.body.appendChild(atomDetails);
  coordinationPlot.on('plotly_click', function(eventData) {
    const point = eventData && eventData.points && eventData.points[0];
    if (!point || !point.data || point.data.mode !== 'markers') return;
    const details = String(point.text || 'No atom details').replace(/<br\\s*\\/?>/gi, '\\n');
    atomDetails.textContent = (point.data.name || 'Selected atom') + '\\n' + details;
  });
}
"""
    if external_residue_network_csv:
        unmapped_lines = "\n".join(external_graph_unmapped) if external_graph_unmapped else "Unmapped edges: none"
        report_text = (
            "External residue-level graph\n"
            f"Mapped edges: {external_graph_input_count - len(external_graph_unmapped)}/{external_graph_input_count}\n"
            f"Unmapped edges: {len(external_graph_unmapped)}\n{unmapped_lines}"
        )
        post_script += f"""
const externalGraphStatus = {json.dumps(report_text, ensure_ascii=False)};
const externalGraphPanel = document.createElement('div');
externalGraphPanel.id = 'external-residue-graph-status';
externalGraphPanel.textContent = externalGraphStatus;
externalGraphPanel.style.cssText = [
  'position: fixed', 'right: 14px', 'top: 14px', 'z-index: 20',
  'max-width: 360px', 'max-height: 30vh', 'overflow: auto',
  'padding: 9px 11px', 'white-space: pre-line',
  'font: 12px Arial, sans-serif', 'line-height: 1.35',
  'color: #17324d', 'background: rgba(255,255,255,0.94)',
  'border: 1px solid #0072B2', 'border-radius: 4px',
  'box-shadow: 0 1px 4px rgba(0,0,0,0.12)', 'pointer-events: none'
].join(';');
document.body.appendChild(externalGraphPanel);
"""
    fig.write_html(
        output_filename,
        config=viewer_config,
        include_plotlyjs="cdn" if compact_html else True,
        post_script=post_script,
        div_id="coordination-network-plot",
    )
    logger.info("Cohesive interactive plot saved as '%s'", output_filename)
