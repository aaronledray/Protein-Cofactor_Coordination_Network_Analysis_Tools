# modules/plotting.py
"""
Lightweight plotting utilities.

Currently exposes:
- plot_evaluation_results(df_results, mode="CA_only", output_csv=None, highlight_labels=None)
"""

from typing import Any, List, Optional, Sequence
import logging
import pandas as pd
import matplotlib.pyplot as plt

from typing import List, Dict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (needed for 3D projection)
from .display import show_matplotlib, show_plotly
from .structure_processing import get_residue_bonds



# modules/plotting.py

from typing import Dict, List, Optional
import numpy as np
import plotly.graph_objects as go
from plotly.colors import sample_colorscale

# pull bond builders from structure_processing
from .structure_processing import generate_all_bonds, generate_residue_bonds
from .moieties import chemical_moieties
from .motif_registry import motif_for_atom

logger = logging.getLogger(__name__)




# --- add near the other imports at top of modules/plotting.py ---
import numpy as np
import plotly.graph_objects as go
from .structure_processing import generate_all_bonds, generate_residue_bonds






# --- add near other imports in modules/plotting.py ---
import numpy as np
import plotly.graph_objects as go
from .structure_processing import generate_all_bonds, generate_residue_bonds










__all__ = ["plot_evaluation_results"]


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







def static_plots_2d(
    cofactor_coords, pcs_coords, scs_coords,
    cofactor_atoms, pcs_atoms, scs_atoms,
    structure, bond_lookup,
    pdb_name: str = "structure.pdb", cofactor_resname="cofactor",
    include_bonds: bool = True,
    focused_bonds: bool = False,
    output_prefix: str = "1_"
):
    """
    Saves static 3D plots (as images) for PCS/SCS mode and Element-Based Coloring using Matplotlib.
    """
    # Set up atom_type_colors; use global if defined or set default
    global atom_type_colors
    if 'atom_type_colors' not in globals():
        atom_type_colors = {
            "C": "black",
            "N": "blue",
            "O": "red",
            "S": "yellow",
            "FE": "orange",
            "MN": "purple",
            "CA": "green",
            "CU": "goldenrod",
            "MO": "teal",
        }
    default_color = "grey"

    def get_color(atom):
        return atom_type_colors.get(atom["element"].upper(), default_color)

    # --- Generate Bonds ---
    all_bonds = []
    for residue in structure.get_residues():
        residue_name = residue.get_resname()
        if residue_name is None:
            continue
        atoms = [{
            'name': atom.get_name(),
            'element': atom.element,
            'coordinates': atom.coord,
            'residue': residue_name
        } for atom in residue]
        if not atoms:
            continue
        all_bonds.extend(get_residue_bonds(atoms, bond_lookup=bond_lookup))

    # Focused bonds: only for residues present in the cofactor/PCS/SCS sets
    allowed_residues = {(a["residue_number"], a["chain"]) for a in (cofactor_atoms + pcs_atoms + scs_atoms)}
    focused_bonds_list = []
    for residue in structure.get_residues():
        resnum = residue.get_id()[1]
        chain_id = residue.get_full_id()[2]
        if (resnum, chain_id) not in allowed_residues:
            continue
        residue_name = residue.get_resname()
        atoms = [{
            'name': atom.get_name(),
            'element': atom.element,
            'coordinates': atom.coord,
            'residue': residue_name
        } for atom in residue]
        if atoms:
            focused_bonds_list.extend(get_residue_bonds(atoms, bond_lookup=bond_lookup))

    # --- Extract Coordinates ---
    cofactor_coords = np.atleast_2d(np.array(cofactor_coords))
    pcs_coords     = np.atleast_2d(np.array(pcs_coords))
    scs_coords     = np.atleast_2d(np.array(scs_coords))
    all_atoms      = cofactor_atoms + pcs_atoms + scs_atoms
    all_coords     = np.atleast_2d(np.array([a['coordinates'] for a in all_atoms])) if all_atoms else np.empty((0,3))

    # --- Fixed (PCS/SCS) Coloring ---
    pcs_scs_colors = []
    for coord in all_coords:
        if cofactor_coords.size and np.any(np.all(coord == cofactor_coords, axis=1)):
            pcs_scs_colors.append('black')
        elif pcs_coords.size and np.any(np.all(coord == pcs_coords, axis=1)):
            pcs_scs_colors.append('blue')
        elif scs_coords.size and np.any(np.all(coord == scs_coords, axis=1)):
            pcs_scs_colors.append('fuchsia')
        else:
            pcs_scs_colors.append('grey')

    # --- Element-Based Coloring ---
    element_colors = [get_color(a) for a in all_atoms]

    # === Plot 1: PCS/SCS Mode ===
    fig1 = plt.figure(figsize=(12, 10))
    ax1 = fig1.add_subplot(111, projection='3d')

    if cofactor_coords.size:
        ax1.scatter(cofactor_coords[:,0], cofactor_coords[:,1], cofactor_coords[:,2],
                    color='black', label='Cofactor Atoms', alpha=0.9, s=10)
    if pcs_coords.size:
        ax1.scatter(pcs_coords[:,0], pcs_coords[:,1], pcs_coords[:,2],
                    color='blue', label='PCS Atoms', alpha=0.7, s=10)
    if scs_coords.size:
        ax1.scatter(scs_coords[:,0], scs_coords[:,1], scs_coords[:,2],
                    color='fuchsia', label='SCS Atoms', alpha=0.7, s=10)

    if include_bonds:
        for bond in all_bonds:
            try:
                x,y,z = zip(*bond)
                ax1.plot(x,y,z, color='gray', linewidth=1)
            except Exception:
                pass
        if focused_bonds:
            for bond in focused_bonds_list:
                try:
                    x,y,z = zip(*bond)
                    ax1.plot(x,y,z, color='black', linewidth=2)
                except Exception:
                    pass

    ax1.set_title(f"PCS/SCS Mode with Bonds for {cofactor_resname} in {pdb_name}", fontsize=16)
    ax1.set_xlabel("X Coordinate"); ax1.set_ylabel("Y Coordinate"); ax1.set_zlabel("Z Coordinate")
    ax1.legend()
    plt.tight_layout()
    out1 = f"{output_prefix}static_pcs_scs_mode_both.png" if focused_bonds else f"{output_prefix}static_pcs_scs_mode.png"
    plt.savefig(out1); logger.info("Saved PCS/SCS plot → %s", out1)
    show_matplotlib(plt)

    # === Plot 2: Element-Based Coloring ===
    fig2 = plt.figure(figsize=(12, 10))
    ax2 = fig2.add_subplot(111, projection='3d')

    # group-wise scatter using element colors
    for coords, atoms, label in zip(
        [cofactor_coords, pcs_coords, scs_coords],
        [cofactor_atoms,  pcs_atoms,  scs_atoms],
        ['Cofactor Atoms','PCS Atoms','SCS Atoms']
    ):
        coords = np.atleast_2d(coords)
        if coords.size:
            colors = [get_color(a) for a in atoms]
            ax2.scatter(coords[:,0], coords[:,1], coords[:,2],
                        c=colors, alpha=0.8, s=10, label=label)

    if include_bonds:
        for bond in all_bonds:
            try:
                x,y,z = zip(*bond)
                ax2.plot(x,y,z, color='gray', linewidth=1)
            except Exception:
                pass
        if focused_bonds:
            for bond in focused_bonds_list:
                try:
                    x,y,z = zip(*bond)
                    ax2.plot(x,y,z, color='black', linewidth=2)
                except Exception:
                    pass

    ax2.set_title(f"Element-Based Coloring with Bonds for {cofactor_resname} in {pdb_name}", fontsize=16)
    ax2.set_xlabel("X Coordinate"); ax2.set_ylabel("Y Coordinate"); ax2.set_zlabel("Z Coordinate")
    plt.tight_layout()
    out2 = f"{output_prefix}static_element_coloring_mode_both.png" if focused_bonds else f"{output_prefix}static_element_coloring_mode.png"
    plt.savefig(out2); logger.info("Saved Element plot → %s", out2)
    show_matplotlib(plt)












# UPDATED FOR COORD DOTTED LINES:





def plot_interactive_modes_with_network(
    structure,
    cofactor_atoms: List[Dict],
    pcs_atoms: List[Dict],
    scs_atoms: List[Dict],
    bond_lookup_table: Dict[str, List[List[str]]],
    pdb_name: str = "structure.pdb",
    cofactor_resname: str = "cofactor",
    atom_type_colors: Optional[Dict[str, str]] = None,
    output_filename: str = "1_template_coordination_network.html",
    # Optional: links overlay (as before)
    links_csv_path: Optional[str] = "Coord_Links.csv",
    links_rows: Optional[List[Dict[str, str]]] = None,
    show: bool = True,
):
    """
    Interactive 3D Plotly viz with:
      • Coloring toggle: Coordination Sphere / Element
      • Backbone sticks toggle: Off (focused) / On (full residue sticks)
      • Dotted link lines for cofactor→PCS and PCS→SCS (if Coord_Links.csv present or rows provided)
      • NEW: Hover shows Moiety for each atom
    """
    import os, csv
    import numpy as np
    import plotly.graph_objects as go

    # ------------------------ Moiety lookup ------------------------
    _BACKBONE_NAMES = {"N","H","CA","HA","C","O","OXT"}
    try:
        # Use your canonical table if available
        from modules.moieties import chemical_moieties as _CHEM_MOIETIES  # type: ignore
    except Exception:
        _CHEM_MOIETIES = {}

    def _moiety_of(atom: Dict) -> str:
        res = str(atom.get("residue",""))
        name = str(atom.get("name",""))
        m = _CHEM_MOIETIES.get((res, name))
        if m:
            return m
        if name in _BACKBONE_NAMES:
            return "backbone"
        return "unknown_moiety"

    # ------------------------ Colors ------------------------
    if atom_type_colors is None:
        atom_type_colors = {
            "C": "black", "N": "blue", "O": "red", "S": "yellow",
            "FE": "orange", "MN": "purple", "CA": "green", "CU": "goldenrod",
            "MO": "teal",
        }
    default_color = "grey"

    def elem_color(atom: Dict) -> str:
        return atom_type_colors.get(str(atom.get("element", "")).upper(), default_color)

    # ------------------------ Bond builders ------------------------
    focused_atoms = (cofactor_atoms or []) + (pcs_atoms or []) + (scs_atoms or [])
    minimal_focused_bonds = generate_residue_bonds(focused_atoms, bond_lookup_table)

    def _compute_full_residue_bonds_for_focused(structure, focused_atoms, bond_lookup):
        allowed_residues = {(a.get("residue_number", None), a.get("chain", None))
                            for a in focused_atoms if a is not None}
        full_bonds = []

        def residue_atoms_as_dicts(residue):
            rname = residue.get_resname()
            if not rname:
                return []
            out = []
            for at in residue:
                out.append({
                    "name": at.get_name(),
                    "element": getattr(at, "element", ""),
                    "coordinates": np.array(at.coord, dtype=float),
                    "residue": rname,
                })
            return out

        for residue in structure.get_residues():
            try:
                resnum = residue.get_id()[1]
                chain_id = residue.get_full_id()[2]
            except Exception:
                continue
            if (resnum, chain_id) not in allowed_residues:
                continue
            atoms_dicts = residue_atoms_as_dicts(residue)
            if not atoms_dicts:
                continue
            bonds = get_residue_bonds(atoms_dicts, bond_lookup=bond_lookup)
            full_bonds.extend(bonds)
        return full_bonds

    full_residue_bonds = _compute_full_residue_bonds_for_focused(structure, focused_atoms, bond_lookup_table)

    # ------------------------ Coordinates (markers) ------------------------
    all_atoms = focused_atoms
    all_coords = np.array([a["coordinates"] for a in all_atoms]) if all_atoms else np.empty((0, 3))
    cof_coords = np.array([a["coordinates"] for a in (cofactor_atoms or [])]) if cofactor_atoms else np.empty((0, 3))
    pcs_coords = np.array([a["coordinates"] for a in (pcs_atoms or [])]) if pcs_atoms else np.empty((0, 3))
    scs_coords = np.array([a["coordinates"] for a in (scs_atoms or [])]) if scs_atoms else np.empty((0, 3))

    pcs_scs_colors = []
    for coord in (all_coords if all_coords.size else []):
        if cof_coords.size and np.any(np.all(coord == cof_coords, axis=1)):
            pcs_scs_colors.append("black")
        elif pcs_coords.size and np.any(np.all(coord == pcs_coords, axis=1)):
            pcs_scs_colors.append("blue")
        elif scs_coords.size and np.any(np.all(coord == scs_coords, axis=1)):
            pcs_scs_colors.append("fuchsia")
        else:
            pcs_scs_colors.append("grey")

    element_colors = [elem_color(a) for a in all_atoms]

    # Build hover text (NOW WITH MOIETY)
    hover_text = [
        f"Name: {a.get('name','?')}<br>"
        f"Residue: {a.get('residue','?')} {a.get('residue_number','?')} {a.get('chain','?')}<br>"
        f"Element: {a.get('element','?')}<br>"
        f"Moiety: { _moiety_of(a) }"
        for a in all_atoms
    ] if all_atoms else []

    # ------------------------ Link segments (optional) ------------------------
    def _akey(a: Dict) -> tuple[str, int, str, str]:
        return (str(a.get("residue","")), int(a.get("residue_number", 0)), str(a.get("chain","")), str(a.get("name","")))

    atom_index: Dict[tuple[str,int,str,str], np.ndarray] = {}
    for a in all_atoms:
        atom_index[_akey(a)] = np.array(a["coordinates"], dtype=float)

    link_segments = []
    if links_rows is None and links_csv_path and os.path.isfile(links_csv_path):
        with open(links_csv_path, "r", newline="") as f:
            links_rows = list(csv.DictReader(f))
    links_rows = links_rows or []

    def _coord_from_row(prefix: str, row: Dict[str,str]) -> Optional[np.ndarray]:
        key = (
            str(row[f"{prefix}_resname"]),
            int(row[f"{prefix}_resnum"]),
            str(row[f"{prefix}_chain"]),
            str(row[f"{prefix}_atom"]),
        )
        return atom_index.get(key)

    missing = 0
    for r in links_rows:
        s = _coord_from_row("src", r)
        d = _coord_from_row("dst", r)
        if s is None or d is None:
            missing += 1
            continue
        link_segments.append((s[0],s[1],s[2], d[0],d[1],d[2], r.get("link_type","link")))
    if missing:
        logger.info("Link rendering: %d links skipped (atoms not present in current marker set).", missing)

    # ------------------------ Figure & traces ------------------------
    fig = go.Figure()

    # Atom markers (single trace)
    if all_coords.size:
        fig.add_trace(
            go.Scatter3d(
                x=all_coords[:, 0], y=all_coords[:, 1], z=all_coords[:, 2],
                mode="markers",
                marker=dict(size=8, color=pcs_scs_colors), # SIZE GOES HERE FOR COORD
                hoverinfo="text",
                text=hover_text,  # <--- includes Moiety now
                showlegend=False,
                name="Atoms",
            )
        )
    else:
        fig.add_trace(go.Scatter3d(x=[], y=[], z=[], mode="markers", marker=dict(size=4), showlegend=False, name="Atoms"))

    def _add_bond_traces(bonds, color, width, visible):
        for bond in bonds:
            try:
                xs, ys, zs = zip(*bond)
            except Exception:
                p1, p2 = bond
                xs, ys, zs = [p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]]
            fig.add_trace(go.Scatter3d(
                x=xs, y=ys, z=zs, mode="lines",
                line=dict(color=color, width=width),
                hoverinfo="skip", showlegend=False, visible=visible
            ))

    # 1) Minimal focused bonds (default visible)
    start_min = len(fig.data)
    _add_bond_traces(minimal_focused_bonds, color="black", width=4, visible=True)
    end_min = len(fig.data)

    # 2) Full residue bonds (initially hidden)
    start_full = len(fig.data)
    _add_bond_traces(full_residue_bonds, color="gray", width=4, visible=False)
    end_full = len(fig.data)

    # # 3) Dotted link traces (always visible)
    # start_links = len(fig.data)
    # for (x1,y1,z1, x2,y2,z2, ltype) in link_segments:
    #     color = "#444" if ltype == "cofactor->pcs" else "#888"
    #     fig.add_trace(go.Scatter3d(
    #         x=[x1, x2], y=[y1, y2], z=[z1, z2],
    #         mode="lines",
    #         line=dict(color=color, width=4, dash="dot"),
    #         hoverinfo="skip",
    #         showlegend=False,
    #         visible=True,
    #         name=ltype,
    #     ))
    # end_links = len(fig.data)


    # 3) Dotted link traces (always visible; lighter & less dense)
    start_links = len(fig.data)
    for (x1, y1, z1, x2, y2, z2, ltype) in link_segments:
        fig.add_trace(
            go.Scatter3d(
                x=[x1, x2], y=[y1, y2], z=[z1, z2],
                mode="lines",
                line=dict(
                    color="grey",       # softer neutral
                    width=2,            # lighter stroke
                    dash="longdashdot"  # wide dash pattern with sparse dots
                ),
                hoverinfo="skip",
                showlegend=False,
                visible=True,
                name=ltype,
            )
        )
    end_links = len(fig.data)




    # ------------------------ UI controls ------------------------
    def _vis_backbone(on: bool):
        vis = [True]                                    # atoms
        vis += ([not on] * (end_min - start_min))       # minimal
        vis += ([on] * (end_full - start_full))         # full
        vis += ([True] * (end_links - start_links))     # links on
        return vis

    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                buttons=[
                    dict(label="By Coordination Sphere", method="restyle", args=[{"marker.color": [pcs_scs_colors]}]),
                    dict(label="By Element", method="restyle", args=[{"marker.color": [element_colors]}]),
                ],
                direction="right", showactive=True, x=0.05, y=1.15, xanchor="left", yanchor="top",
                bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2, font=dict(size=16, color="white"),
            ),
            dict(
                type="buttons",
                buttons=[
                    dict(label="Background On", method="relayout", args=[{
                        "scene.xaxis.visible": True, "scene.yaxis.visible": True, "scene.zaxis.visible": True,
                        "scene.xaxis.showgrid": True, "scene.yaxis.showgrid": True, "scene.zaxis.showgrid": True,
                        "scene.backgroundcolor": "rgba(240,240,240,1)",
                    }]),
                    dict(label="Background Off", method="relayout", args=[{
                        "scene.xaxis.visible": False, "scene.yaxis.visible": False, "scene.zaxis.visible": False,
                        "scene.backgroundcolor": "rgba(255,255,255,1)",
                    }]),
                ],
                direction="right", showactive=True, x=0.05, y=1.05, xanchor="left", yanchor="top",
                bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2, font=dict(size=16, color="white"),
            ),
            dict(
                type="buttons",
                buttons=[
                    dict(label="Backbone Atoms: Off", method="update", args=[{"visible": _vis_backbone(on=False)}]),
                    dict(label="Backbone Atoms: On", method="update", args=[{"visible": _vis_backbone(on=True)}]),
                ],
                direction="right", showactive=True, x=0.05, y=0.95, xanchor="left", yanchor="top",
                bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2, font=dict(size=16, color="white"),
            ),
        ],
        title={"text": f"{cofactor_resname} in {pdb_name}", "x": 0.5, "font": {"size": 22}},
        scene=dict(xaxis_title="X", yaxis_title="Y", zaxis_title="Z"),
        margin=dict(l=0, r=0, t=60, b=0),
    )

    fig.write_html(output_filename)
    logger.info("Interactive plot saved as '%s'", output_filename)
    if show:
        show_plotly(fig)


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

    full_index, network_index, coordinator_index, context_index = 0, 1, 2, 3
    network_motif_index, network_backbone_index, focused_motif_index, focused_backbone_index = 4, 5, 6, 7
    cofactor_bond_index, focused_bond_index, residue_bond_index, full_bond_index, motif_bond_index, backbone_bond_index = 8, 9, 10, 11, 12, 13
    primary_index, secondary_index, tertiary_index = 14, 15, 16
    all_trace_indices = list(range(17))
    network_preset_visibility = [False, True, True, False, False, False, False, False, True, False, True, False, False, False, True, True, True]
    network_motifs_preset_visibility = [False, True, True, False, True, False, False, False, True, False, False, False, True, False, True, True, True]
    structural_context_preset_visibility = [False, True, True, True, False, False, True, True, True, False, False, False, False, True, True, True, True]
    full_protein_preset_visibility = [True, True, True, False, False, False, False, False, True, False, False, True, False, False, True, True, True]
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
                dict(label="Network", method="update", args=[{"visible": network_preset_visibility}, {"annotations[8].text": "View: Network"}]),
                dict(label="Network + motifs", method="update", args=[{"visible": network_motifs_preset_visibility}, {"annotations[8].text": "View: Network + motifs"}]),
                dict(label="Structural context", method="update", args=[{"visible": structural_context_preset_visibility}, {"annotations[8].text": "View: Structural context"}]),
                dict(label="Full protein context", method="update", args=[{"visible": full_protein_preset_visibility}, {"annotations[8].text": "View: Full protein context"}]),
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
            dict(type="dropdown", direction="down", x=0.01, y=0.40, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=12), pad=dict(l=5, r=5, t=3, b=3), buttons=[
                dict(label="All", method="restyle", args=[{"visible": [True, True, True]}, [primary_index, secondary_index, tertiary_index]]),
                dict(label="Primary", method="restyle", args=[{"visible": [True, False, False]}, [primary_index, secondary_index, tertiary_index]]),
                dict(label="Secondary", method="restyle", args=[{"visible": [False, True, False]}, [primary_index, secondary_index, tertiary_index]]),
                dict(label="Tertiary", method="restyle", args=[{"visible": [False, False, True]}, [primary_index, secondary_index, tertiary_index]]),
                dict(label="None", method="restyle", args=[{"visible": [False, False, False]}, [primary_index, secondary_index, tertiary_index]]),
            ]),
            dict(type="dropdown", direction="down", x=0.01, y=0.29, xanchor="left", yanchor="top", showactive=True,
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
            dict(type="dropdown", direction="up", x=0.01, y=0.18, xanchor="left", yanchor="top", showactive=True,
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
            dict(type="dropdown", direction="up", x=0.01, y=0.07, xanchor="left", yanchor="top", showactive=True,
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
            dict(text="Contacts", x=0.01, y=0.445, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Grid", x=0.01, y=0.335, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Background", x=0.01, y=0.225, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Labels", x=0.01, y=0.115, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
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
    fig.write_html(
        output_filename,
        config=viewer_config,
        include_plotlyjs="cdn" if compact_html else True,
        post_script=post_script,
        div_id="coordination-network-plot",
    )
    logger.info("Cohesive interactive plot saved as '%s'", output_filename)



# def plot_interactive_modes_with_network(
#     structure,
#     cofactor_atoms: List[Dict],
#     pcs_atoms: List[Dict],
#     scs_atoms: List[Dict],
#     bond_lookup_table: Dict[str, List[List[str]]],
#     pdb_name: str = "structure.pdb",
#     cofactor_resname: str = "cofactor",
#     atom_type_colors: Optional[Dict[str, str]] = None,
#     output_filename: str = "1_template_coordination_network.html",
#     # NEW: provide either a CSV path or a preloaded list[dict] with link rows
#     links_csv_path: Optional[str] = "Coord_Links.csv",
#     links_rows: Optional[List[Dict[str, str]]] = None,
# ):
#     """
#     Interactive 3D Plotly viz with:
#       • Coloring toggle: Coordination Sphere / Element
#       • Backbone sticks toggle: Off (focused) / On (full residue sticks for residues containing any shown atom)
#       • NEW: Dotted link lines for cofactor→PCS and PCS→SCS based on Coord_Links.csv (or provided rows)
#     Links do not add atom markers; they’re just polylines drawn atop the scene.
#     """
#     import os
#     import csv
#     import numpy as np
#     import plotly.graph_objects as go

#     # ------------------------ Colors ------------------------
#     if atom_type_colors is None:
#         atom_type_colors = {
#             "C": "black", "N": "blue", "O": "red", "S": "yellow",
#             "FE": "orange", "MN": "purple", "CA": "green", "CU": "goldenrod",
#             "MO": "teal",
#         }
#     default_color = "grey"

#     def elem_color(atom: Dict) -> str:
#         return atom_type_colors.get(str(atom.get("element", "")).upper(), default_color)

#     # ------------------------ Bond builders ------------------------
#     focused_atoms = (cofactor_atoms or []) + (pcs_atoms or []) + (scs_atoms or [])
#     minimal_focused_bonds = generate_residue_bonds(focused_atoms, bond_lookup_table)

#     def _compute_full_residue_bonds_for_focused(structure, focused_atoms, bond_lookup):
#         allowed_residues = {(a.get("residue_number", None), a.get("chain", None))
#                             for a in focused_atoms if a is not None}
#         full_bonds = []

#         def residue_atoms_as_dicts(residue):
#             rname = residue.get_resname()
#             if not rname:
#                 return []
#             out = []
#             for at in residue:
#                 out.append({
#                     "name": at.get_name(),
#                     "element": getattr(at, "element", ""),
#                     "coordinates": np.array(at.coord, dtype=float),
#                     "residue": rname,
#                 })
#             return out

#         for residue in structure.get_residues():
#             try:
#                 resnum = residue.get_id()[1]
#                 chain_id = residue.get_full_id()[2]
#             except Exception:
#                 continue
#             if (resnum, chain_id) not in allowed_residues:
#                 continue
#             atoms_dicts = residue_atoms_as_dicts(residue)
#             if not atoms_dicts:
#                 continue
#             bonds = get_residue_bonds(atoms_dicts, bond_lookup=bond_lookup)
#             full_bonds.extend(bonds)
#         return full_bonds

#     full_residue_bonds = _compute_full_residue_bonds_for_focused(structure, focused_atoms, bond_lookup_table)



#     # ------------------------ Coordinates (markers) ------------------------
#     all_atoms = focused_atoms
#     all_coords = np.array([a["coordinates"] for a in all_atoms]) if all_atoms else np.empty((0, 3))
#     cof_coords = np.array([a["coordinates"] for a in (cofactor_atoms or [])]) if cofactor_atoms else np.empty((0, 3))
#     pcs_coords = np.array([a["coordinates"] for a in (pcs_atoms or [])]) if pcs_atoms else np.empty((0, 3))
#     scs_coords = np.array([a["coordinates"] for a in (scs_atoms or [])]) if scs_atoms else np.empty((0, 3))

#     pcs_scs_colors = []
#     for coord in (all_coords if all_coords.size else []):
#         if cof_coords.size and np.any(np.all(coord == cof_coords, axis=1)):
#             pcs_scs_colors.append("black")
#         elif pcs_coords.size and np.any(np.all(coord == pcs_coords, axis=1)):
#             pcs_scs_colors.append("blue")
#         elif scs_coords.size and np.any(np.all(coord == scs_coords, axis=1)):
#             pcs_scs_colors.append("fuchsia")
#         else:
#             pcs_scs_colors.append("grey")

#     element_colors = [elem_color(a) for a in all_atoms]

#     # ------------------------ Build atom index for link lookup ------------------------
#     # Keyed by (resname,resnum,chain,atom)
#     def _akey(a: Dict) -> tuple[str, int, str, str]:
#         return (str(a.get("residue","")), int(a.get("residue_number", 0)), str(a.get("chain","")), str(a.get("name","")))

#     atom_index: Dict[tuple[str,int,str,str], np.ndarray] = {}
#     for a in all_atoms:
#         atom_index[_akey(a)] = np.array(a["coordinates"], dtype=float)

#     # ------------------------ Load links (CSV or provided rows) ------------------------
#     link_segments = []  # each: (x1,y1,z1,x2,y2,z2, link_type)
#     def _try_get_coord(row_prefix: str, row: Dict[str,str]) -> Optional[np.ndarray]:
#         key = (
#             str(row[f"{row_prefix}_resname"]),
#             int(row[f"{row_prefix}_resnum"]),
#             str(row[f"{row_prefix}_chain"]),
#             str(row[f"{row_prefix}_atom"]),
#         )
#         return atom_index.get(key)

#     if links_rows is None and links_csv_path and os.path.isfile(links_csv_path):
#         with open(links_csv_path, "r", newline="") as f:
#             links_rows = list(csv.DictReader(f))
#     links_rows = links_rows or []

#     missing = 0
#     for r in links_rows:
#         src = _try_get_coord("src", r)
#         dst = _try_get_coord("dst", r)
#         if src is None or dst is None:
#             missing += 1
#             continue
#         link_segments.append((src[0],src[1],src[2], dst[0],dst[1],dst[2], r.get("link_type","link")))
#     if missing:
#         print(f"[INFO] Link rendering: {missing} links skipped (atoms not present in current marker set).")

#     # ------------------------ Figure & traces ------------------------
#     fig = go.Figure()

#     # 0) Atom markers (single trace)
#     if all_coords.size:
#         fig.add_trace(
#             go.Scatter3d(
#                 x=all_coords[:, 0], y=all_coords[:, 1], z=all_coords[:, 2],
#                 mode="markers",
#                 marker=dict(size=4, color=pcs_scs_colors),
#                 hoverinfo="text",
#                 text=[
#                     f"Name: {a.get('name','?')}<br>"
#                     f"Residue: {a.get('residue','?')} {a.get('residue_number','?')} {a.get('chain','?')}<br>"
#                     f"Element: {a.get('element','?')}"
#                     for a in all_atoms
#                 ],
#                 showlegend=False,
#                 name="Atoms",
#             )
#         )
#     else:
#         fig.add_trace(go.Scatter3d(x=[], y=[], z=[], mode="markers", marker=dict(size=4), showlegend=False, name="Atoms"))

#     # Helper to add bond traces
#     def _add_bond_traces(bonds, color, width, visible):
#         for bond in bonds:
#             try:
#                 xs, ys, zs = zip(*bond)
#             except Exception:
#                 p1, p2 = bond
#                 xs, ys, zs = [p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]]
#             fig.add_trace(go.Scatter3d(
#                 x=xs, y=ys, z=zs, mode="lines",
#                 line=dict(color=color, width=width),
#                 hoverinfo="skip", showlegend=False, visible=visible
#             ))

#     # 1) Minimal focused bonds (default visible)
#     start_min = len(fig.data)
#     _add_bond_traces(minimal_focused_bonds, color="black", width=2, visible=True)
#     end_min = len(fig.data)

#     # 2) Full residue bonds (initially hidden)
#     start_full = len(fig.data)
#     _add_bond_traces(full_residue_bonds, color="gray", width=2, visible=False)
#     end_full = len(fig.data)

#     # 3) NEW — Dotted link traces (always visible)
#     start_links = len(fig.data)
#     for (x1,y1,z1, x2,y2,z2, ltype) in link_segments:
#         # color by link type: cofactor->pcs darker; pcs->scs lighter
#         color = "#444" if ltype == "cofactor->pcs" else "#888"
#         fig.add_trace(go.Scatter3d(
#             x=[x1, x2], y=[y1, y2], z=[z1, z2],
#             mode="lines",
#             line=dict(color=color, width=4, dash="dot"),
#             hoverinfo="skip",
#             showlegend=False,
#             visible=True,
#             name=ltype,
#         ))
#     end_links = len(fig.data)

#     # ------------------------ UI controls ------------------------
#     def _vis_backbone(on: bool):
#         # atoms (1) + minimal bonds + full bonds + links
#         vis = [True]                                      # atoms
#         vis += ([not on] * (end_min - start_min))         # minimal
#         vis += ([on] * (end_full - start_full))           # full
#         vis += ([True] * (end_links - start_links))       # links always on
#         return vis

#     fig.update_layout(
#         updatemenus=[
#             dict(
#                 type="buttons",
#                 buttons=[
#                     dict(label="By Coordination Sphere", method="restyle", args=[{"marker.color": [pcs_scs_colors]}]),
#                     dict(label="By Element", method="restyle", args=[{"marker.color": [element_colors]}]),
#                 ],
#                 direction="right", showactive=True, x=0.05, y=1.15, xanchor="left", yanchor="top",
#                 bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2, font=dict(size=16, color="white"),
#             ),
#             dict(
#                 type="buttons",
#                 buttons=[
#                     dict(label="Background On", method="relayout", args=[{
#                         "scene.xaxis.visible": True, "scene.yaxis.visible": True, "scene.zaxis.visible": True,
#                         "scene.xaxis.showgrid": True, "scene.yaxis.showgrid": True, "scene.zaxis.showgrid": True,
#                         "scene.backgroundcolor": "rgba(240,240,240,1)",
#                     }]),
#                     dict(label="Background Off", method="relayout", args=[{
#                         "scene.xaxis.visible": False, "scene.yaxis.visible": False, "scene.zaxis.visible": False,
#                         "scene.backgroundcolor": "rgba(255,255,255,1)",
#                     }]),
#                 ],
#                 direction="right", showactive=True, x=0.05, y=1.05, xanchor="left", yanchor="top",
#                 bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2, font=dict(size=16, color="white"),
#             ),
#             dict(
#                 type="buttons",
#                 buttons=[
#                     dict(label="Backbone Atoms: Off", method="update", args=[{"visible": _vis_backbone(on=False)}]),
#                     dict(label="Backbone Atoms: On", method="update", args=[{"visible": _vis_backbone(on=True)}]),
#                 ],
#                 direction="right", showactive=True, x=0.05, y=0.95, xanchor="left", yanchor="top",
#                 bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2, font=dict(size=16, color="white"),
#             ),
#         ],
#         title={"text": f"{cofactor_resname} in {pdb_name}", "x": 0.5, "font": {"size": 22}},
#         scene=dict(xaxis_title="X", yaxis_title="Y", zaxis_title="Z"),
#         margin=dict(l=0, r=0, t=60, b=0),
#     )

#     fig.write_html(output_filename)
#     print(f"Interactive plot saved as '{output_filename}'")
#     fig.show()




# def plot_interactive_modes_with_network(
#     structure,
#     cofactor_atoms: List[Dict],
#     pcs_atoms: List[Dict],
#     scs_atoms: List[Dict],
#     bond_lookup_table: Dict[str, List[List[str]]],
#     pdb_name: str = "structure.pdb",
#     cofactor_resname: str = "cofactor",
#     atom_type_colors: Optional[Dict[str, str]] = None,
#     output_filename: str = "1_template_coordination_network.html",
# ):
#     """
#     Interactive 3D Plotly visualization with two coloring modes:
#       - By coordination sphere (cofactor/PCS/SCS)
#       - By element (C/N/O/S/…)
#     New toggle:
#       - Backbone Atoms: Off (default) -> minimal 'focused' bonds using only the provided atoms
#       - Backbone Atoms: On            -> complete stick connectivity for any residue that contains
#                                          at least one provided cofactor/PCS/SCS atom (includes backbone)
#     Note: The 'On' mode does NOT add new atom markers; it only adds stick lines for those residues.
#     """
#     import numpy as np
#     import plotly.graph_objects as go

#     # ------------------------ Colors ------------------------
#     if atom_type_colors is None:
#         atom_type_colors = {
#             "C": "black",
#             "N": "blue",
#             "O": "red",
#             "S": "yellow",
#             "FE": "orange",
#             "MN": "purple",
#             "CA": "green",
#             "CU": "goldenrod",
#             "MO": "teal",
#         }
#     default_color = "grey"

#     def elem_color(atom: Dict) -> str:
#         return atom_type_colors.get(str(atom.get("element", "")).upper(), default_color)

#     # ------------------------ Bond builders ------------------------
#     # Minimal bonds among only the provided atoms (your prior focused state)
#     # Assumes you already have this helper elsewhere.
#     # If it returns a list of bonds where each bond is a list of 2 points [(x,y,z), (x,y,z)], perfect.
#     focused_atoms = (cofactor_atoms or []) + (pcs_atoms or []) + (scs_atoms or [])
#     minimal_focused_bonds = generate_residue_bonds(focused_atoms, bond_lookup_table)

#     # Helper: for residues that contain any provided atom, compute FULL residue bonds
#     # using the bond_lookup_table (includes backbone sticks). This does not add atom markers.
#     def _compute_full_residue_bonds_for_focused(structure, focused_atoms, bond_lookup):
#         # Gather residues (resnum, chain) that contain at least one focused atom
#         allowed_residues = {(a.get("residue_number", None), a.get("chain", None))
#                             for a in focused_atoms if a is not None}

#         full_bonds = []

#         # Utility: convert a Biopython residue to list[dict] for get_residue_bonds(...)
#         def residue_atoms_as_dicts(residue):
#             rname = residue.get_resname()
#             if rname is None:
#                 return []
#             out = []
#             for at in residue:
#                 out.append({
#                     "name": at.get_name(),
#                     "element": getattr(at, "element", ""),
#                     "coordinates": np.array(at.coord, dtype=float),
#                     "residue": rname,
#                 })
#             return out

#         # We rely on a per-residue bond function (commonly named get_residue_bonds)
#         # which should accept (atoms_as_dicts, bond_lookup=...) and return list of bonds.
#         for residue in structure.get_residues():
#             try:
#                 resnum = residue.get_id()[1]
#             except Exception:
#                 continue
#             try:
#                 chain_id = residue.get_full_id()[2]
#             except Exception:
#                 chain_id = None

#             if (resnum, chain_id) not in allowed_residues:
#                 continue

#             atoms_dicts = residue_atoms_as_dicts(residue)
#             if not atoms_dicts:
#                 continue

#             # IMPORTANT: use the same per-residue bonding function you use elsewhere
#             bonds = get_residue_bonds(atoms_dicts, bond_lookup=bond_lookup)
#             full_bonds.extend(bonds)

#         return full_bonds

#     full_residue_bonds = _compute_full_residue_bonds_for_focused(structure, focused_atoms, bond_lookup_table)

#     # ------------------------ Coordinates (markers) ------------------------
#     all_atoms = focused_atoms
#     all_coords = np.array([a["coordinates"] for a in all_atoms]) if all_atoms else np.empty((0, 3))
#     cof_coords = np.array([a["coordinates"] for a in (cofactor_atoms or [])]) if cofactor_atoms else np.empty((0, 3))
#     pcs_coords = np.array([a["coordinates"] for a in (pcs_atoms or [])]) if pcs_atoms else np.empty((0, 3))
#     scs_coords = np.array([a["coordinates"] for a in (scs_atoms or [])]) if scs_atoms else np.empty((0, 3))

#     # Coordination-sphere coloring (kept as-is with coordinate membership check)
#     pcs_scs_colors = []
#     for coord in (all_coords if all_coords.size else []):
#         if cof_coords.size and np.any(np.all(coord == cof_coords, axis=1)):
#             pcs_scs_colors.append("black")
#         elif pcs_coords.size and np.any(np.all(coord == pcs_coords, axis=1)):
#             pcs_scs_colors.append("blue")
#         elif scs_coords.size and np.any(np.all(coord == scs_coords, axis=1)):
#             pcs_scs_colors.append("fuchsia")
#         else:
#             pcs_scs_colors.append("grey")

#     # Element coloring
#     element_colors = [elem_color(a) for a in all_atoms]

#     # ------------------------ Figure & traces ------------------------
#     fig = go.Figure()

#     # Atom markers (single trace)
#     if all_coords.size:
#         fig.add_trace(
#             go.Scatter3d(
#                 x=all_coords[:, 0],
#                 y=all_coords[:, 1],
#                 z=all_coords[:, 2],
#                 mode="markers",
#                 marker=dict(size=4, color=pcs_scs_colors),
#                 hoverinfo="text",
#                 text=[
#                     f"Name: {a.get('name','?')}<br>"
#                     f"Residue: {a.get('residue','?')}<br>"
#                     f"Residue Number: {a.get('residue_number','?')}<br>"
#                     f"Chain: {a.get('chain','?')}<br>"
#                     f"Element: {a.get('element','?')}"
#                     for a in all_atoms
#                 ],
#                 showlegend=False,
#                 name="Atoms",
#             )
#         )
#     else:
#         # keep an empty atom trace so UI always works
#         fig.add_trace(go.Scatter3d(x=[], y=[], z=[], mode="markers", marker=dict(size=4), showlegend=False, name="Atoms"))

#     # Helper to add bond traces (one trace per bond/polyline)
#     def _add_bond_traces(bonds, color, width, visible):
#         for bond in bonds:
#             try:
#                 xs, ys, zs = zip(*bond)
#             except Exception:
#                 # if bond is just two points
#                 p1, p2 = bond
#                 xs, ys, zs = [p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]]
#             fig.add_trace(
#                 go.Scatter3d(
#                     x=xs, y=ys, z=zs,
#                     mode="lines",
#                     line=dict(color=color, width=width),
#                     hoverinfo="skip",
#                     showlegend=False,
#                     visible=visible,
#                 )
#             )

#     # 1) Minimal focused bonds (default visible) -> "Backbone Atoms: Off"
#     start_min = len(fig.data)
#     _add_bond_traces(minimal_focused_bonds, color="black", width=2, visible=True)
#     end_min = len(fig.data)

#     # 2) Full residue bonds for any residue containing focused atoms (initially hidden) -> "Backbone Atoms: On"
#     start_full = len(fig.data)
#     _add_bond_traces(full_residue_bonds, color="gray", width=2, visible=False)
#     end_full = len(fig.data)

#     # ------------------------ UI controls ------------------------
#     # Visibility arrays for the toggle (atoms always True)
#     def _vis_backbone(on: bool):
#         vis = [True]  # atom markers
#         # minimal range
#         vis += ([not on] * (end_min - start_min))
#         # full range
#         vis += ([on] * (end_full - start_full))
#         return vis

#     fig.update_layout(
#         updatemenus=[
#             # Coloring mode
#             dict(
#                 type="buttons",
#                 buttons=[
#                     dict(
#                         label="By Coordination Sphere",
#                         method="restyle",
#                         args=[{"marker.color": [pcs_scs_colors]}],
#                     ),
#                     dict(
#                         label="By Element",
#                         method="restyle",
#                         args=[{"marker.color": [element_colors]}],
#                     ),
#                 ],
#                 direction="right",
#                 showactive=True,
#                 x=0.05,
#                 y=1.15,
#                 xanchor="left",
#                 yanchor="top",
#                 bgcolor="rgba(50, 50, 50, 0.8)",
#                 bordercolor="black",
#                 borderwidth=2,
#                 font=dict(size=16, color="white"),
#             ),
#             # Background toggle
#             dict(
#                 type="buttons",
#                 buttons=[
#                     dict(
#                         label="Background On",
#                         method="relayout",
#                         args=[{
#                             "scene.xaxis.visible": True,
#                             "scene.yaxis.visible": True,
#                             "scene.zaxis.visible": True,
#                             "scene.xaxis.showgrid": True,
#                             "scene.yaxis.showgrid": True,
#                             "scene.zaxis.showgrid": True,
#                             "scene.backgroundcolor": "rgba(240,240,240,1)",
#                         }],
#                     ),
#                     dict(
#                         label="Background Off",
#                         method="relayout",
#                         args=[{
#                             "scene.xaxis.visible": False,
#                             "scene.yaxis.visible": False,
#                             "scene.zaxis.visible": False,
#                             "scene.backgroundcolor": "rgba(255,255,255,1)",
#                         }],
#                     ),
#                 ],
#                 direction="right",
#                 showactive=True,
#                 x=0.05,
#                 y=1.05,
#                 xanchor="left",
#                 yanchor="top",
#                 bgcolor="rgba(50, 50, 50, 0.8)",
#                 bordercolor="black",
#                 borderwidth=2,
#                 font=dict(size=16, color="white"),
#             ),
#             # Backbone Atoms toggle (replaces prior bond visibility buttons)
#             dict(
#                 type="buttons",
#                 buttons=[
#                     dict(
#                         label="Backbone Atoms: Off",
#                         method="update",
#                         args=[{"visible": _vis_backbone(on=False)}],
#                     ),
#                     dict(
#                         label="Backbone Atoms: On",
#                         method="update",
#                         args=[{"visible": _vis_backbone(on=True)}],
#                     ),
#                 ],
#                 direction="right",
#                 showactive=True,
#                 x=0.05,
#                 y=0.95,
#                 xanchor="left",
#                 yanchor="top",
#                 bgcolor="rgba(50, 50, 50, 0.8)",
#                 bordercolor="black",
#                 borderwidth=2,
#                 font=dict(size=16, color="white"),
#             ),
#         ],
#         title={"text": f"{cofactor_resname} in {pdb_name}", "x": 0.5, "font": {"size": 22}},
#         scene=dict(xaxis_title="X", yaxis_title="Y", zaxis_title="Z"),
#         margin=dict(l=0, r=0, t=60, b=0),
#     )

#     # ------------------------ Save & show ------------------------
#     fig.write_html(output_filename)
#     print(f"Interactive plot saved as '{output_filename}'")
#     fig.show()

















# --- paste anywhere in modules/plotting.py (e.g., after plot_interactive_modes_with_network) ---
def plot_interactive_modes_with_roi(
    structure, 
    cofactor_atoms, 
    roi_atoms, 
    bond_lookup_table, 
    pdb_name="structure.pdb", 
    cofactor_resname="cofactor",
    atom_type_colors: dict = None, 
    output_filename: str = "residues_of_interest_coordination_network.html"
):
    """
    Interactive Plotly 3D viewer that toggles coloring between:
      - Residues of Interest (ROI) vs Cofactor
      - Element-based coloring
    with buttons to show focused bonds (cofactor+ROI) or the whole-protein bonds.
    cofactor_atoms and roi_atoms are lists of dicts with keys:
      'coordinates', 'name', 'element', 'residue', 'residue_number', 'chain'
    """
    if atom_type_colors is None:
        atom_type_colors = {
            "C": "black",  # Carbon
            "N": "blue",   # Nitrogen
            "O": "red",    # Oxygen
            "S": "yellow", # Sulfur
            "FE": "orange",
            "MN": "purple",
            "CA": "green",
            "CU": "goldenrod",
            "MO": "teal",
        }
    default_color = "grey"

    def get_atom_coordinate(atom):
        return atom['coordinates']

    def get_atom_element(atom):
        return atom['element']

    def get_atom_name(atom):
        return atom['name']

    def get_atom_residue(atom):
        return atom['residue']

    def get_atom_residue_number(atom):
        return atom['residue_number']

    def get_atom_chain(atom):
        return atom['chain']

    def get_element_color(atom):
        return atom_type_colors.get(str(get_atom_element(atom)).upper(), default_color)

    # Bonds
    all_bonds = generate_all_bonds(structure, bond_lookup_table)
    focused_atoms = cofactor_atoms + roi_atoms
    focused_bonds = generate_residue_bonds(focused_atoms, bond_lookup_table)

    # Coordinates
    all_coords = np.array([get_atom_coordinate(a) for a in focused_atoms], dtype=float)
    cofactor_coords = np.array([get_atom_coordinate(a) for a in cofactor_atoms], dtype=float) if cofactor_atoms else np.empty((0,3))
    roi_coords = np.array([get_atom_coordinate(a) for a in roi_atoms], dtype=float) if roi_atoms else np.empty((0,3))

    # Fixed (ROI/cofactor) coloring
    fixed_colors = []
    for coord in all_coords:
        if roi_coords.size and np.any(np.all(coord == roi_coords, axis=1)):
            fixed_colors.append('fuchsia')
        elif cofactor_coords.size and np.any(np.all(coord == cofactor_coords, axis=1)):
            fixed_colors.append('black')
        else:
            fixed_colors.append('grey')

    # Element coloring
    element_colors = [get_element_color(a) for a in focused_atoms]

    # Hover text
    hover_text = [
        f"Name: {get_atom_name(a)}<br>Residue: {get_atom_residue(a)}<br>"
        f"Residue Number: {get_atom_residue_number(a)}<br>Element: {get_atom_element(a)}"
        for a in focused_atoms
    ]

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=all_coords[:, 0], y=all_coords[:, 1], z=all_coords[:, 2],
        mode='markers',
        marker=dict(size=4, color=fixed_colors),
        hoverinfo='text', text=hover_text, showlegend=False
    ))

    # Add bonds
    def add_bonds(bonds, visible):
        for bond in bonds:
            x_coords, y_coords, z_coords = zip(*bond)
            fig.add_trace(go.Scatter3d(
                x=x_coords, y=y_coords, z=z_coords,
                mode='lines', line=dict(color='black', width=2),
                hoverinfo='skip', showlegend=False, visible=visible
            ))

    add_bonds(focused_bonds, True)
    add_bonds(all_bonds, False)

    num_focused_bonds = len(focused_bonds)
    num_all_bonds = len(all_bonds)
    focused_bonds_visible = [True] * num_focused_bonds + [False] * num_all_bonds
    all_bonds_visible = [False] * num_focused_bonds + [True] * num_all_bonds

    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                buttons=[
                    dict(label="By Residues of Interest",
                         method="restyle",
                         args=[{"marker.color": [fixed_colors]},
                               {"title": f"Residues of Interest Coloring for {cofactor_resname} in {pdb_name}"}]),
                    dict(label="By Element",
                         method="restyle",
                         args=[{"marker.color": [element_colors]},
                               {"title": f"Element Coloring for {cofactor_resname} in {pdb_name}"}]),
                ],
                direction="right", showactive=True, active=0,
                x=0.05, y=1.15, xanchor="left", yanchor="top",
                bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2,
                font=dict(size=20, color="white")
            ),
            dict(
                type="buttons",
                buttons=[
                    dict(label="Background On", method="relayout",
                         args=[{"scene.xaxis.visible": True, "scene.yaxis.visible": True, "scene.zaxis.visible": True,
                                "scene.xaxis.showgrid": True, "scene.yaxis.showgrid": True, "scene.zaxis.showgrid": True,
                                "scene.backgroundcolor": "rgba(240,240,240,1)"}]),
                    dict(label="Background Off", method="relayout",
                         args=[{"scene.xaxis.visible": False, "scene.yaxis.visible": False, "scene.zaxis.visible": False,
                                "scene.backgroundcolor": "rgba(255,255,255,1)"}]),
                ],
                direction="right", showactive=True, active=0,
                x=0.05, y=1.05, xanchor="left", yanchor="top",
                bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2,
                font=dict(size=20, color="white")
            ),
            dict(
                type="buttons",
                buttons=[
                    dict(label="Focused Bonds", method="update",
                         args=[{"visible": [True] + focused_bonds_visible},
                               {"title": f"Focused Bonds for {cofactor_resname} in {pdb_name}"}]),
                    dict(label="Whole Protein Bonds", method="update",
                         args=[{"visible": [True] + all_bonds_visible},
                               {"title": f"Whole Protein Bonds for {cofactor_resname} in {pdb_name}"}]),
                ],
                direction="right", showactive=True, active=0,
                x=0.05, y=0.95, xanchor="left", yanchor="top",
                bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2,
                font=dict(size=20, color="white")
            ),
        ],
        title={"text": f"{cofactor_resname} in {pdb_name}", "x": 0.5, "font": {"size": 24}},
        scene=dict(xaxis_title="X Coordinate", yaxis_title="Y Coordinate", zaxis_title="Z Coordinate"),
    )

    fig.write_html(output_filename)
    logger.info("Interactive plot saved as '%s'", output_filename)
    show_plotly(fig)







# --- paste somewhere below your other plotting funcs ---
def plot_template_heatmap_interactive(
    structure, 
    cofactor_atoms, pcs_atoms, scs_atoms,
    bond_lookup_table, 
    pdb_name="structure.pdb", 
    cofactor_resname="cofactor",
    highly_conserved_template_resnum_list=None,
    atom_type_colors: dict = None
):
    """
    Interactive Plotly 3D viewer that can color by:
      - coordination sphere (cofactor / PCS / SCS / other),
      - element type,
      - and optionally highlight a set of 'highly conserved' residues (dark orange).
    Also includes buttons to toggle focused vs whole-protein bonds and background grid.
    """
    # default element colors
    if atom_type_colors is None:
        atom_type_colors = {
            "C": "black", "N": "blue", "O": "red", "S": "yellow",
            "FE": "orange", "MN": "purple", "CA": "green", "CU": "goldenrod",
        }
    default_color = "grey"

    def get_element_color(atom):
        return atom_type_colors.get(str(atom["element"]).upper(), default_color)

    # Bonds
    all_bonds = generate_all_bonds(structure, bond_lookup_table)
    focused_atoms = cofactor_atoms + pcs_atoms + scs_atoms
    focused_bonds = generate_residue_bonds(focused_atoms, bond_lookup_table)

    # Coordinates
    all_atoms = focused_atoms
    all_coords = np.array([a['coordinates'] for a in all_atoms], dtype=float) if all_atoms else np.empty((0,3))
    cofactor_coords = np.array([a['coordinates'] for a in cofactor_atoms], dtype=float) if cofactor_atoms else np.empty((0,3))
    pcs_coords      = np.array([a['coordinates'] for a in pcs_atoms], dtype=float) if pcs_atoms else np.empty((0,3))
    scs_coords      = np.array([a['coordinates'] for a in scs_atoms], dtype=float) if scs_atoms else np.empty((0,3))

    if highly_conserved_template_resnum_list is not None:
        hc_set = set(highly_conserved_template_resnum_list)
        highly_conserved_coords = np.array(
            [a['coordinates'] for a in all_atoms if a.get("residue_number") in hc_set],
            dtype=float
        )
    else:
        highly_conserved_coords = np.empty((0,3))

    # Coord-sphere coloring (with HC override)
    pcs_scs_colors = []
    for coord in all_coords:
        if cofactor_coords.size and np.any(np.all(coord == cofactor_coords, axis=1)):
            pcs_scs_colors.append('black')
        elif highly_conserved_coords.size and np.any(np.all(coord == highly_conserved_coords, axis=1)):
            pcs_scs_colors.append('darkorange')
        elif pcs_coords.size and np.any(np.all(coord == pcs_coords, axis=1)):
            pcs_scs_colors.append('blue')
        elif scs_coords.size and np.any(np.all(coord == scs_coords, axis=1)):
            pcs_scs_colors.append('fuchsia')
        else:
            pcs_scs_colors.append('grey')

    # Element coloring
    element_colors = [get_element_color(a) for a in all_atoms]

    # Hover text
    hover_text = [
        f"Name: {a['name']}<br>Residue: {a['residue']}<br>"
        f"Residue Number: {a['residue_number']}<br>Element: {a['element']}"
        for a in all_atoms
    ]

    fig = go.Figure()
    if all_coords.size:
        fig.add_trace(go.Scatter3d(
            x=all_coords[:,0], y=all_coords[:,1], z=all_coords[:,2],
            mode='markers',
            marker=dict(size=4, color=pcs_scs_colors),
            hoverinfo='text', text=hover_text, showlegend=False
        ))

    # bonds as separate traces
    def add_bonds(bonds, visible=True):
        for bond in bonds:
            x, y, z = zip(*bond)
            fig.add_trace(go.Scatter3d(
                x=x, y=y, z=z, mode='lines',
                line=dict(color='black', width=2),
                hoverinfo='skip', showlegend=False, visible=visible
            ))
    add_bonds(focused_bonds, True)
    add_bonds(all_bonds, False)

    # visibility toggles for bond traces
    num_focused = len(focused_bonds)
    num_all     = len(all_bonds)
    focused_vis = [True] * num_focused + [False] * num_all
    all_vis     = [False] * num_focused + [True] * num_all

    fig.update_layout(
        updatemenus=[
            # coloring mode
            dict(
                type="buttons",
                buttons=[
                    dict(label="By Coordination Sphere",
                         method="restyle",
                         args=[{"marker.color": [pcs_scs_colors]},
                               {"title": f"PCS & SCS Coloring for {cofactor_resname} in {pdb_name}"}]),
                    dict(label="By Element",
                         method="restyle",
                         args=[{"marker.color": [element_colors]},
                               {"title": f"Element Coloring for {cofactor_resname} in {pdb_name}"}]),
                ],
                direction="right", showactive=True, active=0,
                x=0.05, y=1.15, xanchor="left", yanchor="top",
                bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2,
                font=dict(size=20, color="white"),
            ),
            # background
            dict(
                type="buttons",
                buttons=[
                    dict(label="Background On", method="relayout",
                         args=[{"scene.xaxis.visible": True, "scene.yaxis.visible": True, "scene.zaxis.visible": True,
                                "scene.xaxis.showgrid": True, "scene.yaxis.showgrid": True, "scene.zaxis.showgrid": True,
                                "scene.backgroundcolor": "rgba(240,240,240,1)"}]),
                    dict(label="Background Off", method="relayout",
                         args=[{"scene.xaxis.visible": False, "scene.yaxis.visible": False, "scene.zaxis.visible": False,
                                "scene.backgroundcolor": "rgba(255,255,255,1)"}]),
                ],
                direction="right", showactive=True, active=0,
                x=0.05, y=1.05, xanchor="left", yanchor="top",
                bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2,
                font=dict(size=20, color="white"),
            ),
            # bonds
            dict(
                type="buttons",
                buttons=[
                    dict(label="Focused Bonds", method="update",
                         args=[{"visible": [True] + focused_vis},
                               {"title": f"Focused Bonds for {cofactor_resname} in {pdb_name}"}]),
                    dict(label="Whole Protein Bonds", method="update",
                         args=[{"visible": [True] + all_vis},
                               {"title": f"Whole Protein Bonds for {cofactor_resname} in {pdb_name}"}]),
                ],
                direction="right", showactive=True, active=0,
                x=0.05, y=0.95, xanchor="left", yanchor="top",
                bgcolor="rgba(50,50,50,0.8)", bordercolor="black", borderwidth=2,
                font=dict(size=20, color="white"),
            ),
        ],
        title={"text": f"{cofactor_resname} in {pdb_name}", "x": 0.5, "font": {"size": 24}},
        scene=dict(xaxis_title="X Coordinate", yaxis_title="Y Coordinate", zaxis_title="Z Coordinate"),
    )

    fig.write_html("1_conserved_interactive_modes_with_bonds.html")
    logger.info("Interactive plot saved as '1_conserved_interactive_modes_with_bonds.html'")
    show_plotly(fig)
