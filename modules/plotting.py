# modules/plotting.py
"""
Lightweight plotting utilities.

Currently exposes:
- plot_evaluation_results(df_results, mode="CA_only", output_csv=None, highlight_labels=None)
"""

from typing import List, Optional
import pandas as pd
import matplotlib.pyplot as plt

from typing import List, Dict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (needed for 3D projection)
from .structure_processing import get_residue_bonds



# modules/plotting.py

from typing import Dict, List, Optional
import numpy as np
import plotly.graph_objects as go

# pull bond builders from structure_processing
from .structure_processing import generate_all_bonds, generate_residue_bonds




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

    plt.show()

    # Optional CSV
    if output_csv:
        df_results.to_csv(output_csv, index=False)
        print(f"Results written to {output_csv}")







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
    plt.savefig(out1); print(f"Saved PCS/SCS plot → {out1}")
    plt.show()

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
    plt.savefig(out2); print(f"Saved Element plot → {out2}")
    plt.show()












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
        print(f"[INFO] Link rendering: {missing} links skipped (atoms not present in current marker set).")

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
    print(f"Interactive plot saved as '{output_filename}'")
    fig.show()



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
    print(f"Interactive plot saved as '{output_filename}'")
    fig.show()







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
    print("Interactive plot saved as '1_conserved_interactive_modes_with_bonds.html'")
    fig.show()



