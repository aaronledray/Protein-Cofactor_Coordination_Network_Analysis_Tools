"""Standalone interactive Plotly viewer for cofactor-to-target wires."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import pandas as pd
import plotly.graph_objects as go


EDGE_COLORS = {
    "cofactor_contact": "#333333",
    "polar_contact": "#2166ac",
    "water_mediated_candidate": "#1b9e77",
    "aromatic_redox_contact": "#7b3294",
    "sulfur_contact": "#d95f02",
    "through_space_contact": "#999999",
}


def _records(table: Any) -> List[Dict[str, Any]]:
    if table is None:
        return []
    if hasattr(table, "to_dict"):
        return [dict(row) for row in table.to_dict("records")]
    return [dict(row) for row in table]


def _json_list(value: Any) -> List[Any]:
    if isinstance(value, (list, tuple)):
        return list(value)
    try:
        decoded = json.loads(str(value))
        return decoded if isinstance(decoded, list) else []
    except (TypeError, ValueError, json.JSONDecodeError):
        return []


def _node_label(node: Mapping[str, Any]) -> str:
    return (
        f"{node.get('residue_name', '?')} {node.get('residue_number', '?')}"
        f" {node.get('chain', '?')}:{node.get('atom_name', '?')}"
    )


def plot_interactive_protein_wire_network(
    nodes: Any,
    edges: Any,
    paths: Any,
    *,
    output_filename: str = "wire_network.html",
    pdb_name: str = "structure",
    cofactor_resname: str = "cofactor",
    wire_mode: str = "generic",
) -> None:
    """Write a self-contained interactive protein-wire viewer.

    The base graph remains visible while the path selector can highlight one
    ranked cofactor-to-target route at a time.  Edge and node controls are
    intentionally independent so users can inspect a path without losing its
    surrounding relay network.
    """
    node_rows = _records(nodes)
    edge_rows = _records(edges)
    path_rows = _records(paths)
    node_by_id = {str(row.get("node_id")): row for row in node_rows}

    role_groups = {"source": [], "relay": [], "target": []}
    for row in node_rows:
        role = str(row.get("node_role", "relay"))
        if "source" in role:
            role_groups["source"].append(row)
        elif "target" in role:
            role_groups["target"].append(row)
        else:
            role_groups["relay"].append(row)

    def marker_trace(rows: Sequence[Mapping[str, Any]], role: str, color: str, size: int) -> go.Scatter3d:
        return go.Scatter3d(
            x=[row.get("x") for row in rows],
            y=[row.get("y") for row in rows],
            z=[row.get("z") for row in rows],
            mode="markers",
            marker=dict(size=size, color=color, opacity=0.95, line=dict(color="#222222", width=1)),
            text=[
                (
                    f"{role.title()}<br>"
                    f"{_node_label(row)}<br>"
                    f"Motif: {row.get('motif', '?')}<br>"
                    f"Element: {row.get('element', '?')}<br>"
                    f"Capabilities: {row.get('relay_classes', '') or 'none'}<br>"
                    f"Target: {row.get('target_labels', '')}"
                )
                for row in rows
            ],
            customdata=[_node_label(row) for row in rows],
            hoverinfo="text",
            name=f"{role.title()} nodes",
            visible=True,
            showlegend=False,
        )

    fig = go.Figure()
    fig.add_trace(marker_trace(role_groups["source"], "source", "#111111", 9))
    fig.add_trace(marker_trace(role_groups["relay"], "relay", "#2166ac", 6))
    fig.add_trace(marker_trace(role_groups["target"], "target", "#d73027", 10))

    edge_groups: Dict[str, List[Mapping[str, Any]]] = {}
    for row in edge_rows:
        edge_groups.setdefault(str(row.get("interaction_type", "through_space_contact")), []).append(row)

    edge_trace_indices: Dict[str, int] = {}
    for edge_type in sorted(edge_groups):
        rows = edge_groups[edge_type]
        xs: List[Any] = []
        ys: List[Any] = []
        zs: List[Any] = []
        hover: List[Any] = []
        for row in rows:
            first = node_by_id.get(str(row.get("src_node_id")))
            second = node_by_id.get(str(row.get("dst_node_id")))
            if first is None or second is None:
                continue
            label = (
                f"{edge_type}<br>"
                f"{_node_label(first)} → {_node_label(second)}<br>"
                f"Geometry: {row.get('geometry_status', '?')}<br>"
                f"Support P/E/PCET: {float(row.get('proton_support', 0.0)):.2f} / "
                f"{float(row.get('electron_support', 0.0)):.2f} / "
                f"{float(row.get('pcet_support', 0.0)):.2f}<br>"
                f"Distance: {float(row.get('distance_A', 0.0)):.3f} Å"
            )
            xs.extend([first.get("x"), second.get("x"), None])
            ys.extend([first.get("y"), second.get("y"), None])
            zs.extend([first.get("z"), second.get("z"), None])
            hover.extend([label, label, None])
        edge_trace_indices[edge_type] = len(fig.data)
        fig.add_trace(go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="lines",
            line=dict(color=EDGE_COLORS.get(edge_type, "#777777"), width=4),
            text=hover,
            hoverinfo="text",
            name=edge_type.replace("_", " ").title(),
            visible=True,
            showlegend=False,
        ))

    path_trace_indices: List[int] = []
    path_labels: List[str] = []
    path_colors = ["#e41a1c", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf"]
    for path_index, path in enumerate(path_rows):
        node_ids = _json_list(path.get("node_ids"))
        edge_types = _json_list(path.get("edge_types"))
        xs: List[Any] = []
        ys: List[Any] = []
        zs: List[Any] = []
        hover: List[Any] = []
        for index in range(len(node_ids) - 1):
            first = node_by_id.get(str(node_ids[index]))
            second = node_by_id.get(str(node_ids[index + 1]))
            if first is None or second is None:
                continue
            edge_type = edge_types[index] if index < len(edge_types) else "wire hop"
            label = (
                f"{path.get('path_id', f'path_{path_index + 1}')}<br>"
                f"{_node_label(first)} → {_node_label(second)}<br>"
                f"Hop: {edge_type}"
            )
            xs.extend([first.get("x"), second.get("x"), None])
            ys.extend([first.get("y"), second.get("y"), None])
            zs.extend([first.get("z"), second.get("z"), None])
            hover.extend([label, label, None])
        path_trace_indices.append(len(fig.data))
        path_label = f"{path.get('path_id', f'path_{path_index + 1}')}: {path.get('target_label', 'target')}"
        path_labels.append(path_label)
        fig.add_trace(go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="lines",
            line=dict(color=path_colors[path_index % len(path_colors)], width=8),
            text=hover,
            hoverinfo="text",
            name=path_label,
            visible=False,
            showlegend=False,
        ))

    node_indices = [0, 1, 2]
    all_edge_indices = list(edge_trace_indices.values())
    all_trace_indices = list(range(len(fig.data)))
    base_visibility = [True] * len(fig.data)
    for index in path_trace_indices:
        base_visibility[index] = False

    path_buttons = [
        dict(
            label="Network edges",
            method="update",
            args=[{"visible": base_visibility}, {"annotations[5].text": "Wire view: network edges"}],
        )
    ]
    for path_index, path_label in enumerate(path_labels):
        selected_visibility = list(base_visibility)
        selected_visibility[path_trace_indices[path_index]] = True
        path_buttons.append(
            dict(
                label=path_label,
                method="update",
                args=[{"visible": selected_visibility}, {"annotations[5].text": f"Wire view: {path_label}"}],
            )
        )

    edge_type_buttons = [
        dict(label="All edge types", method="restyle", args=[{"visible": [True] * len(all_edge_indices)}, all_edge_indices]),
    ]
    for edge_type, trace_index in sorted(edge_trace_indices.items()):
        edge_type_buttons.append(
            dict(
                label=edge_type.replace("_", " ").title(),
                method="restyle",
                args=[{"visible": [index == trace_index for index in all_edge_indices]}, all_edge_indices],
            )
        )

    fig.update_layout(
        title={"text": f"Protein wire network: {cofactor_resname} → targets in {pdb_name} ({wire_mode})", "x": 0.60, "font": {"size": 22}},
        scene=dict(
            xaxis=dict(title="X (Å)", showgrid=True, gridcolor="#d9e1ec", showbackground=True, backgroundcolor="#e5ecf6"),
            yaxis=dict(title="Y (Å)", showgrid=True, gridcolor="#d9e1ec", showbackground=True, backgroundcolor="#e5ecf6"),
            zaxis=dict(title="Z (Å)", showgrid=True, gridcolor="#d9e1ec", showbackground=True, backgroundcolor="#e5ecf6"),
            bgcolor="#e5ecf6",
            aspectmode="data",
            domain=dict(x=[0.20, 1.0], y=[0.0, 1.0]),
        ),
        margin=dict(l=0, r=0, t=70, b=0),
        updatemenus=[
            dict(type="dropdown", direction="down", x=0.01, y=0.95, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, buttons=path_buttons),
            dict(type="dropdown", direction="down", x=0.01, y=0.84, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, buttons=edge_type_buttons),
            dict(type="dropdown", direction="down", x=0.01, y=0.73, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, buttons=[
                     dict(label="All nodes", method="restyle", args=[{"visible": [True, True, True]}, node_indices]),
                     dict(label="Source nodes", method="restyle", args=[{"visible": [True, False, False]}, node_indices]),
                     dict(label="Relay nodes", method="restyle", args=[{"visible": [False, True, False]}, node_indices]),
                     dict(label="Target nodes", method="restyle", args=[{"visible": [False, False, True]}, node_indices]),
                 ]),
            dict(type="dropdown", direction="down", x=0.01, y=0.62, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, buttons=[
                     dict(label="Labels on", method="restyle", args=[{"mode": ["markers+text", "markers+text", "markers+text"]}, node_indices]),
                     dict(label="Labels off", method="restyle", args=[{"mode": ["markers", "markers", "markers"]}, node_indices]),
                 ]),
            dict(type="dropdown", direction="down", x=0.01, y=0.51, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, buttons=[
                     dict(label="Grid on", method="relayout", args=[{"scene.xaxis.showgrid": True, "scene.yaxis.showgrid": True, "scene.zaxis.showgrid": True}]),
                     dict(label="Grid off", method="relayout", args=[{"scene.xaxis.showgrid": False, "scene.yaxis.showgrid": False, "scene.zaxis.showgrid": False}]),
                 ]),
            dict(type="dropdown", direction="up", x=0.01, y=0.40, xanchor="left", yanchor="top", showactive=True,
                 bgcolor="white", bordercolor="#b8c4d6", borderwidth=1, buttons=[
                     dict(label="Background on", method="relayout", args=[{"scene.bgcolor": "#e5ecf6", "scene.xaxis.showbackground": True, "scene.yaxis.showbackground": True, "scene.zaxis.showbackground": True}]),
                     dict(label="Background off", method="relayout", args=[{"scene.bgcolor": "rgba(0,0,0,0)", "scene.xaxis.showbackground": False, "scene.yaxis.showbackground": False, "scene.zaxis.showbackground": False}]),
                 ]),
        ],
        annotations=[
            dict(text="Wire paths", x=0.01, y=0.995, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Edge types", x=0.01, y=0.885, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Nodes", x=0.01, y=0.775, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Labels", x=0.01, y=0.665, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Grid", x=0.01, y=0.555, xref="paper", yref="paper", showarrow=False, font=dict(size=11, color="#52627a"), xanchor="left"),
            dict(text="Wire view: network edges", x=0.21, y=0.995, xref="paper", yref="paper", showarrow=False, bgcolor="rgba(255,255,255,0.8)", bordercolor="#b8c4d6", borderwidth=1, font=dict(size=10, color="#52627a"), xanchor="left"),
        ],
    )
    fig.write_html(
        output_filename,
        include_plotlyjs=True,
        config={
            "displayModeBar": True,
            "modeBarButtonsToAdd": ["orbitRotation", "tableRotation", "pan3d", "zoom3d", "resetCameraDefault3d", "toImage"],
            "displaylogo": False,
            "responsive": True,
        },
        div_id="protein-wire-plot",
    )


__all__ = ["plot_interactive_protein_wire_network"]
