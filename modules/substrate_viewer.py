"""Interactive 3D view of a substrate point and its chains through the network."""

from typing import Any, Mapping, Optional

import pandas as pd
import plotly.graph_objects as go

SHELL_COLORS = {"Cofactor": "#111111", "PCS": "#2166ac", "SCS": "#c51b7d", "TCS": "#1b7837"}
FALLBACK_COLOR = "#7f7f7f"
CHAIN_COLORS = ["#e66101", "#5e3c99", "#018571", "#a6611a", "#d01c8b", "#4dac26"]
VISIBLE_CHAINS_BY_DEFAULT = 3


def _label(row: Mapping[str, Any]) -> str:
    return f'{row["residue_name"]}{row["residue_number"]}:{row["atom_name"]}'


def build_substrate_figure(
    tables: Mapping[str, pd.DataFrame],
    pdb_name: str = "",
    site_id: Optional[str] = None,
) -> go.Figure:
    """Build a figure for one site: network atoms, the point, and ranked chains.

    ``tables`` is the result of ``analyze_substrate_seed``. The first few chains
    are visible initially; the legend toggles each chain and the contact lines.
    """
    atoms, seed = tables["atoms"], tables["seed"]
    if seed.empty:
        raise ValueError("no substrate seed to draw")
    site = site_id if site_id is not None else str(seed["site_id"].iloc[0])
    site_atoms = atoms[atoms["site_id"] == site]
    point = seed[seed["site_id"] == site].iloc[0]
    figure = go.Figure()

    for shell, group in site_atoms.groupby("shell", sort=False):
        figure.add_trace(go.Scatter3d(
            x=group["x"], y=group["y"], z=group["z"], mode="markers", name=str(shell),
            marker={"size": 4, "color": SHELL_COLORS.get(str(shell), FALLBACK_COLOR), "opacity": 0.7},
            text=[f"{_label(r)} ({shell}, {r['motif']})" for r in group.to_dict("records")],
            hoverinfo="text",
        ))

    contacts = tables["substrate_contacts"]
    contacts = contacts[contacts["site_id"] == site]
    if not contacts.empty:
        lookup = {(str(r["residue_name"]), str(r["residue_number"]), str(r["chain"]), str(r["atom_name"])): r
                  for r in site_atoms.to_dict("records")}
        xs, ys, zs = [], [], []
        for row in contacts.to_dict("records"):
            atom = lookup.get((str(row["residue_name"]), str(row["residue_number"]),
                               str(row["chain"]), str(row["atom_name"])))
            if atom is not None:
                xs += [point["x"], atom["x"], None]
                ys += [point["y"], atom["y"], None]
                zs += [point["z"], atom["z"], None]
        figure.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="lines", name="point contacts", visible="legendonly",
            line={"color": "#999999", "width": 2, "dash": "dot"}, hoverinfo="skip",
        ))

    chains = tables["chains"][tables["chains"]["site_id"] == site]
    steps = tables["steps"]
    for index, chain in enumerate(chains.to_dict("records")):
        path = steps[steps["chain_id"] == chain["chain_id"]].sort_values("step")
        coordinates = []
        for row in path.to_dict("records"):
            if row["node_type"] == "substrate_point":
                coordinates.append((point["x"], point["y"], point["z"], "substrate point"))
            else:
                match = site_atoms[
                    (site_atoms["residue_name"] == row["residue_name"])
                    & (site_atoms["residue_number"] == row["residue_number"])
                    & (site_atoms["chain"] == row["chain"])
                    & (site_atoms["atom_name"] == row["atom_name"])
                ].iloc[0]
                coordinates.append((match["x"], match["y"], match["z"], _label(row)))
        figure.add_trace(go.Scatter3d(
            x=[c[0] for c in coordinates], y=[c[1] for c in coordinates],
            z=[c[2] for c in coordinates], mode="lines+markers+text",
            name=f'#{chain["rank"]} {chain["entry_residue"]} ({chain["n_hops"]} hops, {chain["total_length_A"]} Å)',
            text=[c[3] for c in coordinates], textposition="top center",
            line={"color": CHAIN_COLORS[index % len(CHAIN_COLORS)], "width": 6},
            marker={"size": 5, "color": CHAIN_COLORS[index % len(CHAIN_COLORS)]},
            visible=True if index < VISIBLE_CHAINS_BY_DEFAULT else "legendonly",
            hovertext=[f'{c[3]}' for c in coordinates], hoverinfo="text",
        ))

    figure.add_trace(go.Scatter3d(
        x=[point["x"]], y=[point["y"]], z=[point["z"]], mode="markers", name="substrate point",
        marker={"size": 9, "color": "#ffd92f", "symbol": "diamond",
                "line": {"color": "#000000", "width": 2}},
        hovertext=[f'hypothetical substrate point ({point["x"]:.2f}, {point["y"]:.2f}, {point["z"]:.2f})'],
        hoverinfo="text",
    ))
    figure.update_layout(
        title=f"{pdb_name} {site}: substrate-point chains (structural hypothesis)",
        scene={"xaxis_title": "x (Å)", "yaxis_title": "y (Å)", "zaxis_title": "z (Å)", "aspectmode": "data"},
        legend={"itemsizing": "constant"}, margin={"l": 0, "r": 0, "b": 0, "t": 40},
    )
    return figure


def plot_interactive_substrate_chains(
    tables: Mapping[str, pd.DataFrame],
    output_filename: str,
    pdb_name: str = "",
    site_id: Optional[str] = None,
    compact_html: bool = False,
) -> None:
    figure = build_substrate_figure(tables, pdb_name=pdb_name, site_id=site_id)
    figure.write_html(output_filename, include_plotlyjs="cdn" if compact_html else True)
