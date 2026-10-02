"""Static and interactive plots for coordination-network comparisons."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np

from .network_comparison import compare_network_signatures


def build_similarity_matrix(
    records: Sequence[Mapping[str, Any]],
    *,
    distance_tolerance_A: float = 0.5,
) -> Tuple[List[str], np.ndarray, List[Dict[str, Any]]]:
    """Return labels, a symmetric overall-score matrix, and pairwise details."""
    if not records:
        raise ValueError("At least one signature record is required")
    labels: List[str] = []
    seen: Dict[str, int] = {}
    for record in records:
        signature = record["signature"]
        base = str(signature.get("structure_id") or "structure")
        seen[base] = seen.get(base, 0) + 1
        labels.append(base if seen[base] == 1 else f"{base}__{seen[base]}")

    matrix = np.eye(len(records), dtype=float)
    details: List[Dict[str, Any]] = []
    for left_index in range(len(records)):
        for right_index in range(left_index + 1, len(records)):
            comparison = compare_network_signatures(
                records[left_index]["signature"],
                records[right_index]["signature"],
                distance_tolerance_A=distance_tolerance_A,
            )
            score = float(comparison["scores"]["overall"])
            matrix[left_index, right_index] = score
            matrix[right_index, left_index] = score
            details.append({
                "left_id": labels[left_index],
                "right_id": labels[right_index],
                "left_path": records[left_index].get("structure_path", ""),
                "right_path": records[right_index].get("structure_path", ""),
                **comparison,
            })
    return labels, matrix, details


def write_similarity_heatmap(
    records: Sequence[Mapping[str, Any]],
    output_dir: Path,
    *,
    stem: str = "network_similarity_heatmap",
    title: str = "Coordination-network similarity",
    distance_tolerance_A: float = 0.5,
) -> Dict[str, Any]:
    """Write CSV, PNG, and interactive HTML similarity heatmaps."""
    labels, matrix, details = build_similarity_matrix(
        records,
        distance_tolerance_A=distance_tolerance_A,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_csv = output_dir / f"{stem}.csv"
    with matrix_csv.open("w", encoding="utf-8") as handle:
        handle.write(",".join(["structure_id", *labels]) + "\n")
        for label, row in zip(labels, matrix):
            values = ",".join(f"{value:.6f}" for value in row)
            handle.write(f"{label},{values}\n")

    # Matplotlib gives us a dependable static artifact for reports and GitHub
    # Pages, while Plotly provides hoverable pairwise values.
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    figure_size = max(6.0, min(14.0, 1.0 + 0.55 * len(labels)))
    figure, axis = plt.subplots(figsize=(figure_size, figure_size))
    image = axis.imshow(matrix, vmin=0.0, vmax=1.0, cmap="viridis")
    axis.set_xticks(range(len(labels)), labels, rotation=45, ha="right")
    axis.set_yticks(range(len(labels)), labels)
    axis.set_title(title)
    axis.set_xlabel("Structure")
    axis.set_ylabel("Structure")
    colorbar = figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    colorbar.set_label("Overall similarity")
    if len(labels) <= 20:
        for row_index in range(len(labels)):
            for column_index in range(len(labels)):
                text_color = "white" if matrix[row_index, column_index] < 0.55 else "black"
                axis.text(
                    column_index,
                    row_index,
                    f"{matrix[row_index, column_index]:.2f}",
                    ha="center",
                    va="center",
                    color=text_color,
                    fontsize=8,
                )
    figure.tight_layout()
    png_path = output_dir / f"{stem}.png"
    figure.savefig(png_path, dpi=180)
    plt.close(figure)

    import plotly.graph_objects as go

    text = [[f"{value:.3f}" for value in row] for row in matrix]
    plot = go.Figure(
        data=go.Heatmap(
            z=matrix.tolist(),
            x=labels,
            y=labels,
            zmin=0.0,
            zmax=1.0,
            colorscale="Viridis",
            text=text,
            texttemplate="%{text}",
            hovertemplate="%{y} vs %{x}<br>Similarity: %{z:.3f}<extra></extra>",
            colorbar={"title": "Similarity"},
        )
    )
    plot.update_layout(
        title=title,
        xaxis_title="Structure",
        yaxis_title="Structure",
        width=max(700, min(1400, 180 + 65 * len(labels))),
        height=max(650, min(1400, 180 + 65 * len(labels))),
        margin={"l": 90, "r": 30, "t": 70, "b": 120},
    )
    html_path = output_dir / f"{stem}.html"
    plot.write_html(str(html_path), include_plotlyjs=True)
    return {
        "labels": labels,
        "matrix": matrix.tolist(),
        "details": details,
        "csv": str(matrix_csv),
        "png": str(png_path),
        "html": str(html_path),
    }


def _unique_structure_labels(records: Sequence[Mapping[str, Any]]) -> List[str]:
    labels: List[str] = []
    seen: Dict[str, int] = {}
    for record in records:
        signature = record.get("signature", {})
        base = str(signature.get("structure_id") or "structure")
        seen[base] = seen.get(base, 0) + 1
        labels.append(base if seen[base] == 1 else f"{base}__{seen[base]}")
    return labels


def _display_position(position: str) -> str:
    """Shorten a stable position key for plot axes without losing its identity."""
    if position.startswith("ref:"):
        fields = position[4:].split(":")
        if len(fields) >= 5:
            _, residue, chain, number, insertion = fields[:5]
            location = f"{chain}:{number}{insertion}"
            return f"{residue} {location}" if chain or number else residue
        return position[4:]
    return position


def _template_residue_number(position: str) -> str:
    """Return the template residue number used on the legacy-style y-axis."""
    if position.startswith("ref:"):
        fields = position[4:].split(":")
        if len(fields) >= 5:
            return f"{fields[3]}{fields[4]}"
    return position


def build_residue_conservation_matrix(
    profile_result: Mapping[str, Any],
) -> Dict[str, Any]:
    """Build a template-numbered residue-type frequency matrix.

    This follows the legacy heatmap convention: rows are template residue
    numbers, columns are observed residue types, and each cell is the fraction
    of reference structures carrying that residue type at that template
    position.  Full position keys and structure membership remain available in
    ``cells`` for unambiguous downstream use.
    """
    records = list(profile_result.get("references", []))
    if not records:
        raise ValueError("Profile result has no reference records")
    labels = _unique_structure_labels(records)
    template = profile_result.get("template", {})

    features_by_record: List[Dict[str, str]] = []
    positions = set()
    for record in records:
        by_position: Dict[str, str] = {}
        for feature in record.get("signature", {}).get("residues", []):
            position = str(feature.get("position") or "").strip()
            residue = str(feature.get("residue") or "").strip().upper()
            if not position:
                continue
            positions.add(position)
            if not residue:
                continue
            existing = by_position.get(position)
            if existing and residue not in existing.split("/"):
                by_position[position] = "/".join(sorted({existing, residue}))
            else:
                by_position[position] = residue
        features_by_record.append(by_position)
    if not positions:
        raise ValueError(
            "Reference signatures do not contain aligned residue positions; "
            "enable alignment to build a 2D conservation map"
        )

    template_positions = [
        str(position) for position in template.get("position_order", []) if str(position)
    ]
    # Use the template's encounter order and numbering first.  Any mapped
    # positions not present in the template metadata are appended explicitly
    # rather than silently dropped.
    position_labels = [position for position in template_positions if position in positions]
    position_labels.extend(sorted(positions - set(template_positions)))
    if not position_labels:
        position_labels = sorted(positions)
    counts_by_position: Dict[str, Dict[str, int]] = {}
    first_seen_types: List[str] = []
    for position in position_labels:
        counts: Dict[str, int] = {}
        for by_position in features_by_record:
            residue = by_position.get(position)
            if not residue:
                continue
            counts[residue] = counts.get(residue, 0) + 1
            if residue not in first_seen_types:
                first_seen_types.append(residue)
        counts_by_position[position] = counts
    residue_types = sorted(
        first_seen_types,
        key=lambda residue: (
            -sum(counts.get(residue, 0) for counts in counts_by_position.values()),
            first_seen_types.index(residue),
        ),
    )

    frequency_matrix: List[List[float]] = []
    cells: List[Dict[str, Any]] = []
    for position in position_labels:
        counts = counts_by_position[position]
        frequency_row = [round(counts.get(residue, 0) / len(records), 6) for residue in residue_types]
        frequency_matrix.append(frequency_row)
        for residue, frequency in zip(residue_types, frequency_row):
            matching_structures = [
                labels[index]
                for index, by_position in enumerate(features_by_record)
                if by_position.get(position) == residue
            ]
            cells.append({
                "position": position,
                "template_residue_number": _template_residue_number(position),
                "position_label": _display_position(position),
                "residue_type": residue,
                "count": len(matching_structures),
                "support": frequency,
                "structure_ids": matching_structures,
            })
    text_matrix = [[f"{value:.2f}" for value in row] for row in frequency_matrix]

    return {
        "structure_labels": labels,
        "template_structure_id": template.get("signature", {}).get("structure_id", ""),
        "template_structure_path": template.get("structure_path", ""),
        "positions": position_labels,
        "position_labels": [_template_residue_number(position) for position in position_labels],
        "residue_types": residue_types,
        "frequency_matrix": frequency_matrix,
        "support_matrix": frequency_matrix,
        "text_matrix": text_matrix,
        "cells": cells,
    }


def write_residue_conservation_heatmap(
    profile_result: Mapping[str, Any],
    output_dir: Path,
    *,
    stem: str = "residue_conservation_map",
    title: str = "Family residue conservation",
) -> Dict[str, Any]:
    """Write CSV, PNG, and interactive template-numbered frequency maps."""
    result = build_residue_conservation_matrix(profile_result)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{stem}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Template residue number", *result["residue_types"]])
        for position_label, row in zip(result["position_labels"], result["frequency_matrix"]):
            writer.writerow([position_label, *[f"{value:.6f}" for value in row]])

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    matrix = np.asarray(result["frequency_matrix"], dtype=float)
    figure_width = max(9.0, min(22.0, 3.0 + 0.55 * len(result["residue_types"])))
    figure_height = max(5.0, min(16.0, 2.0 + 0.30 * len(result["positions"])))
    figure, axis = plt.subplots(figsize=(figure_width, figure_height))
    image = axis.imshow(matrix, vmin=0.0, vmax=1.0, cmap="magma", aspect="auto")
    axis.set_xticks(range(len(result["residue_types"])), result["residue_types"], rotation=45, ha="right")
    axis.set_yticks(range(len(result["positions"])), result["position_labels"])
    axis.set_title(title)
    axis.set_xlabel("Residue type")
    axis.set_ylabel("Template residue number")
    colorbar = figure.colorbar(image, ax=axis, fraction=0.025, pad=0.02)
    colorbar.set_label("Frequency")
    for row_index, text_row in enumerate(result["text_matrix"]):
        for column_index, text in enumerate(text_row):
            value = matrix[row_index, column_index]
            text_color = "white" if value < 0.55 else "black"
            axis.text(
                column_index,
                row_index,
                text,
                ha="center",
                va="center",
                color=text_color,
                fontsize=7,
            )
    figure.tight_layout()
    png_path = output_dir / f"{stem}.png"
    figure.savefig(png_path, dpi=180)
    plt.close(figure)

    import plotly.graph_objects as go

    cells_by_key = {
        (cell["position"], cell["residue_type"]): cell
        for cell in result["cells"]
    }
    hover_text = []
    for row_index, position in enumerate(result["positions"]):
        row = []
        for residue_type in result["residue_types"]:
            cell = cells_by_key[(position, residue_type)]
            row.append(
                f"Template residue {cell['template_residue_number']}"
                f" ({cell['position_label']})<br>Residue: {residue_type}<br>"
                f"Frequency: {cell['support']:.1%} ({cell['count']}/{len(result['structure_labels'])})<br>"
                f"Structures: {', '.join(cell['structure_ids']) or 'none'}"
            )
        hover_text.append(row)
    plot = go.Figure(
        data=go.Heatmap(
            z=matrix.tolist(),
            x=result["residue_types"],
            y=result["position_labels"],
            zmin=0.0,
            zmax=1.0,
            colorscale="Magma",
            text=result["text_matrix"],
            texttemplate="%{text}",
            customdata=hover_text,
            hovertemplate="%{customdata}<extra></extra>",
            colorbar={"title": "Frequency"},
        )
    )
    plot.update_layout(
        title=title,
        xaxis_title="Residue type",
        yaxis_title="Template residue number",
        width=max(900, min(1800, 360 + 48 * len(result["residue_types"]))),
        height=max(520, min(1400, 180 + 28 * len(result["positions"]))),
        margin={"l": 110, "r": 90, "t": 70, "b": 130},
    )
    html_path = output_dir / f"{stem}.html"
    plot.write_html(str(html_path), include_plotlyjs=True)
    result.update({"csv": str(csv_path), "png": str(png_path), "html": str(html_path)})
    return result


__all__ = [
    "build_similarity_matrix",
    "write_similarity_heatmap",
    "build_residue_conservation_matrix",
    "write_residue_conservation_heatmap",
]
