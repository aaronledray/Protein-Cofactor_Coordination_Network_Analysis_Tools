"""Tests for the cohesive interactive single-structure viewer."""

from pathlib import Path
from tempfile import TemporaryDirectory
import re
import unittest
import json
import csv

from modules.coordination_api import analyze_structure
from modules.io_utils import unpack_pdb_file
from modules.moieties import bond_lookup
from modules.plotting import plot_interactive_cohesive_network


def read_html(path):
    """Read viewer HTML with JSON \\uXXXX escapes decoded.

    Newer plotly versions escape non-ASCII characters (Å, –, ·) in the embedded
    figure JSON, older ones write them literally; the viewer is the same.
    """
    text = Path(path).read_text(encoding="utf-8")
    return re.sub(r"\\u([0-9a-fA-F]{4})", lambda m: chr(int(m.group(1), 16)), text)


ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"
MYOGLOBIN = ROOT / "reference_structures/3_Myoglobin/1a6m.pdb"


class InteractiveViewerTests(unittest.TestCase):
    def _render_external_graph(self, csv_path, output_path):
        tables = analyze_structure(
            MYOGLOBIN,
            "HEM",
            exclude_moieties=["alanine_sidechain"],
            shells=3,
        )
        structure, _ = unpack_pdb_file(str(MYOGLOBIN))
        shell_groups = {}
        for row in tables["atoms"].to_dict("records"):
            shell_groups.setdefault(row["shell"], []).append(row)
        plot_interactive_cohesive_network(
            structure=structure,
            cofactor_atoms=shell_groups["Cofactor"],
            pcs_atoms=shell_groups["PCS"],
            scs_atoms=shell_groups["SCS"],
            focused_atoms_by_shell=shell_groups,
            contacts=tables["contacts"],
            bond_lookup_table=bond_lookup,
            pdb_name=MYOGLOBIN.name,
            cofactor_resname="HEM",
            output_filename=str(output_path),
            external_residue_network_csv=str(csv_path),
        )
        html_text = read_html(output_path)
        marker = 'Plotly.newPlot(                        "coordination-network-plot",'
        data_start = html_text.index(marker) + len(marker)
        payload = html_text[data_start:].lstrip()
        traces, end = json.JSONDecoder().raw_decode(payload)
        layout, _ = json.JSONDecoder().raw_decode(payload[end:].lstrip().lstrip(",").lstrip())
        return html_text, traces, layout

    def test_external_residue_graph_maps_cas_and_reports_unmapped_edges(self):
        with TemporaryDirectory(prefix="sscna-external-graph-small-") as work_dir:
            csv_path = Path(work_dir) / "small.csv"
            output = Path(work_dir) / "small.html"
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "source_chain", "source_residue_number", "target_chain",
                        "target_residue_number", "source_insertion_code",
                        "target_insertion_code", "edge_group", "color", "score",
                    ],
                )
                writer.writeheader()
                writer.writerow({
                    "source_chain": "A", "source_residue_number": "10",
                    "target_chain": "A", "target_residue_number": "11",
                    "edge_group": "example", "color": "#123456", "score": "0.91",
                })
                writer.writerow({
                    "source_chain": "A", "source_residue_number": "10",
                    "target_chain": "Z", "target_residue_number": "999",
                    "edge_group": "example", "score": "0.12",
                })

            html_text, traces, layout = self._render_external_graph(csv_path, output)
            external = [trace for trace in traces if trace["name"].startswith("External residue graph")]
            self.assertEqual(len(external), 1)
            self.assertEqual(sum(sum(value is None for value in trace["x"]) for trace in external), 1)
            structure, _ = unpack_pdb_file(str(MYOGLOBIN))
            source_ca = structure[0]["A"][10]["CA"].coord
            target_ca = structure[0]["A"][11]["CA"].coord
            self.assertEqual(external[0]["line"]["color"], "#123456")
            self.assertEqual(external[0]["x"][:2], [float(source_ca[0]), float(target_ca[0])])
            self.assertEqual(external[0]["y"][:2], [float(source_ca[1]), float(target_ca[1])])
            self.assertEqual(external[0]["z"][:2], [float(source_ca[2]), float(target_ca[2])])
            rendered = "\n".join(str(value) for trace in external for value in trace.get("text", []))
            self.assertIn("External residue-level graph", rendered)
            self.assertIn("not a CNI coordination contact", rendered)
            self.assertIn("score: 0.91", rendered)
            self.assertIn("External residue graph", html_text)
            self.assertIn("External graph (residue level)", html_text)
            self.assertIn("Mapped edges: 1/2", html_text)
            self.assertIn("Unmapped edges: 1", html_text)
            self.assertIn("A:10 → Z:999", html_text)
            preset = next(
                button for button in layout["updatemenus"][0]["buttons"]
                if button["label"] == "External residue graph"
            )
            self.assertTrue(preset["args"][0]["visible"][17])
            external_control = next(
                menu for menu in layout["updatemenus"]
                if [button["label"] for button in menu["buttons"]] == ["On", "Off"]
                and menu["buttons"][0]["args"][1] == [17]
            )
            self.assertEqual(external_control["buttons"][0]["args"][0]["visible"], [True])
            self.assertEqual(external_control["buttons"][1]["args"][0]["visible"], [False])
            self.assertIn("toImage", html_text)

    def test_external_residue_graph_renders_all_18_pairs_in_its_view_preset(self):
        structure, _ = unpack_pdb_file(str(MYOGLOBIN))
        ca_residues = [
            residue for residue in structure[0]["A"]
            if residue.child_dict.get("CA") is not None
        ][:19]
        self.assertEqual(len(ca_residues), 19)
        with TemporaryDirectory(prefix="sscna-external-graph-18-") as work_dir:
            csv_path = Path(work_dir) / "graph18.csv"
            output = Path(work_dir) / "graph18.html"
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "source_chain", "source_residue_number", "target_chain",
                        "target_residue_number", "edge_group", "evidence",
                    ],
                )
                writer.writeheader()
                expected_pairs = []
                for index in range(18):
                    source_number = ca_residues[index].id[1]
                    target_number = ca_residues[index + 1].id[1]
                    expected_pairs.append((source_number, target_number))
                    writer.writerow({
                        "source_chain": "A",
                        "source_residue_number": source_number,
                        "target_chain": "A",
                        "target_residue_number": target_number,
                        "edge_group": "sequence-balanced",
                        "evidence": f"manuscript edge {index + 1}",
                    })

            html_text, traces, layout = self._render_external_graph(csv_path, output)
            external = [trace for trace in traces if trace["name"].startswith("External residue graph")]
            self.assertEqual(len(external), 1)
            graph = external[0]
            self.assertEqual(sum(value is None for value in graph["x"]), 18)
            self.assertEqual(sum(value is None for value in graph["y"]), 18)
            self.assertEqual(sum(value is None for value in graph["z"]), 18)
            self.assertEqual(graph["name"], "External residue graph — sequence-balanced (18 edges)")
            hover_text = "\n".join(str(value) for value in graph["text"] if value)
            self.assertIn("evidence: manuscript edge 18", hover_text)
            edge_hovers = {str(value) for value in graph["text"] if value}
            self.assertTrue(all(
                any(
                    re.search(rf"Source: .* A:{source}<br>Target: .* A:{target}<br>", hover)
                    for hover in edge_hovers
                )
                for source, target in expected_pairs
            ))
            self.assertIn("External residue graph", html_text)
            self.assertIn("External residue graph (residue level)", html_text)
            self.assertIn("Mapped edges: 18/18", html_text)
            self.assertIn("Unmapped edges: none", html_text)
            preset = next(
                button for button in layout["updatemenus"][0]["buttons"]
                if button["label"] == "External residue graph"
            )
            self.assertTrue(preset["args"][0]["visible"][17])
            self.assertIn("toImage", html_text)

    def test_cohesive_viewer_embeds_all_network_representations(self):
        tables = analyze_structure(
            PLASTOCYANIN,
            "CU",
            exclude_moieties=["alanine_sidechain"],
            shells=3,
        )
        structure, _ = unpack_pdb_file(str(PLASTOCYANIN))
        shell_groups = {}
        for row in tables["atoms"].to_dict("records"):
            shell_groups.setdefault(row["shell"], []).append(row)

        with TemporaryDirectory(prefix="sscna-viewer-") as work_dir:
            output = Path(work_dir) / "cohesive.html"
            plot_interactive_cohesive_network(
                structure=structure,
                cofactor_atoms=shell_groups["Cofactor"],
                pcs_atoms=shell_groups["PCS"],
                scs_atoms=shell_groups["SCS"],
                focused_atoms_by_shell=shell_groups,
                contacts=tables["contacts"],
                bond_lookup_table=bond_lookup,
                pdb_name=PLASTOCYANIN.name,
                cofactor_resname="CU",
                output_filename=str(output),
            )

            html = read_html(output)
            self.assertIn("Full protein", html)
            self.assertIn("Coordination Network atoms", html)
            self.assertIn("Active-site context", html)
            self.assertIn("Actual coordinators", html)
            self.assertIn("primary_coordinator", html)
            self.assertIn("active_site_component", html)
            self.assertIn("Focused motif atoms", html)
            self.assertIn("Focused backbone atoms", html)
            self.assertIn("Network atoms + motifs", html)
            self.assertIn("Network atoms + backbone", html)
            self.assertIn("Network motifs only", html)
            self.assertIn("Network backbone only", html)
            self.assertIn("Focused motifs only", html)
            self.assertIn("Focused backbone only", html)
            self.assertIn("Backbone bonds", html)
            self.assertIn("Cofactor bonds", html)
            self.assertIn("Network residue bonds", html)
            self.assertIn("Residues", html)
            self.assertIn("Color by", html)
            self.assertIn("Coordination shells", html)
            self.assertIn("dropdown", html)
            self.assertIn("Atom layers", html)
            self.assertIn("Bonds", html)
            self.assertIn("Contacts", html)
            self.assertIn("Grid", html)
            self.assertIn("Background", html)
            self.assertIn("scene.bgcolor", html)
            self.assertIn("#e5ecf6", html)
            self.assertIn("rgba(0,0,0,0)", html)
            self.assertIn("showbackground", html)
            self.assertIn("Labels", html)
            self.assertIn("showticklabels", html)
            self.assertIn("scene.xaxis.title.text", html)
            self.assertIn("showgrid", html)
            self.assertIn("Display full protein structure", html)
            self.assertIn("View preset", html)
            self.assertIn("Network + motifs", html)
            self.assertIn("Structural context", html)
            self.assertIn("Full protein context", html)
            self.assertIn("View: Network", html)
            self.assertIn("X (Å)", html)
            self.assertIn("coordination-atom-details", html)
            self.assertIn("Click an atom to inspect its details.", html)
            self.assertIn("toImage", html)
            self.assertNotIn("coordinationOrbitTarget", html)
            self.assertNotIn("coordinationRightOrbit", html)
            self.assertIn("modeBarButtonsToAdd", html)
            self.assertIn("orbitRotation", html)
            self.assertIn("tableRotation", html)
            self.assertIn("Primary coordination", html)
            self.assertIn("Secondary coordination", html)
            self.assertIn("Tertiary coordination", html)
            self.assertIn("Tertiary", html)
            self.assertIn('"dash":"dot"', html)
            self.assertIn("imidazole", html)
            self.assertIn("ND1", html)

            # The third-shell trace must contain actual line segments, not
            # merely a label/control for an empty trace.
            marker = 'Plotly.newPlot(                        "coordination-network-plot",'
            data_start = html.index(marker) + len(marker)
            traces, _ = json.JSONDecoder().raw_decode(html[data_start:].lstrip())
            tertiary = next(trace for trace in traces if trace["name"] == "Tertiary coordination")
            self.assertGreaterEqual(len(tertiary["x"]), 3)
            self.assertEqual(tertiary["line"]["width"], 3)

    def test_heme_viewer_does_not_draw_proximity_pairs_as_coordination(self):
        tables = analyze_structure(
            MYOGLOBIN,
            "HEM",
            exclude_moieties=["alanine_sidechain"],
            shells=3,
        )
        structure, _ = unpack_pdb_file(str(MYOGLOBIN))
        shell_groups = {}
        for row in tables["atoms"].to_dict("records"):
            shell_groups.setdefault(row["shell"], []).append(row)

        with TemporaryDirectory(prefix="sscna-heme-viewer-") as work_dir:
            output = Path(work_dir) / "heme.html"
            plot_interactive_cohesive_network(
                structure=structure,
                cofactor_atoms=shell_groups["Cofactor"],
                pcs_atoms=shell_groups["PCS"],
                scs_atoms=shell_groups["SCS"],
                focused_atoms_by_shell=shell_groups,
                contacts=tables["contacts"],
                bond_lookup_table=bond_lookup,
                pdb_name=MYOGLOBIN.name,
                cofactor_resname="HEM",
                output_filename=str(output),
            )

            html = read_html(output)
            self.assertIn("HEM 154:FE", html)
            self.assertNotIn("HEM 154:C1A", html)
            self.assertIn("OXY 157:O1", html)
            self.assertIn("HIS 64:NE2", html)
            self.assertIn("heme_propionate", html)
            self.assertIn("Inferred O–H···O hydrogen bond", html)
            self.assertNotIn("HOH 1029:O [water] → HOH 1184:O [water]", html)


if __name__ == "__main__":
    unittest.main()
