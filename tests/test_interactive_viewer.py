"""Tests for the cohesive interactive single-structure viewer."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import json

from modules.coordination_api import analyze_structure
from modules.io_utils import unpack_pdb_file
from modules.moieties import bond_lookup
from modules.plotting import plot_interactive_cohesive_network


ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"
MYOGLOBIN = ROOT / "reference_structures/3_Myoglobin/1a6m.pdb"


class InteractiveViewerTests(unittest.TestCase):
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

            html = output.read_text(encoding="utf-8")
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

            html = output.read_text(encoding="utf-8")
            self.assertIn("HEM 154:FE", html)
            self.assertNotIn("HEM 154:C1A", html)
            self.assertIn("OXY 157:O1", html)
            self.assertIn("HIS 64:NE2", html)
            self.assertIn("heme_propionate", html)
            self.assertIn("Inferred O–H···O hydrogen bond", html)
            self.assertNotIn("HOH 1029:O [water] → HOH 1184:O [water]", html)


if __name__ == "__main__":
    unittest.main()
