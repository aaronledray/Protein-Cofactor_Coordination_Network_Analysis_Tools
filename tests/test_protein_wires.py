"""Tests for the opt-in cofactor-to-target protein-wire graph."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from modules.protein_wires import _shortest_paths, build_protein_wire_network
from modules.wire_viewer import plot_interactive_protein_wire_network


def _atom(residue, number, chain, name, element, coordinates, model_id=0):
    return {
        "residue": residue,
        "residue_number": number,
        "chain": chain,
        "name": name,
        "element": element,
        "coordinates": coordinates,
        "model_id": model_id,
    }


class ProteinWireTests(unittest.TestCase):
    def test_builds_ranked_atom_level_cofactor_to_target_path(self):
        atoms = [
            _atom("HEM", 1500, "A", "FE", "FE", [0.0, 0.0, 0.0]),
            _atom("TRP", 107, "A", "NE1", "N", [2.5, 0.0, 0.0]),
            _atom("TRP", 107, "A", "CE2", "C", [2.7, 0.0, 0.0]),
            _atom("HIS", 108, "A", "NE2", "N", [5.5, 0.0, 0.0]),
            _atom("TYR", 200, "A", "OH", "O", [8.4, 0.0, 0.0]),
        ]

        result = build_protein_wire_network(
            atoms,
            "HEM",
            [{
                "label": "catalytic tyrosine",
                "residue": "TYR",
                "residue_number": 200,
                "chain": "A",
                "atom": "OH",
            }],
            max_hop_distance=3.6,
        )

        self.assertEqual(len(result["paths"]), 1)
        path = result["paths"].iloc[0]
        self.assertEqual(path["target_label"], "catalytic tyrosine")
        self.assertEqual(path["hops"], 3)
        node_ids = json.loads(path["node_ids"])
        self.assertEqual([node_id.split("|")[1] for node_id in node_ids], ["HEM", "TRP", "HIS", "TYR"])
        self.assertGreaterEqual(len(result["edges"]), 3)
        self.assertIn("cofactor_contact", set(result["edges"]["interaction_type"]))
        self.assertIn("aromatic_redox_contact", set(result["edges"]["interaction_type"]))

    def test_does_not_connect_atoms_within_one_residue(self):
        atoms = [
            _atom("HEM", 1500, "A", "FE", "FE", [0.0, 0.0, 0.0]),
            _atom("TRP", 107, "A", "NE1", "N", [2.5, 0.0, 0.0]),
            _atom("TRP", 107, "A", "CE2", "C", [2.6, 0.0, 0.0]),
            _atom("TYR", 200, "A", "OH", "O", [5.3, 0.0, 0.0]),
        ]

        result = build_protein_wire_network(
            atoms,
            "HEM",
            [{"residue": "TYR", "residue_number": 200, "chain": "A", "atom": "OH"}],
        )

        self.assertEqual(result["paths"].iloc[0]["hops"], 2)
        self.assertGreaterEqual(len(result["edges"]), 2)
        self.assertTrue(all(
            not (row.src_node_id.split("|")[1] == row.dst_node_id.split("|")[1] == "TRP")
            for row in result["edges"].itertuples()
        ))

    def test_water_is_optional_relay(self):
        atoms = [
            _atom("HEM", 1500, "A", "FE", "FE", [0.0, 0.0, 0.0]),
            _atom("HOH", 1, "A", "O", "O", [2.5, 0.0, 0.0]),
            _atom("TYR", 200, "A", "OH", "O", [5.3, 0.0, 0.0]),
        ]
        target = [{"residue": "TYR", "residue_number": 200, "chain": "A", "atom": "OH"}]

        with_water = build_protein_wire_network(atoms, "HEM", target, include_water=True)
        without_water = build_protein_wire_network(atoms, "HEM", target, include_water=False)

        self.assertEqual(with_water["paths"].iloc[0]["hops"], 2)
        self.assertTrue(without_water["paths"].empty)

    def test_explicit_hydrogen_geometry_is_retained(self):
        atoms = [
            _atom("HEM", 1500, "A", "FE", "FE", [-2.5, 0.0, 0.0]),
            _atom("SER", 10, "A", "OG", "O", [0.0, 0.0, 0.0]),
            _atom("SER", 10, "A", "HG", "H", [0.96, 0.0, 0.0]),
            _atom("TYR", 20, "A", "OH", "O", [1.92, 0.0, 0.0]),
        ]

        result = build_protein_wire_network(
            atoms,
            "HEM",
            [{"residue": "TYR", "residue_number": 20, "chain": "A", "atom": "OH"}],
            wire_mode="proton",
        )

        self.assertIn("validated_hydrogen_bond", set(result["edges"]["geometry_status"]))
        self.assertTrue((result["edges"]["wire_mode"] == "proton").all())
        self.assertIn("edge_geometry", result["paths"].columns)

    def test_electron_and_proton_modes_have_distinct_graph_policies(self):
        atoms = [
            _atom("HEM", 1500, "A", "FE", "FE", [0.0, 0.0, 0.0]),
            _atom("TRP", 107, "A", "NE1", "N", [2.5, 0.0, 0.0]),
            _atom("TRP", 107, "A", "CE2", "C", [2.7, 0.0, 0.0]),
            _atom("TYR", 200, "A", "OH", "O", [5.3, 0.0, 0.0]),
        ]
        target = [{"residue": "TYR", "residue_number": 200, "chain": "A", "atom": "OH"}]

        proton = build_protein_wire_network(atoms, "HEM", target, wire_mode="proton")
        electron = build_protein_wire_network(atoms, "HEM", target, wire_mode="electron")

        self.assertEqual(len(proton["nodes"]), len(electron["nodes"]))
        self.assertEqual(len(proton["edges"]), len(electron["edges"]))
        self.assertIn("redox_capable", electron["nodes"].columns)
        self.assertIn("protonatable", proton["nodes"].columns)
        self.assertNotEqual(
            float(proton["paths"].iloc[0]["path_cost"]),
            float(electron["paths"].iloc[0]["path_cost"]),
        )

    def test_wire_viewer_contains_path_and_geometry_controls(self):
        atoms = [
            _atom("HEM", 1500, "A", "FE", "FE", [0.0, 0.0, 0.0]),
            _atom("TRP", 107, "A", "NE1", "N", [2.5, 0.0, 0.0]),
            _atom("TYR", 200, "A", "OH", "O", [5.3, 0.0, 0.0]),
        ]
        result = build_protein_wire_network(
            atoms,
            "HEM",
            [{"label": "target", "residue": "TYR", "residue_number": 200, "chain": "A", "atom": "OH"}],
        )
        with TemporaryDirectory(prefix="sscna-wire-viewer-") as work_dir:
            output = Path(work_dir) / "wire.html"
            plot_interactive_protein_wire_network(
                result["nodes"], result["edges"], result["paths"],
                output_filename=str(output),
                pdb_name="wire_test.pdb",
                cofactor_resname="HEM",
                wire_mode="generic",
            )
            html = output.read_text(encoding="utf-8")
            self.assertIn("Protein wire network", html)
            self.assertIn("Wire paths", html)
            self.assertIn("Edge types", html)
            self.assertIn("Labels on", html)
            self.assertIn("Background off", html)
            self.assertIn("Geometry:", html)

    def test_redox_mode_allows_wider_aromatic_hops_and_keeps_residue_paths_unique(self):
        atoms = [
            _atom("HEM", 1500, "A", "FE", "FE", [0.0, 0.0, 0.0]),
            _atom("TRP", 90, "A", "NE1", "N", [4.5, 0.0, 0.0]),
            _atom("TYR", 200, "A", "OH", "O", [8.5, 0.0, 0.0]),
        ]
        result = build_protein_wire_network(
            atoms,
            "HEM",
            [{"residue": "TYR", "residue_number": 200, "chain": "A", "atom": "OH"}],
            max_hop_distance=3.6,
            redox_hop_distance=5.0,
            wire_mode="redox",
        )

        self.assertEqual(len(result["paths"]), 1)
        residue_path = json.loads(result["paths"].iloc[0]["residue_path"])
        residue_keys = [
            (row["residue"], str(row["residue_number"]), row["chain"])
            for row in residue_path
        ]
        self.assertEqual(len(residue_keys), len(set(residue_keys)))
        self.assertIn("redox_capable", result["nodes"].columns)

    def test_target_status_reports_numbering_mismatch(self):
        atoms = [
            _atom("HEM", 1500, "A", "FE", "FE", [0.0, 0.0, 0.0]),
            _atom("TRP", 91, "A", "NE1", "N", [2.5, 0.0, 0.0]),
            _atom("ALA", 93, "A", "CA", "C", [12.0, 0.0, 0.0]),
        ]
        result = build_protein_wire_network(
            atoms,
            "HEM",
            [
                {"label": "mapped tryptophan", "residue": "TRP", "residue_number": 91, "chain": "A", "atom": "NE1"},
                {"label": "wrong numbering", "residue": "TRP", "residue_number": 93, "chain": "A", "atom": "NE1"},
            ],
        )

        status = result["target_status"].set_index("target_label")
        self.assertEqual(status.loc["mapped tryptophan", "status"], "matched")
        self.assertFalse(status.loc["wrong numbering", "matched"])
        self.assertEqual(status.loc["wrong numbering", "status"], "residue_name_mismatch")
        self.assertEqual(status.loc["wrong numbering", "observed_residues"], "ALA")

    def test_unique_residue_constraint_rejects_atom_level_residue_revisit(self):
        def edge(src, dst):
            return {
                "interaction_type": "aromatic_redox_contact",
                "geometry_status": "aromatic_face_to_face",
                "edge_cost": 1.0,
            }

        adjacency = {
            "source": [("his_first", edge("source", "his_first"))],
            "his_first": [("tyr", edge("his_first", "tyr"))],
            "tyr": [("his_second", edge("tyr", "his_second"))],
            "his_second": [("target", edge("his_second", "target"))],
        }
        residue_keys = {
            "source": ("", "HEM", "1500", "A", "", ""),
            "his_first": ("", "HIS", "270", "A", "", ""),
            "tyr": ("", "TYR", "229", "A", "", ""),
            "his_second": ("", "HIS", "270", "A", "", ""),
            "target": ("", "TRP", "321", "A", "", ""),
        }
        paths = _shortest_paths(
            adjacency,
            {"source"},
            {"target"},
            {"target": ("TRP:A:321:NE1",)},
            max_hops=8,
            wire_mode="redox",
            residue_key_by_node=residue_keys,
            enforce_unique_residues=True,
        )
        self.assertEqual(paths, [])


if __name__ == "__main__":
    unittest.main()
