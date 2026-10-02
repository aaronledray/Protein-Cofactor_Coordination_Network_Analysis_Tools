"""Tests for pairwise similarity matrix and heatmap artifacts."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from modules.comparison_plots import (
    build_residue_conservation_matrix,
    build_similarity_matrix,
    write_residue_conservation_heatmap,
    write_similarity_heatmap,
)


def _record(structure_id, edge_key):
    return {
        "structure_path": f"{structure_id}.pdb",
        "signature": {
            "structure_id": structure_id,
            "cofactor": {"key": "OEX"},
            "edges": {
                "primary": [{"match_key": edge_key}],
                "secondary": [],
                "tertiary": [],
            },
            "residues": [],
        },
    }


class ComparisonPlotTests(unittest.TestCase):
    def test_similarity_matrix_is_symmetric(self):
        labels, matrix, details = build_similarity_matrix([
            _record("a", "edge-a"),
            _record("b", "edge-a"),
            _record("c", "edge-c"),
        ])

        self.assertEqual(labels, ["a", "b", "c"])
        self.assertEqual(matrix.shape, (3, 3))
        self.assertTrue((matrix == matrix.T).all())
        self.assertEqual(matrix[0, 1], 1.0)
        self.assertLess(matrix[0, 2], 1.0)
        self.assertEqual(len(details), 3)

    def test_heatmap_writes_static_and_interactive_artifacts(self):
        with TemporaryDirectory(prefix="sscna-heatmap-") as work_dir:
            result = write_similarity_heatmap(
                [_record("a", "edge-a"), _record("b", "edge-b")],
                Path(work_dir),
            )

            self.assertTrue(Path(result["csv"]).is_file())
            self.assertTrue(Path(result["png"]).is_file())
            self.assertTrue(Path(result["html"]).is_file())
            self.assertIn("Coordination-network similarity", Path(result["html"]).read_text())

    def test_residue_conservation_map_tracks_aligned_positions(self):
        profile_result = {
            "template": {
                "structure_path": "template.pdb",
                "signature": {"structure_id": "template"},
                "position_order": [
                    "ref:SCS:GLU:A:64:",
                    "ref:PCS:HIS:A:37:",
                ],
            },
            "references": [
                {
                    "signature": {
                        "structure_id": "a",
                        "residues": [
                            {
                                "position": "ref:PCS:HIS:A:37:",
                                "match_key": "his-37",
                                "residue": "HIS",
                                "shell": "PCS",
                                "motif": "imidazole",
                            },
                        ],
                    },
                },
                {
                    "signature": {
                        "structure_id": "b",
                        "residues": [
                            {
                                "position": "ref:PCS:HIS:A:37:",
                                "match_key": "his-37",
                                "residue": "HIS",
                                "shell": "PCS",
                                "motif": "imidazole",
                            },
                            {
                                "position": "ref:SCS:GLU:A:64:",
                                "match_key": "glu-64",
                                "residue": "GLU",
                                "shell": "SCS",
                                "motif": "carboxylate",
                            },
                        ],
                    },
                },
            ],
            "profile": {
                "residues": {
                    "his-37": {"support": 1.0},
                    "glu-64": {"support": 0.5},
                },
            },
        }

        matrix = build_residue_conservation_matrix(profile_result)
        self.assertEqual(matrix["structure_labels"], ["a", "b"])
        self.assertEqual(
            matrix["positions"],
            ["ref:SCS:GLU:A:64:", "ref:PCS:HIS:A:37:"],
        )
        self.assertEqual(matrix["position_labels"], ["64", "37"])
        self.assertEqual(matrix["residue_types"], ["HIS", "GLU"])
        self.assertEqual(matrix["frequency_matrix"][0], [0.0, 0.5])
        self.assertEqual(matrix["frequency_matrix"][1], [1.0, 0.0])

        with TemporaryDirectory(prefix="sscna-conservation-map-") as work_dir:
            artifacts = write_residue_conservation_heatmap(profile_result, Path(work_dir))
            self.assertTrue(Path(artifacts["csv"]).is_file())
            self.assertTrue(Path(artifacts["png"]).is_file())
            self.assertTrue(Path(artifacts["html"]).is_file())
            self.assertIn("Family residue conservation", Path(artifacts["html"]).read_text())


if __name__ == "__main__":
    unittest.main()
