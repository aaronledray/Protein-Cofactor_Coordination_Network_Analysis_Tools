"""Tests for geometry-aware network alignment and position mapping."""

import unittest

import numpy as np

from modules.network_alignment import align_analysis_tables


def _atom(residue, number, atom, element, shell, x, y, z, motif="unknown_motif"):
    return {
        "residue_name": residue,
        "residue_number": number,
        "chain": "A",
        "insertion_code": "",
        "atom_name": atom,
        "element": element,
        "shell": shell,
        "motif": motif,
        "x": x,
        "y": y,
        "z": z,
    }


class NetworkAlignmentTests(unittest.TestCase):
    def test_kabsch_alignment_maps_rotated_network_to_reference_positions(self):
        reference_atoms = [
            _atom("OEX", 1, "MN1", "MN", "Cofactor", 0, 0, 0, "metal_cluster"),
            _atom("OEX", 1, "MN2", "MN", "Cofactor", 1, 0, 0, "metal_cluster"),
            _atom("OEX", 1, "MN3", "MN", "Cofactor", 0, 1, 0, "metal_cluster"),
            _atom("HIS", 42, "ND1", "N", "PCS", 0, 0, 1, "imidazole"),
        ]
        rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        translation = np.asarray([8.0, -3.0, 4.0])
        query_atoms = []
        for atom in reference_atoms:
            coordinate = rotation @ np.asarray([atom["x"], atom["y"], atom["z"]]) + translation
            query_atoms.append({**atom, "x": coordinate[0], "y": coordinate[1], "z": coordinate[2]})

        alignment = align_analysis_tables(
            {"atoms": reference_atoms},
            {"atoms": query_atoms},
            match_cutoff_A=0.1,
        )

        self.assertEqual(alignment["metrics"]["status"], "aligned")
        self.assertEqual(alignment["metrics"]["method"], "cofactor_kabsch")
        self.assertEqual(alignment["metrics"]["anchor_count"], 3)
        self.assertAlmostEqual(alignment["metrics"]["anchor_rmsd_A"], 0.0)
        self.assertEqual(alignment["metrics"]["matched_atom_count"], 4)
        self.assertEqual(alignment["metrics"]["residue_mapping_count"], 2)
        self.assertEqual(
            alignment["query_residue_map"][("HIS", "A", "42", "")],
            alignment["reference_residue_map"][("HIS", "A", "42", "")],
        )

    def test_insufficient_anchors_returns_explicit_fallback_status(self):
        reference = {"atoms": [_atom("CU", 1, "CU", "CU", "Cofactor", 0, 0, 0)]}
        query = {"atoms": [_atom("CU", 9, "CU", "CU", "Cofactor", 5, 5, 5)]}

        alignment = align_analysis_tables(reference, query)

        self.assertEqual(alignment["metrics"]["status"], "insufficient_anchors")
        self.assertEqual(alignment["metrics"]["method"], "none")
        self.assertEqual(alignment["reference_residue_map"], {})
        self.assertEqual(alignment["query_residue_map"], {})


if __name__ == "__main__":
    unittest.main()
