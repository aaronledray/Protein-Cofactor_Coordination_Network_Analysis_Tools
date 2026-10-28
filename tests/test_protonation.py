"""Tests for hydrogen-bond geometry helpers; PROPKA and reduce are exercised only when installed."""

import importlib.util
import shutil
import unittest

import numpy as np
import pandas as pd

from modules import protonation as pr


class GeometryTests(unittest.TestCase):
    def test_directional_hbond_needs_angle_and_distance(self):
        donor, hydrogen = np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])
        self.assertTrue(pr.directional_hbonds([(donor, hydrogen)], np.array([[2.8, 0.0, 0.0]])))
        self.assertFalse(pr.directional_hbonds([(donor, hydrogen)], np.array([[0.5, 2.5, 0.0]])))   # hydrogen points away
        self.assertFalse(pr.directional_hbonds([(donor, hydrogen)], np.array([[4.5, 0.0, 0.0]])))   # too far

    def test_virtual_hydrogens_only_on_unprotonated_ring_nitrogens(self):
        theta = np.linspace(0, 2 * np.pi, 6)[:-1]
        ring = {n: np.array([1.2 * np.cos(t), 1.2 * np.sin(t), 0.0]) for n, t in zip(pr.HIS_RING, theta)}
        both = pr.virtual_histidine_hydrogens(ring, [])
        self.assertEqual(len(both), 2)
        for n, h in both:
            self.assertAlmostEqual(float(np.linalg.norm(h - n)), 1.0, places=6)
        ne2_h = ring["NE2"] + np.array([0.0, 0.0, 1.0])
        self.assertEqual(len(pr.virtual_histidine_hydrogens(ring, [ne2_h])), 1)

    def test_carboxylate_distances(self):
        atoms = pd.DataFrame([("s", 399, "ASP", "OD1", 0.0, 0.0, 0.0), ("s", 399, "ASP", "OD2", 1.0, 0.0, 0.0),
                              ("s", 5, "ARG", "NH1", 0.0, 0.0, 0.0)], columns=["structure_id", "position", "resname", "atom", "x", "y", "z"])
        acceptors = pd.DataFrame([("s", 1.0, 2.6, 0.0)], columns=["structure_id", "x", "y", "z"])
        d = pr.carboxylate_distances(atoms, acceptors)
        self.assertEqual(list(d["position"]), [399])
        self.assertAlmostEqual(d["oo_distance"].iloc[0], 2.6, places=6)


@unittest.skipUnless(importlib.util.find_spec("propka"), "PROPKA not installed")
class PropkaTests(unittest.TestCase):
    def test_reference_structure_has_titratable_groups(self):
        from pathlib import Path

        cif = Path(__file__).resolve().parents[1] / "reference_structures/0_Plastocyanin/1ag6.cif"
        from Bio.PDB import MMCIFParser, PDBIO
        import tempfile

        structure = MMCIFParser(QUIET=True).get_structure("x", str(cif))
        with tempfile.TemporaryDirectory() as td:
            pdb = Path(td) / "x.pdb"
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(pdb))
            table = pr.run_propka(str(pdb))
        self.assertGreater(len(table), 5)
        self.assertIn("HIS", set(table["residue_name"]))


@unittest.skipUnless(shutil.which("reduce"), "reduce not installed")
class ReduceTests(unittest.TestCase):
    def test_hydrogens_are_added(self):
        from pathlib import Path
        import tempfile

        from Bio.PDB import MMCIFParser, PDBIO

        cif = Path(__file__).resolve().parents[1] / "reference_structures/0_Plastocyanin/1ag6.cif"
        structure = MMCIFParser(QUIET=True).get_structure("x", str(cif))
        with tempfile.TemporaryDirectory() as td:
            pdb = Path(td) / "x.pdb"
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(pdb))
            out = pr.add_hydrogens(str(pdb), str(Path(td) / "xH.pdb"))
            with open(out) as handle:
                n_h = sum(1 for line in handle if line.startswith(("ATOM", "HETATM")) and line[76:78].strip() == "H")
        self.assertGreater(n_h, 100)


if __name__ == "__main__":
    unittest.main()
