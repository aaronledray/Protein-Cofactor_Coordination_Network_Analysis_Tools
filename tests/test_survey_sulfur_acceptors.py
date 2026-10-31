"""Offline checks for scripts/survey_sulfur_acceptors.py: hydrogen placement and the 1AG6 thiolate."""

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CIF = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("survey_sulfur_acceptors", ROOT / "scripts/survey_sulfur_acceptors.py")
survey = importlib.util.module_from_spec(spec)
spec.loader.exec_module(survey)


class SulfurSurveyTests(unittest.TestCase):
    def test_amide_hydrogen_is_on_the_bisector_at_nh_length(self):
        n, c_prev, ca = np.zeros(3), np.array([-1.33, 0.0, 0.0]), np.array([0.5, 1.3, 0.0])
        h = survey.amide_hydrogen(n, c_prev, ca)
        self.assertAlmostEqual(float(np.linalg.norm(h - n)), survey.NH_LENGTH, places=6)
        to_c, to_ca = (c_prev - n) / np.linalg.norm(c_prev - n), (ca - n) / np.linalg.norm(ca - n)
        self.assertAlmostEqual(float(np.dot(h - n, to_c)), float(np.dot(h - n, to_ca)), places=6)      # equidistant in angle
        self.assertLess(float(np.dot(h - n, to_c + to_ca)), 0.0)                                      # points away from both neighbors

    def test_geometry_of_a_linear_hydrogen_bond(self):
        n, s = np.zeros(3), np.array([3.4, 0.0, 0.0])
        h = np.array([1.01, 0.0, 0.0])
        d, angle = survey.nh_geometry(h, n, s)
        self.assertAlmostEqual(d, 2.39, places=2)
        self.assertAlmostEqual(angle, 180.0, places=3)

    @staticmethod
    def stub(resname, coords):
        class Atom:
            def __init__(self, xyz):
                self.coord = np.array(xyz, float)

        class Residue(dict):
            pass

        residue = Residue({name: Atom(xyz) for name, xyz in coords.items()})
        residue.resname = resname
        return residue

    @staticmethod
    def angle(a, b, c):
        u, v = a - b, c - b
        return float(np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))))

    def test_ring_nh_hydrogen_is_on_the_exocyclic_bisector(self):
        trp = self.stub("TRP", {"NE1": (0, 0, 0), "CD1": (-1.1, 1.1, 0), "CE2": (1.1, 1.1, 0)})
        (h,) = survey.planar_hydrogens(trp, "NE1")
        self.assertAlmostEqual(float(np.linalg.norm(h)), survey.NH_LENGTH, places=6)
        self.assertAlmostEqual(self.angle(trp["CD1"].coord, trp["NE1"].coord, h), self.angle(trp["CE2"].coord, trp["NE1"].coord, h), places=6)
        self.assertAlmostEqual(float(h[2]), 0.0, places=6)                  # in the ring plane
        self.assertLess(float(h[1]), 0.0)                                   # away from the ring

    def test_amide_hydrogens_are_planar_at_120_degrees(self):
        gln = self.stub("GLN", {"NE2": (1.33, 0, 0), "CD": (0, 0, 0), "OE1": (-0.6, 1.1, 0)})
        hydrogens = survey.planar_hydrogens(gln, "NE2")
        self.assertEqual(len(hydrogens), 2)
        for h in hydrogens:
            self.assertAlmostEqual(float(np.linalg.norm(h - gln["NE2"].coord)), survey.NH_LENGTH, places=6)
            self.assertAlmostEqual(self.angle(gln["CD"].coord, gln["NE2"].coord, h), 120.0, places=4)
            self.assertAlmostEqual(float(h[2]), 0.0, places=6)              # in the amide plane
        self.assertGreater(float(np.linalg.norm(hydrogens[0] - hydrogens[1])), 1.5)
        self.assertIsNone(survey.planar_hydrogens(self.stub("SER", {"OG": (0, 0, 0)}), "OG"))      # hydroxyl H rotates: not placed

    def test_1ag6_cys84_thiolate_takes_the_asn38_amide(self):
        rows = pd.DataFrame(survey.sulfur_contacts(CIF, "A", min_separation=2))
        donors = rows[rows["donor"] != "sulfur"]
        cys = donors[(donors["kind"] == "Cys metal-bound") & (donors["donor_atom"] == "ASN N")]
        self.assertEqual(len(cys), 1)
        self.assertAlmostEqual(float(cys["distance"].iloc[0]), 3.446, places=2)
        self.assertTrue(bool(survey.hbond_like(cys).iloc[0]))
        sulfurs = rows[rows["donor"] == "sulfur"]["kind"].tolist()
        self.assertEqual(sorted(sulfurs), ["Cys metal-bound", "Met free", "Met metal-bound"])   # Cys84, a second Met, Met92 (2.88 A from Cu)
        self.assertFalse(((donors["distance"] < 3.0) & (donors["kind"] == "Cys metal-bound")).any())
        ser85 = donors[(donors["donor_atom"] == "SER N") & (donors["kind"] == "Cys metal-bound")]
        self.assertTrue(ser85.empty)                                   # the i+1 neighbor is excluded by the sequence rule (TODO D)

    def test_the_next_residue_amide_is_a_neighbor_whose_nh_does_not_point_at_the_sulfur(self):
        rows = pd.DataFrame(survey.sulfur_contacts(CIF, "A"))                       # min_separation=1 includes the neighbors
        ser85 = rows[(rows["kind"] == "Cys metal-bound") & (rows["donor_atom"] == "SER N")]
        self.assertEqual(ser85["offset"].tolist(), [1])
        self.assertAlmostEqual(float(ser85["distance"].iloc[0]), 3.58, places=2)
        self.assertFalse(bool(survey.hbond_like(ser85).iloc[0]))                    # H...S 3.63 A at 79 degrees


if __name__ == "__main__":
    unittest.main()
