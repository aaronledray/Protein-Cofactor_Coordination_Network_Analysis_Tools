"""Offline checks for scripts/survey_cu_type1.py (no network): site recognition and distances on 1AG6."""

import importlib.util
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CIF = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"
spec = importlib.util.spec_from_file_location("survey_cu_type1", ROOT / "scripts/survey_cu_type1.py")
survey = importlib.util.module_from_spec(spec)
spec.loader.exec_module(survey)


class SurveyTests(unittest.TestCase):
    def test_1ag6_is_one_type1_site_with_the_axial_methionine(self):
        sites = survey.type1_sites(CIF)
        self.assertEqual(len(sites), 1)
        row = survey.characterise(sites[0])
        self.assertAlmostEqual(row["cys_s"], 2.15, delta=0.01)
        self.assertAlmostEqual(row["met_s_axial"], 2.88, delta=0.01)       # Met92 is His87 + 5
        self.assertIsNone(row["met_s_other"])
        self.assertLess(row["his_n_near"], 2.1)
        self.assertLess(row["his_n_far"], 2.1)

    def test_a_met_at_another_offset_is_not_called_axial(self):
        site = {"pdb": "TEST", "chain": "A", "donors": [
            {"resname": "CYS", "atom": "SG", "chain": "A", "resnum": 84, "distance": 2.2},
            {"resname": "HIS", "atom": "ND1", "chain": "A", "resnum": 37, "distance": 2.0},
            {"resname": "HIS", "atom": "ND1", "chain": "A", "resnum": 87, "distance": 2.0},
            {"resname": "MET", "atom": "SD", "chain": "A", "resnum": 120, "distance": 3.0}]}
        row = survey.characterise(site)
        self.assertIsNone(row["met_s_axial"])
        self.assertEqual(row["met_s_other"], 3.0)

    def test_per_protein_takes_the_median_over_entries(self):
        rows = [{"pdb": p, "chain": "A", "cys_s": c, "his_n_near": 2.0, "his_n_far": 2.1, "met_s_axial": m,
                 "met_s_other": None, "amide_o": None, "backbone_o": None}
                for p, c, m in (("1AAA", 2.1, 2.7), ("1BBB", 2.3, 2.9), ("1CCC", 2.2, None))]
        groups = {(p, "A"): {"protein": "P1" if p != "1CCC" else "P2", "name": "Plastocyanin", "organism": "x"}
                  for p in ("1AAA", "1BBB", "1CCC")}
        table = survey.per_protein(rows, groups).set_index("protein")
        self.assertEqual(len(table), 2)
        self.assertAlmostEqual(table.loc["P1", "cys_s"], 2.2)
        self.assertAlmostEqual(table.loc["P1", "met_s_axial"], 2.8)
        self.assertTrue(pd.isna(table.loc["P2", "met_s_axial"]))
        self.assertEqual(table.loc["P1", "family"], "plastocyanin")


if __name__ == "__main__":
    unittest.main()
