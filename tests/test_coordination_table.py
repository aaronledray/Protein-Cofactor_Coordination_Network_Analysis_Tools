import unittest
from pathlib import Path

import pandas as pd

from modules import coordination_table as ct

SURVEY = Path(__file__).resolve().parents[1] / "docs/scope_checks/02_cu_type1_survey_proteins.csv"


class CoordinationTableTests(unittest.TestCase):
    def test_plastocyanin_ligand_distances(self):
        self.assertEqual(ct.coordination_tier("CU", "HIS", "ND1", 1.99).tier, "bond")
        self.assertEqual(ct.coordination_tier("CU", "CYS", "SG", 2.15).tier, "bond")
        self.assertEqual(ct.coordination_tier("CU", "MET", "SD", 2.88).tier, "weak")

    def test_boundaries_are_inclusive_and_beyond_is_none(self):
        self.assertEqual(ct.coordination_tier("CU", "MET", "SD", 2.7).tier, "bond")
        self.assertEqual(ct.coordination_tier("CU", "MET", "SD", 2.71).tier, "weak")
        self.assertEqual(ct.coordination_tier("CU", "MET", "SD", 3.3).tier, "weak")
        self.assertEqual(ct.coordination_tier("CU", "MET", "SD", 3.31).tier, "none")
        self.assertEqual(ct.coordination_tier("CU", "HIS", "NE2", 2.41).tier, "none")   # no weak range for His
        self.assertEqual(ct.coordination_tier("CU", "CYS", "SG", 2.51).tier, "none")

    def test_selenomethionine_counts_as_methionine(self):
        self.assertEqual(ct.coordination_tier("CU", "MSE", "SE", 2.6).tier, "bond")

    def test_unlisted_pairs_follow_the_global_cutoff_and_say_so(self):
        for args, tier in ((("ZN", "HIS", "NE2", 2.0), "default"), (("ZN", "HIS", "NE2", 3.0), "none"),
                           (("CU", "PRO", "O", 2.7), "default"),     # backbone O is not a Cu donor class yet
                           (("CU", "ASN", "N", 2.5), "default")):   # backbone N is never bond or weak
            call = ct.coordination_tier(*args)
            self.assertEqual(call.tier, tier)
            self.assertEqual(call.donor_class, "")
            self.assertIn("not in table", call.rule)
        self.assertEqual(ct.coordination_tier("CU", "PRO", "O", 2.9, default_max=3.0).tier, "default")

    def test_rule_text_states_the_cutoffs(self):
        self.assertEqual(ct.coordination_tier("CU", "MET", "SD", 2.88).rule, "Cu-Met S: bond <= 2.7, weak <= 3.3")

    def test_strongest_prefers_tier_then_distance(self):
        weak = (2.9, ct.coordination_tier("CU", "MET", "SD", 2.9))
        bond = (2.1, ct.coordination_tier("CU", "CYS", "SG", 2.1))
        self.assertEqual(ct.strongest([weak, bond]).donor_class, "Cys S")
        self.assertEqual(ct.strongest([]).tier, "none")

    def test_rows_cover_the_survey_they_cite(self):
        survey = pd.read_csv(SURVEY)
        by_class = {r.donor_class: r for r in ct.RULES if r.metal == "CU"}
        self.assertLessEqual(survey["cys_s"].max(), by_class["Cys S"].bond_max)
        self.assertLessEqual(survey["his_n_far"].max(), by_class["His N"].bond_max)
        axial = survey["met_s_axial"].dropna()
        self.assertGreaterEqual(by_class["Met S"].weak_max, axial.quantile(0.9))     # the weak range reaches most axial Met
        self.assertGreater(axial.max(), by_class["Met S"].weak_max)                   # and is a convention: some sit beyond it


if __name__ == "__main__":
    unittest.main()
