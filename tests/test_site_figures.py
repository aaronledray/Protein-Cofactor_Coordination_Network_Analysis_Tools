"""The standard human-check figures for a metal site (plastocyanin 1AG6 as the fixture)."""

import tempfile
import unittest
from pathlib import Path

from modules import site_figures

CIF = Path(__file__).resolve().parents[1] / "reference_structures/0_Plastocyanin/1ag6.cif"


class SiteFigureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.site = site_figures.collect_site(CIF, "CU")

    def test_roles_follow_the_distance_and_shell_rules(self):
        atoms = self.site["atoms"]
        role = {(r.position, r.atom): r.role for r in atoms.itertuples() if r.element in ("N", "O", "S")}
        for key in ((37, "ND1"), (84, "SG"), (87, "ND1")):
            self.assertEqual(role[key], "coordination")
        self.assertEqual(role[(92, "SD")], "weak")                 # 2.88 A: beyond the Cu-Met bond cutoff (2.7), inside its weak cutoff (3.3)
        self.assertEqual(role[(36, "O")], "primary")
        self.assertEqual(role[(38, "N")], "secondary")
        self.assertEqual(role[(85, "N")], "secondary")

    def test_sulfur_table_reports_how_the_contact_rule_treats_each_partner(self):
        table = site_figures.sulfur_partners(self.site).set_index(["sulfur", "partner"])
        asn38 = table.loc[("Cys84 SG", "Asn38 N")]
        self.assertAlmostEqual(asn38["distance"], 3.446, places=2)
        self.assertEqual(asn38["pair_listed"], "listed")
        # polar needs an N/O pair, so a contact to a sulfur is never typed polar; the opt-in sulfur rule types it.
        self.assertFalse(asn38["atom_pair_polar"])
        self.assertTrue(asn38["atom_pair_polar_s"])
        self.assertEqual(table.loc[("Cys84 SG", "Ser85 N"), "pair_listed"], "excluded: sequence neighbor")

    def test_writes_the_four_figures_and_the_checklist(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = site_figures.site_check_figures(CIF, "CU", tmp)
            for name, path in paths.items():
                self.assertTrue(path.exists() and path.stat().st_size > 1000, name)
            text = paths["checklist"].read_text()
            self.assertIn("Cys84 SG", text)
            self.assertIn("Question:", text)
            self.assertIn("## Rules in force", text)                   # the cutoffs behind every number are written next to it
            for name in ("polar_s_geom", "coord_tier: Cu Met S", "hydrophobic"):
                self.assertIn(name, text)

    def test_non_metal_cofactor_is_rejected(self):
        with self.assertRaises(Exception):
            site_figures.collect_site(CIF, "HOH")


if __name__ == "__main__":
    unittest.main()
