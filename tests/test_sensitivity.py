"""Tests for label-sensitivity comparison."""

from pathlib import Path
import unittest

import pandas as pd

from modules.coordination_api import analyze_structure
from modules.ml_export import export_residue_labels
from modules.sensitivity import compare_labels, sensitivity_table

PLASTOCYANIN = Path(__file__).resolve().parents[1] / "reference_structures/0_Plastocyanin/1ag6.cif"


class SensitivityTests(unittest.TestCase):
    def test_identical_labels_are_fully_stable(self):
        labels = export_residue_labels(analyze_structure(PLASTOCYANIN, "CU", site_mode="per-site"))
        summary = compare_labels(labels, labels)
        self.assertEqual(summary["primary_jaccard"], 1.0)
        self.assertEqual((summary["shell_changed"], summary["gained"], summary["lost"]), (0, 0, 0))

    def test_cutoff_changes_membership_in_expected_direction(self):
        table = sensitivity_table(PLASTOCYANIN, "CU").set_index("variant")
        self.assertGreater(table.loc["cutoff_3.2", "lost"], 0)
        self.assertEqual(table.loc["cutoff_3.2", "gained"], 0)
        self.assertGreater(table.loc["cutoff_4.0", "gained"], 0)
        self.assertEqual(table.loc["cutoff_4.0", "lost"], 0)
        self.assertGreater(table.loc["carbon_seeds", "gained"], 0)

    def test_failing_variant_is_reported_not_raised(self):
        table = sensitivity_table(PLASTOCYANIN, "CU", variants={"bad": {"shells": 0}})
        self.assertIn("ValueError", table.set_index("variant").loc["bad", "error"])

    def test_empty_labels_compare_cleanly(self):
        empty = export_residue_labels(analyze_structure(PLASTOCYANIN, "ZZZ", site_mode="per-site"))
        self.assertEqual(compare_labels(empty, empty)["primary_jaccard"], 1.0)


if __name__ == "__main__":
    unittest.main()
