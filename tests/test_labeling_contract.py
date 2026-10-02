"""Golden-fixture tests that freeze labeling contract 1.0.

If one of these fails, either a regression changed shell labels or the
contract changed deliberately. In the latter case bump
``LABELING_CONTRACT_VERSION`` and regenerate tests/fixtures/labels/.
"""

import json
from pathlib import Path
import unittest

import pandas as pd

from modules.ml_export import (
    CONTRACT_DEFAULTS,
    LABELING_CONTRACT_VERSION,
    MixedCofactorError,
    analyze_labeling_example,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/fixtures/labels"
CASES = {
    "1ag6": (ROOT / "reference_structures/0_Plastocyanin/1ag6.cif", "CU"),
    "1a6m": (ROOT / "reference_structures/3_Myoglobin/1a6m.pdb", "HEM"),
    "1ca2": (ROOT / "reference_structures/3_Metals/1ca2.pdb", "ZN"),
}


class LabelingContractTests(unittest.TestCase):
    def test_version_and_defaults_are_frozen(self):
        self.assertEqual(LABELING_CONTRACT_VERSION, "1.0")
        self.assertEqual(
            CONTRACT_DEFAULTS,
            {
                "distance_cutoff": 3.6, "shells": 3, "include_carbon_seeds": False,
                "direct_coordination": True, "direct_coordination_cutoff": 2.6,
                "expand_residues": False, "site_mode": "per-site",
                "site_model_mode": "per-model",
            },
        )

    def test_labels_and_contract_match_golden_fixtures(self):
        for name, (path, cofactor) in CASES.items():
            with self.subTest(structure=name):
                labels, contract = analyze_labeling_example(path, cofactor)
                golden = pd.read_csv(FIXTURES / f"{name}_labels.csv", keep_default_na=False)
                pd.testing.assert_frame_equal(
                    labels.astype(str).reset_index(drop=True),
                    golden.astype(str).reset_index(drop=True),
                )
                expected = json.loads((FIXTURES / f"{name}_contract.json").read_text())
                self.assertEqual(json.loads(json.dumps(contract)), expected)

    def test_mixed_cofactor_families_are_rejected_unless_allowed(self):
        path = CASES["1a6m"][0]
        with self.assertRaises(MixedCofactorError):
            analyze_labeling_example(path, ["HEM", "ZN"])
        _, contract = analyze_labeling_example(path, ["HEM", "ZN"], allow_mixed_cofactors=True)
        self.assertTrue(contract["mixed_cofactor_classes"])

    def test_policy_overrides_are_refused(self):
        path, cofactor = CASES["1ag6"]
        for option in ({"site_mode": "union"}, {"first_model_only": True}):
            with self.assertRaises(ValueError):
                analyze_labeling_example(path, cofactor, **option)

    def test_overrides_are_recorded_in_contract(self):
        path, cofactor = CASES["1ag6"]
        _, contract = analyze_labeling_example(path, cofactor, distance_cutoff=3.2)
        self.assertEqual(contract["analysis_parameters"]["distance_cutoff"], 3.2)
        self.assertEqual(contract["effective_cutoff_A"], 3.2)

    def test_exclude_solvent_is_opt_in_and_recorded(self):
        path, cofactor = CASES["1a6m"]
        default, default_contract = analyze_labeling_example(path, cofactor)
        dry, contract = analyze_labeling_example(path, cofactor, exclude_solvent=True)
        self.assertEqual((default["residue_name"] == "HOH").sum(), 18)
        self.assertEqual(len(dry), len(default) - 18)
        self.assertFalse((dry["residue_name"] == "HOH").any())
        self.assertTrue(contract["exclude_solvent"])
        self.assertNotIn("exclude_solvent", default_contract)
        # non-water labels are unchanged by the filter
        kept = default[default["residue_name"] != "HOH"].reset_index(drop=True)
        self.assertTrue(kept.equals(dry.reset_index(drop=True)))

    def test_class_cutoff_changes_effective_cutoff(self):
        path, cofactor = CASES["1ag6"]
        _, contract = analyze_labeling_example(
            path, cofactor, cofactor_class_cutoffs={"metal": 2.8}
        )
        self.assertEqual(contract["effective_cutoff_A"], 2.8)


if __name__ == "__main__":
    unittest.main()
