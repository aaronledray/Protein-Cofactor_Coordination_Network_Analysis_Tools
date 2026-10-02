"""Tests for named cofactor families and cutoff configuration."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from modules.cofactor_classes import (
    load_cofactor_class_config,
    merge_cofactor_class_configs,
    resolve_cofactor_classes,
    resolve_effective_distance_cutoff,
)


class CofactorClassTests(unittest.TestCase):
    def test_reference_cofactors_resolve_to_named_families(self):
        self.assertEqual(resolve_cofactor_classes(["CU"]), ("metal_ion",))
        self.assertEqual(resolve_cofactor_classes(["HEM"]), ("heme",))
        self.assertEqual(resolve_cofactor_classes(["SF4"]), ("iron_sulfur_cluster",))
        self.assertEqual(
            resolve_cofactor_classes(["ICS", "CLF", "HCA"]),
            ("metallo_cluster",),
        )

    def test_legacy_flat_cutoffs_and_aliases_remain_supported(self):
        self.assertEqual(
            resolve_effective_distance_cutoff(["CU"], 3.6, {"metal": 2.8}),
            2.8,
        )
        self.assertEqual(
            resolve_effective_distance_cutoff(["HEM"], 3.6, {"heme": 3.2}),
            3.2,
        )
        self.assertEqual(
            resolve_effective_distance_cutoff(["HEM"], 3.6, {"HEM": 3.1}),
            3.1,
        )
        self.assertEqual(
            resolve_effective_distance_cutoff(["HEM"], 3.6, {"organic": 3.1}),
            3.1,
        )

    def test_structured_rules_support_custom_residue_families(self):
        config = {
            "cofactor_classes": {
                "custom_cluster": {
                    "residues": ["ABC", "ABD"],
                    "cutoff_A": 2.4,
                }
            }
        }
        self.assertEqual(
            resolve_cofactor_classes(["ABC"], config),
            ("custom_cluster",),
        )
        self.assertEqual(
            resolve_effective_distance_cutoff(["ABC"], 3.6, config),
            2.4,
        )

    def test_mixed_families_keep_fallback_for_unconfigured_classes(self):
        self.assertEqual(
            resolve_effective_distance_cutoff(["CU", "HEM"], 3.6, {"metal": 2.8}),
            3.6,
        )
        self.assertEqual(
            resolve_effective_distance_cutoff(
                ["CU", "HEM"],
                3.6,
                {"metal": 2.8, "heme": 3.1},
            ),
            3.1,
        )

    def test_json_rules_load_and_merge_with_cli_overrides(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "cofactor_classes.json"
            path.write_text(
                json.dumps({"cofactor_classes": {"heme": {"cutoff_A": 3.4}}}),
                encoding="utf-8",
            )
            loaded = load_cofactor_class_config(path)
            merged = merge_cofactor_class_configs(loaded, {"heme": 3.2})
            self.assertEqual(
                resolve_effective_distance_cutoff(["HEM"], 3.6, merged),
                3.2,
            )


if __name__ == "__main__":
    unittest.main()
