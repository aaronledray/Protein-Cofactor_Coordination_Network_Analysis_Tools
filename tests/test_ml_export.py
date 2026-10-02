"""Tests for the deterministic ML-facing residue label export."""

from pathlib import Path
import unittest

import pandas as pd

from modules.coordination_api import analyze_structure
from modules.ml_export import (
    LABEL_COLUMNS,
    LABELING_CONTRACT_VERSION,
    export_residue_labels,
    labeling_contract,
    shell_depth,
)


ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


class MlExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables = analyze_structure(
            PLASTOCYANIN, "CU", shells=3, direct_coordination=True, site_mode="per-site"
        )

    def test_labels_match_residue_table_and_schema(self):
        labels = export_residue_labels(self.tables)
        self.assertEqual(list(labels.columns), LABEL_COLUMNS)
        self.assertTrue((labels["contract_version"] == LABELING_CONTRACT_VERSION).all())
        residues = self.tables["residues"]
        residues = residues[residues["shell"] != "Cofactor"]
        self.assertEqual(len(labels), len(residues))
        primary = labels[labels["is_primary"]]
        self.assertEqual(
            sorted(primary["residue_number"].astype(int)), [36, 37, 84, 87, 92]
        )
        self.assertEqual(set(primary["shell"]), {"PCS"})

    def test_direct_coordination_and_motifs_are_labeled(self):
        labels = export_residue_labels(self.tables).set_index("residue_number")
        labels.index = labels.index.astype(int)
        self.assertTrue(labels.loc[37, "direct_coordination"])
        self.assertEqual(labels.loc[37, "motifs"], "imidazole")
        self.assertEqual(labels.loc[84, "motifs"], "thiol")
        self.assertEqual(labels.loc[37, "cofactor_residue_name"], "CU")

    def test_export_is_deterministic_and_input_order_independent(self):
        first = export_residue_labels(self.tables)
        shuffled = {
            key: frame.sample(frac=1, random_state=3).reset_index(drop=True)
            for key, frame in self.tables.items()
        }
        pd.testing.assert_frame_equal(first, export_residue_labels(shuffled))

    def test_selectors_and_errors(self):
        labels = export_residue_labels(self.tables, model_id=0, site_id="site_1")
        self.assertEqual(set(labels["site_id"]), {"site_1"})
        with self.assertRaises(ValueError):
            export_residue_labels(self.tables, model_id=99)
        with self.assertRaises(ValueError):
            export_residue_labels(self.tables, site_id="missing")

    def test_multiple_pooled_models_are_rejected(self):
        tables = {key: frame.copy() for key, frame in self.tables.items()}
        extra = tables["atoms"].copy()
        extra["model_id"] = 1
        extra["site_id"] = "all"
        tables["atoms"] = pd.concat([tables["atoms"].assign(site_id="all"), extra])
        with self.assertRaises(ValueError):
            export_residue_labels(tables)
        self.assertFalse(export_residue_labels(tables, model_id=1).empty)

    def test_contract_records_parameters_and_shell_depth(self):
        contract = labeling_contract(distance_cutoff=3.6, shells=3)
        self.assertEqual(contract["contract_version"], LABELING_CONTRACT_VERSION)
        self.assertEqual(contract["analysis_parameters"]["shells"], 3)
        self.assertEqual([shell_depth(s) for s in ("PCS", "SCS", "TCS", "Shell4")], [1, 2, 3, 4])
        with self.assertRaises(ValueError):
            shell_depth("bogus")

    def test_empty_tables_return_empty_labels(self):
        empty = analyze_structure(PLASTOCYANIN, "ZZZ", site_mode="per-site")
        self.assertEqual(list(export_residue_labels(empty).columns), LABEL_COLUMNS)


if __name__ == "__main__":
    unittest.main()
