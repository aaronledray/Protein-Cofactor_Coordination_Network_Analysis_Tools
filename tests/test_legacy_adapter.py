"""Tests for mapping legacy CSV artifacts into canonical analysis tables."""

import json
from pathlib import Path
import unittest

from modules.legacy_adapter import legacy_csvs_to_analysis, legacy_csvs_to_signature


ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = ROOT / "tests/fixtures/baseline/1ag6"


class LegacyAdapterTests(unittest.TestCase):
    def test_mixed_legacy_breakdown_maps_to_current_table_contract(self):
        tables = legacy_csvs_to_analysis(
            FIXTURE_DIR / "Coord_Breakdown.csv",
            FIXTURE_DIR / "Coord_Links.csv",
            structure_id="1ag6",
        )

        self.assertEqual(set(tables), {"atoms", "residues", "links", "contacts"})
        self.assertGreater(len(tables["atoms"]), 0)
        self.assertEqual(set(tables["atoms"]["shell"]), {"Cofactor", "PCS", "SCS"})
        self.assertEqual(len(tables["links"]), 7)
        self.assertEqual(len(tables["contacts"]), 7)
        self.assertEqual(set(tables["contacts"]["contact_role"]), {
            "legacy_primary_contact",
            "legacy_secondary_contact",
        })

    def test_legacy_signature_preserves_geometric_primary_provenance(self):
        signature = legacy_csvs_to_signature(
            FIXTURE_DIR / "Coord_Breakdown.csv",
            FIXTURE_DIR / "Coord_Links.csv",
            structure_id="1ag6",
        )

        self.assertEqual(signature["schema_version"], "1.0")
        self.assertEqual(signature["cofactor"]["key"], "CU")
        self.assertEqual(len(signature["edges"]["primary"]), 5)
        self.assertEqual(len(signature["edges"]["secondary"]), 2)
        self.assertEqual(signature["source"]["format"], "legacy_csv")
        json.dumps(signature)


if __name__ == "__main__":
    unittest.main()
