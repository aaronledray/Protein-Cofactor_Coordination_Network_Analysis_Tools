"""Tests for canonical network signatures and family-profile scoring."""

from contextlib import redirect_stdout
from copy import deepcopy
from io import StringIO
import json
from pathlib import Path
import unittest

from modules.coordination_api import analyze_structure
from modules.network_comparison import (
    build_network_profile,
    build_network_signature,
    compare_network_signatures,
    score_signature_against_profile,
)


ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


class NetworkComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with redirect_stdout(StringIO()):
            analysis = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
                shells=3,
            )
        cls.signature = build_network_signature(analysis, structure_id="1ag6")

    def test_signature_is_stable_and_json_friendly(self):
        signature = self.signature

        self.assertEqual(signature["schema_version"], "1.0")
        self.assertEqual(signature["cofactor"]["key"], "CU")
        self.assertEqual(len(signature["edges"]["primary"]), 3)
        self.assertEqual(len(signature["edges"]["secondary"]), 2)
        self.assertEqual(len(signature["edges"]["tertiary"]), 1)
        self.assertTrue(all("distance_A" in edge for edge in signature["edges"]["primary"]))
        self.assertTrue(all("residue_number" not in atom for atom in signature["atoms"]))

        # This is the contract the future upload/PDB-ID web endpoint can return.
        json.dumps(signature)

    def test_identical_signatures_score_as_identical(self):
        comparison = compare_network_signatures(self.signature, self.signature)

        self.assertTrue(comparison["cofactor_compatible"])
        self.assertEqual(comparison["matched_edge_count"], 6)
        self.assertEqual(comparison["scores"]["overall"], 1.0)
        self.assertEqual(comparison["scores"]["primary"], 1.0)
        self.assertEqual(comparison["unmatched"]["primary"]["reference_only"], [])

    def test_unmatched_edges_are_reported_for_duplicate_features(self):
        reference = {
            "structure_id": "reference",
            "cofactor": {"key": "CU"},
            "edges": {
                "primary": [{"match_key": "same"}, {"match_key": "same"}],
                "secondary": [],
                "tertiary": [],
            },
            "residues": [],
        }
        query = deepcopy(reference)
        query["structure_id"] = "query"
        query["edges"]["primary"] = [{"match_key": "same"}]

        comparison = compare_network_signatures(reference, query)

        self.assertEqual(len(comparison["unmatched"]["primary"]["reference_only"]), 1)
        self.assertEqual(comparison["unmatched"]["primary"]["query_only"], [])
        self.assertLess(comparison["scores"]["primary"], 1.0)

    def test_profile_preserves_feature_support_and_multiplicity(self):
        profile = build_network_profile([self.signature, deepcopy(self.signature)])

        self.assertEqual(profile["n_signatures"], 2)
        self.assertEqual(profile["dominant_cofactor_key"], "CU")
        primary_entries = profile["edges"]["primary"]
        self.assertTrue(all(entry["required"] for entry in primary_entries.values()))
        self.assertTrue(all(entry["support"] == 1.0 for entry in primary_entries.values()))
        self.assertTrue(any(entry["multiplicity"]["max"] == 2 for entry in primary_entries.values()))

        score = score_signature_against_profile(self.signature, profile)
        self.assertTrue(score["cofactor_compatible"])
        self.assertEqual(score["edges"]["primary"]["matched_required_count"], 2)
        self.assertEqual(score["residues"]["matched_required_count"], len(profile["residues"]))
        self.assertGreater(score["scores"]["residue"], 0.95)
        self.assertGreater(score["scores"]["overall"], 0.95)

    def test_profile_score_explains_missing_required_features(self):
        profile = build_network_profile([self.signature])
        query = deepcopy(self.signature)
        missing_key = query["edges"]["primary"][0]["match_key"]
        query["edges"]["primary"] = [
            edge for edge in query["edges"]["primary"]
            if edge["match_key"] != missing_key
        ]

        score = score_signature_against_profile(query, profile)

        self.assertLess(score["scores"]["overall"], 1.0)
        self.assertIn(
            missing_key,
            score["edges"]["primary"]["missing_required_features"],
        )


if __name__ == "__main__":
    unittest.main()
