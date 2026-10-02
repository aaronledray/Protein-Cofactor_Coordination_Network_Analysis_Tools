"""Tests for file-oriented pairwise and profile comparison workflows."""

from pathlib import Path
import unittest

from modules.comparison_runner import (
    collect_signatures,
    compare_reference_to_queries,
    profile_reference_set,
)


ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


class ComparisonRunnerTests(unittest.TestCase):
    def test_pairwise_runner_attaches_paths_and_sorts_results(self):
        result = compare_reference_to_queries(
            PLASTOCYANIN,
            [PLASTOCYANIN],
            "CU",
            analysis_options={
                "exclude_moieties": ["alanine_sidechain"],
                "shells": 3,
            },
        )

        self.assertEqual(result["mode"], "pairwise")
        self.assertEqual(len(result["comparisons"]), 1)
        comparison = result["comparisons"][0]
        self.assertEqual(comparison["reference_path"], str(PLASTOCYANIN))
        self.assertEqual(comparison["query_path"], str(PLASTOCYANIN))
        self.assertEqual(comparison["scores"]["overall"], 1.0)

    def test_profile_runner_scores_queries_and_isolates_bad_inputs(self):
        result = profile_reference_set(
            [PLASTOCYANIN],
            "CU",
            query_inputs=[PLASTOCYANIN, ROOT / "missing-query.pdb"],
            analysis_options={
                "exclude_moieties": ["alanine_sidechain"],
                "shells": 3,
            },
        )

        self.assertEqual(result["profile"]["n_signatures"], 1)
        self.assertEqual(len(result["scores"]), 1)
        self.assertEqual(result["scores"][0]["scores"]["residue"], 1.0)
        self.assertEqual(len(result["errors"]), 1)
        self.assertIn("missing-query.pdb", result["errors"][0]["structure_path"])

    def test_profile_template_supplies_map_numbering_anchor(self):
        result = profile_reference_set(
            [PLASTOCYANIN],
            "CU",
            template_path=PLASTOCYANIN,
            analysis_options={
                "exclude_moieties": ["alanine_sidechain"],
                "shells": 3,
            },
        )

        self.assertEqual(result["template"]["signature"]["structure_id"], "1ag6")
        self.assertTrue(result["template"]["position_order"])
        self.assertEqual(result["references"][0]["signature"]["structure_id"], "1ag6")

    def test_collect_signatures_accepts_directories(self):
        records, errors = collect_signatures(
            [ROOT / "reference_structures/0_Plastocyanin"],
            "CU",
            analysis_options={"exclude_moieties": ["alanine_sidechain"]},
        )

        self.assertEqual({record["signature"]["structure_id"] for record in records}, {"1ag6"})
        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main()
