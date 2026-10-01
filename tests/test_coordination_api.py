"""Tests for the side-effect-free and batch coordination APIs."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from modules.coordination_api import analyze_structure, batch_analyze


ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"
OEX = ROOT / "reference_structures/1_OEX/0_Aligned_Reduced/4ub6.pdb"


class CoordinationApiTests(unittest.TestCase):
    def test_analyze_structure_returns_tables_without_files(self):
        with TemporaryDirectory(prefix="sscna-api-") as work_dir:
            with redirect_stdout(StringIO()):
                result = analyze_structure(
                    PLASTOCYANIN,
                    "CU",
                    exclude_moieties=["alanine_sidechain"],
                )

            self.assertEqual(set(result), {"residues", "atoms", "links"})
            self.assertEqual(len(result["residues"]), 8)
            self.assertEqual(len(result["atoms"]), 8)
            self.assertEqual(len(result["links"]), 7)
            self.assertEqual(
                list(result["residues"].columns),
                [
                    "structure_id", "site_id", "cofactor_residue_name",
                    "cofactor_residue_number", "cofactor_chain",
                    "cofactor_insertion_code", "shell", "residue_name",
                    "residue_number", "chain", "insertion_code", "hetero_flag",
                    "atoms_involved", "minimum_distance_A",
                ],
            )
            self.assertEqual(result["residues"]["site_id"].unique().tolist(), ["all"])
            self.assertEqual(
                result["residues"]
                .loc[result["residues"]["shell"] == "PCS", "minimum_distance_A"]
                .isna()
                .sum(),
                0,
            )
            self.assertEqual(list(Path(work_dir).iterdir()), [])

    def test_batch_analyze_isolates_bad_inputs_and_combines_results(self):
        result = batch_analyze(
            [PLASTOCYANIN, OEX, ROOT / "does-not-exist.pdb"],
            "CU,OEX",
            workers=2,
            exclude_moieties=["alanine_sidechain"],
        )

        self.assertEqual(len(result["errors"]), 1)
        self.assertEqual(
            set(result["residues"]["structure_id"]),
            {"1ag6", "4ub6"},
        )
        self.assertEqual(set(result["links"]["structure_id"]), {"1ag6", "4ub6"})

    def test_shells_three_adds_tcs_and_adjacent_links(self):
        with redirect_stdout(StringIO()):
            result = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
                shells=3,
            )

        self.assertIn("TCS", set(result["residues"]["shell"]))
        self.assertIn("scs->tcs", set(result["links"]["link_type"]))

    def test_per_site_mode_assigns_stable_site_ids(self):
        with redirect_stdout(StringIO()):
            result = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
                site_mode="per-site",
            )

        self.assertEqual(set(result["residues"]["site_id"]), {"site_1"})
        self.assertEqual(set(result["atoms"]["site_id"]), {"site_1"})
        self.assertIn("site_id", result["links"].columns)

    def test_definition_options_are_opt_in(self):
        with redirect_stdout(StringIO()):
            default = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
            )
            carbon = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
                include_carbon_seeds=True,
            )
            direct = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
                direct_coordination=True,
            )

        default_pcs = len(default["atoms"].query("shell == 'PCS'"))
        carbon_pcs = len(carbon["atoms"].query("shell == 'PCS'"))
        self.assertEqual(default_pcs, 5)
        self.assertGreater(carbon_pcs, default_pcs)
        self.assertEqual(int(default["links"]["direct_coordination"].sum()), 0)
        self.assertEqual(int(direct["links"]["direct_coordination"].sum()), 3)


if __name__ == "__main__":
    unittest.main()
