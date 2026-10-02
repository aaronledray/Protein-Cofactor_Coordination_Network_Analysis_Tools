"""Tests for the opt-in model_id column in the residue table."""

from pathlib import Path
import tempfile
import unittest

from modules.coordination_api import (
    RESIDUE_COLUMNS,
    RESIDUE_COLUMNS_WITH_MODEL,
    analyze_structure,
    batch_analyze,
)

ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


def _line(serial, name, residue, chain, number, x, element, record="ATOM  "):
    return (
        f"{record}{serial:5d} {name:>4s} {residue:>3s} {chain}{number:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00          {element:>2s}\n"
    )


TWO_MODELS = (
    "MODEL        1\n"
    + _line(1, "CU", "CU", "A", 200, 0.0, "CU", "HETATM")
    + _line(2, "ND1", "HIS", "A", 52, 2.0, "N")
    + "ENDMDL\nMODEL        2\n"
    + _line(1, "CU", "CU", "A", 200, 0.0, "CU", "HETATM")
    + _line(2, "ND1", "HIS", "A", 52, 2.5, "N")
    + "ENDMDL\nEND\n"
)


class ResidueModelIdTests(unittest.TestCase):
    def test_default_schema_is_unchanged(self):
        residues = analyze_structure(PLASTOCYANIN, "CU")["residues"]
        self.assertEqual(list(residues.columns), RESIDUE_COLUMNS)
        self.assertNotIn("model_id", residues.columns)

    def test_opt_in_adds_column_without_changing_single_model_values(self):
        plain = analyze_structure(PLASTOCYANIN, "CU", shells=3)["residues"]
        tagged = analyze_structure(PLASTOCYANIN, "CU", shells=3, include_model_id=True)["residues"]
        self.assertEqual(list(tagged.columns), RESIDUE_COLUMNS_WITH_MODEL)
        self.assertEqual(set(tagged["model_id"]), {0})
        self.assertTrue(plain.equals(tagged.drop(columns="model_id")))

    def test_per_model_sites_carry_their_own_model_and_distances(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "two_models.pdb"
            path.write_text(TWO_MODELS)
            tagged = analyze_structure(
                path, "CU", site_mode="per-site", site_model_mode="per-model",
                include_model_id=True,
            )["residues"]
        his = tagged[tagged["residue_name"] == "HIS"].set_index("model_id")
        self.assertEqual(sorted(his.index), [0, 1])
        self.assertAlmostEqual(his.loc[0, "minimum_distance_A"], 2.0, places=3)
        self.assertAlmostEqual(his.loc[1, "minimum_distance_A"], 2.5, places=3)
        self.assertEqual(his.loc[0, "site_id"], "site_1_model_0")
        self.assertEqual(his.loc[1, "site_id"], "site_2_model_1")
        cofactors = tagged[tagged["shell"] == "Cofactor"]
        self.assertEqual(sorted(cofactors["model_id"]), [0, 1])

    def test_pooled_mode_keeps_only_the_nearest_copy_across_models(self):
        """Pinned legacy semantics: pooled mode picks one representative atom per
        residue identity and moiety across all models (the closest to the pooled
        cofactor), so a same-identity residue in a later model is not re-seeded.
        Use site_model_mode='per-model' to keep each model's residues."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "two_models.pdb"
            path.write_text(TWO_MODELS)
            pooled = analyze_structure(path, "CU", include_model_id=True)["residues"]
        his = pooled[pooled["residue_name"] == "HIS"]
        self.assertEqual(his["model_id"].tolist(), [0])
        self.assertAlmostEqual(his["minimum_distance_A"].iloc[0], 2.0, places=3)

    def test_empty_and_batch_paths_honor_the_flag(self):
        empty = analyze_structure(PLASTOCYANIN, "ZZZ", site_mode="per-site", include_model_id=True)
        self.assertEqual(list(empty["residues"].columns), RESIDUE_COLUMNS_WITH_MODEL)
        batch = batch_analyze([PLASTOCYANIN], "CU", include_model_id=True)
        self.assertIn("model_id", batch["residues"].columns)
        none = batch_analyze([PLASTOCYANIN.with_name("missing.cif")], "CU", include_model_id=True)
        self.assertEqual(list(none["residues"].columns), RESIDUE_COLUMNS_WITH_MODEL)


if __name__ == "__main__":
    unittest.main()
