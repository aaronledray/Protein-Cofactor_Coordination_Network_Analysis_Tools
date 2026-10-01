"""Known-site checks used to document expected chemistry."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import unittest

from modules.coordination_api import analyze_structure


ROOT = Path(__file__).resolve().parents[1]


def analyze(path, cofactor):
    with redirect_stdout(StringIO()):
        return analyze_structure(path, cofactor, exclude_moieties=["alanine_sidechain"])


class ChemistryValidationTests(unittest.TestCase):
    def test_plastocyanin_known_primary_ligands(self):
        result = analyze(ROOT / "reference_structures/0_Plastocyanin/1ag6.cif", "CU")
        pcs = set(
            zip(
                result["residues"].query("shell == 'PCS'")["residue_name"],
                result["residues"].query("shell == 'PCS'")["residue_number"],
            )
        )
        self.assertTrue({("HIS", 37), ("CYS", 84), ("HIS", 87), ("MET", 92)} <= pcs)

    def test_myoglobin_proximal_and_distal_histidines(self):
        result = analyze(ROOT / "reference_structures/3_Myoglobin/1a6m.pdb", "HEM")
        shells = {
            (row.residue_name, row.residue_number): row.shell
            for row in result["residues"].itertuples()
        }
        self.assertEqual(shells[("HIS", 93)], "PCS")
        self.assertEqual(shells[("HIS", 64)], "SCS")

    def test_single_zinc_site_has_expected_histidine_ligands(self):
        result = analyze(ROOT / "reference_structures/3_Metals/1ca2.pdb", "ZN")
        pcs = set(
            zip(
                result["residues"].query("shell == 'PCS'")["residue_name"],
                result["residues"].query("shell == 'PCS'")["residue_number"],
            )
        )
        self.assertTrue({("HIS", 94), ("HIS", 96), ("HIS", 119)} <= pcs)

    def test_sf4_reference_runs_and_finds_cysteine_ligands(self):
        result = analyze(ROOT / "reference_structures/3_Metals/2zvs.pdb", "SF4")
        self.assertEqual(len(result["residues"].query("shell == 'Cofactor'")), 6)
        self.assertGreaterEqual(
            len(result["residues"].query("shell == 'PCS' and residue_name == 'CYS'")),
            12,
        )


if __name__ == "__main__":
    unittest.main()
