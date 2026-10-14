"""Tests for the opt-in exclusion of deposited hydrogens."""

from pathlib import Path
import tempfile
import unittest

from modules.coordination_api import analyze_structure

ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


def _line(serial, name, residue, number, x, element, record="ATOM  "):
    return (f"{record}{serial:5d} {name:>4s} {residue:>3s} A{number:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00          {element:>2s}\n")


# Cu at the origin; His ND1 and a hydrogen of a neighboring Cys within 3.6 A.
WITH_HYDROGENS = (
    _line(1, "CU", "CU", 200, 0.0, "CU", "HETATM")
    + _line(2, "ND1", "HIS", 52, 2.0, "N")
    + _line(3, "HD1", "HIS", 52, 2.9, "H")
    + _line(4, "HG", "CYS", 84, 3.0, "H")
    + _line(5, "SG", "CYS", 84, 5.5, "S")
    + "END\n"
)


class ExcludeHydrogenTests(unittest.TestCase):
    def _atoms(self, **options):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "h.pdb"
            path.write_text(WITH_HYDROGENS)
            return analyze_structure(path, "CU", **options)["atoms"]

    def test_default_keeps_hydrogens_as_possible_seeds(self):
        shell = self._atoms()
        self.assertIn("H", set(shell[shell["shell"] == "PCS"]["element"]))

    def test_opt_in_removes_hydrogens_before_shells_are_built(self):
        atoms = self._atoms(exclude_hydrogens=True)
        self.assertNotIn("H", set(atoms["element"]))
        self.assertEqual(set(atoms[atoms["shell"] == "PCS"]["atom_name"]), {"ND1"})

    def test_structures_without_hydrogens_are_unaffected(self):
        plain = analyze_structure(PLASTOCYANIN, "CU", shells=3)
        stripped = analyze_structure(PLASTOCYANIN, "CU", shells=3, exclude_hydrogens=True)
        for name in ("residues", "atoms", "links", "contacts"):
            self.assertTrue(plain[name].reset_index(drop=True).equals(stripped[name].reset_index(drop=True)), name)


if __name__ == "__main__":
    unittest.main()
