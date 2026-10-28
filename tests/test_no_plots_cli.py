"""Integration test for the headless `sscna analyze --no-plots` option."""

from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import sys
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STRUCTURE = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


class NoPlotsCliTests(unittest.TestCase):
    def test_no_plots_writes_tables_without_rendered_artifacts(self):
        with TemporaryDirectory(prefix="sscna-no-plots-") as work_dir:
            completed = subprocess.run(
                [sys.executable, str(ROOT / "sscna_cli.py"), "analyze", "--input", str(STRUCTURE),
                 "--cofactor", "CU", "--distance", "3.6", "--no-plots"],
                cwd=work_dir, check=True, capture_output=True, text=True,
            )
            self.assertEqual(completed.stdout, "")
            output_dir = Path(work_dir) / "SSCNA_output"
            self.assertEqual(
                sorted(p.name for p in output_dir.iterdir()),
                ["1ag6.cif_Coord_Contacts.csv", "1ag6.cif_Coord_Links.csv",
                 "1ag6.cif_Coordination_Atoms.csv", "1ag6.cif_Coordination_Residues.csv"],
            )
            residues = pd.read_csv(output_dir / "1ag6.cif_Coordination_Residues.csv")
            self.assertGreater(len(residues), 0)
            self.assertIn("shell", residues.columns)


if __name__ == "__main__":
    unittest.main()
