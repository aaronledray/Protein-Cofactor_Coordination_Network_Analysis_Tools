"""Integration test for the additive headless SSCNA option."""

from pathlib import Path
from tempfile import TemporaryDirectory
import os
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SSCNA = ROOT / "1_Single_Structure_Cofactor_Network_Analysis_SSCNA_v0.0.2.py"
STRUCTURE = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"
FIXTURE_DIR = ROOT / "tests/fixtures/baseline/1ag6"


class NoPlotsCliTests(unittest.TestCase):
    def test_no_plots_writes_csvs_without_rendered_artifacts(self):
        with TemporaryDirectory(prefix="sscna-no-plots-") as work_dir:
            env = os.environ.copy()
            env["MPLCONFIGDIR"] = str(Path(work_dir) / "mpl-cache")
            env["XDG_CACHE_HOME"] = str(Path(work_dir) / "cache")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SSCNA),
                    "--template",
                    str(STRUCTURE),
                    "--cofactor",
                    "CU",
                    "--distance",
                    "3.6",
                    "--exclude-moieties",
                    "alanine_sidechain",
                    "--mode",
                    "Coord_Network",
                    "--no-plots",
                ],
                cwd=work_dir,
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )

            self.assertEqual(completed.stdout, "")

            output_dir = Path(work_dir) / "SSCNA_output"
            actual_breakdown = output_dir / "1ag6.cif_Coord_Breakdown.csv"
            actual_links = output_dir / "1ag6.cif_Coord_Links.csv"
            self.assertEqual(
                actual_breakdown.read_bytes(),
                (FIXTURE_DIR / "Coord_Breakdown.csv").read_bytes(),
            )
            self.assertEqual(
                actual_links.read_bytes(),
                (FIXTURE_DIR / "Coord_Links.csv").read_bytes(),
            )
            self.assertEqual(
                sorted(p.suffix for p in output_dir.iterdir()),
                [".csv", ".csv", ".csv"],
            )


if __name__ == "__main__":
    unittest.main()
