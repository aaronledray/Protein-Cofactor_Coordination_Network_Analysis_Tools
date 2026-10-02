"""Tests for the substrate-point viewer and `sscna substrate` command."""

import json
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import unittest

from modules.substrate_seeds import analyze_substrate_seed
from modules.substrate_viewer import VISIBLE_CHAINS_BY_DEFAULT, build_substrate_figure

ROOT = Path(__file__).resolve().parents[1]
MYOGLOBIN = ROOT / "reference_structures/3_Myoglobin/1a6m.pdb"
OFFSET = ("1.775", "-0.841", "0.378")


class SubstrateViewerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = analyze_substrate_seed(
            MYOGLOBIN, "HEM", {"cofactor_atom": "FE", "offset": tuple(map(float, OFFSET))}
        )

    def test_figure_has_shells_point_contacts_and_ranked_chains(self):
        figure = build_substrate_figure(self.result, pdb_name="1a6m")
        names = [trace.name for trace in figure.data]
        for expected in ("Cofactor", "PCS", "SCS", "TCS", "point contacts", "substrate point"):
            self.assertIn(expected, names)
        chain_traces = [t for t in figure.data if str(t.name).startswith("#")]
        self.assertEqual(len(chain_traces), len(self.result["chains"]))
        visible = [t.visible for t in chain_traces]
        self.assertEqual(visible[:VISIBLE_CHAINS_BY_DEFAULT], [True] * min(VISIBLE_CHAINS_BY_DEFAULT, len(visible)))
        # the distal-histidine chain ends on a cofactor atom and starts at the point
        distal = next(t for t in chain_traces if "HIS64" in t.name)
        self.assertEqual(distal.text[0], "substrate point")
        self.assertTrue(distal.text[-1].startswith("HEM154:"))
        self.assertIn("structural hypothesis", figure.layout.title.text)

    def test_empty_seed_is_rejected(self):
        empty = dict(self.result)
        empty["seed"] = self.result["seed"].iloc[0:0]
        with self.assertRaises(ValueError):
            build_substrate_figure(empty)


class SubstrateCliTests(unittest.TestCase):
    def _run(self, directory, *extra):
        return subprocess.run(
            [sys.executable, str(ROOT / "sscna_cli.py"), "substrate", "--input", str(MYOGLOBIN),
             "--cofactor", "HEM", "--output-dir", str(Path(directory) / "out"), *extra],
            cwd=directory,
            env={**os.environ, "MPLBACKEND": "Agg", "MPLCONFIGDIR": str(Path(directory) / "mpl")},
            capture_output=True, text=True, check=False,
        )

    def test_substrate_command_writes_tables_summary_and_viewer(self):
        with tempfile.TemporaryDirectory() as directory:
            done = self._run(directory, "--from-atom", "FE", "--offset", *OFFSET, "--compact-html")
            self.assertEqual(done.returncode, 0, done.stderr)
            out = Path(directory) / "out"
            for name in ("substrate_seed.csv", "substrate_contacts.csv", "substrate_chains.csv",
                         "substrate_steps.csv", "substrate_summary.json", "substrate_chains.html"):
                self.assertTrue((out / name).is_file(), name)
            summary = json.loads((out / "substrate_summary.json").read_text())
            self.assertEqual(summary["chain_count"], 3)
            self.assertTrue(summary["html_written"])
            self.assertIn("hypotheses", summary["note"])

    def test_explicit_point_no_html_and_argument_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            done = self._run(directory, "--point", "500", "500", "500", "--no-html")
            self.assertEqual(done.returncode, 0, done.stderr)
            summary = json.loads((Path(directory) / "out/substrate_summary.json").read_text())
            self.assertEqual((summary["chain_count"], summary["html_written"]), (0, False))
            self.assertNotEqual(self._run(directory, "--from-atom", "FE").returncode, 0)
            self.assertNotEqual(self._run(directory, "--offset", "0", "0", "1").returncode, 0)
            self.assertNotEqual(self._run(directory).returncode, 0)


if __name__ == "__main__":
    unittest.main()
