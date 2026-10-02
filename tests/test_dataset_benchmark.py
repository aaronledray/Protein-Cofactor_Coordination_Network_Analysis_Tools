"""Keeps the dataset-scale benchmark helpers working."""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
# A real import (not spec loading) so process-pool workers can unpickle run_one.
sys.path.insert(0, str(ROOT / "benchmarks"))
import dataset_benchmark as bench  # noqa: E402
REF = ROOT / "reference_structures"


class DatasetBenchmarkTests(unittest.TestCase):
    def test_sampling_is_even_deterministic_and_bounded(self):
        paths = [Path(f"s{i:03d}.pdb") for i in range(100)]
        sample = bench.sample_paths(reversed(paths), 10)
        self.assertEqual(sample, bench.sample_paths(paths, 10))
        self.assertEqual(len(sample), 10)
        self.assertEqual(sample[0], paths[0])
        self.assertEqual(len(bench.sample_paths(paths, 0)), 100)
        self.assertEqual(len(bench.sample_paths(paths, 500)), 100)

    def test_run_one_classifies_ok_missing_cofactor_and_failure(self):
        ok = bench.run_one((str(REF / "0_Plastocyanin/1ag6.cif"), "CU"))
        self.assertEqual((ok["status"], ok["n_sites"]), ("ok", 1))
        none = bench.run_one((str(REF / "0_Plastocyanin/1ag6.cif"), "ZZZ"))
        self.assertEqual(none["status"], "no_cofactor")
        bad = bench.run_one((str(REF / "missing.pdb"), "CU"))
        self.assertEqual(bad["status"], "failed")
        self.assertIn("FileNotFoundError", bad["error"])

    def test_pass_and_summary_with_two_workers(self):
        paths = [REF / "0_Plastocyanin/1ag6.cif", REF / "3_Metals/1ca2.pdb"]
        frame, wall = bench.run_pass(paths, "CU", workers=2)
        summary = bench.summarize(frame, wall, 2)
        self.assertEqual(summary["structures"], 2)
        self.assertEqual(summary["status_counts"], {"ok": 1, "no_cofactor": 1})


if __name__ == "__main__":
    unittest.main()
