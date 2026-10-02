"""Keeps the batch smoke benchmark runnable."""

import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("batch_smoke_benchmark", ROOT / "benchmarks/batch_smoke_benchmark.py")
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)


class BenchmarkSmokeTests(unittest.TestCase):
    def test_case_records_counts_and_isolates_failures(self):
        cases = {c["name"]: c for c in bench.CASES}
        ok = bench.run_case(cases["cu_plastocyanin"], shells=2)
        self.assertEqual(ok["status"], "ok")
        self.assertEqual(ok["n_sites"], 1)
        self.assertGreater(ok["n_contacts"], 0)
        failed = bench.run_case(cases["missing_input_isolated"], shells=2)
        self.assertEqual(failed["status"], "failed")
        self.assertIn("FileNotFoundError", failed["error"])

    def test_multiprocessing_matches_serial_and_isolates_errors(self):
        result = bench.multiprocessing_check(workers=2)
        self.assertTrue(result["serial_equals_parallel"])
        self.assertEqual(result["errors_isolated"], 1)


if __name__ == "__main__":
    unittest.main()
