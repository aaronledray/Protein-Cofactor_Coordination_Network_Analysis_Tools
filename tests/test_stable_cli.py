"""Regression tests for the installed/stable command wrapper."""

from pathlib import Path
import os
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


class StableCliTests(unittest.TestCase):
    @staticmethod
    def _pdb_atom(serial, record, name, residue, chain, number, x, y, z, element):
        return (
            f"{record:<6}{serial:5d} {name:>4s} {residue:>3s} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}\n"
        )

    def test_analyze_wrapper_preserves_legacy_no_plots_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            completed = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "sscna_cli.py"),
                    "analyze",
                    "--input",
                    str(PLASTOCYANIN),
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
                cwd=directory,
                env={
                    **os.environ,
                    "MPLCONFIGDIR": str(Path(directory) / "mpl-cache"),
                    "XDG_CACHE_HOME": str(Path(directory) / "cache"),
                },
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            output_dir = Path(directory) / "SSCNA_output"
            self.assertTrue((output_dir / "1ag6.cif_Coord_Breakdown.csv").is_file())
            self.assertTrue((output_dir / "1ag6.cif_Coord_Links.csv").is_file())
            self.assertFalse((output_dir / "1ag6_coordination_network.html").exists())

    def _run_analyze(self, directory, *extra):
        completed = subprocess.run(
            [
                sys.executable, str(ROOT / "sscna_cli.py"), "analyze",
                "--input", str(PLASTOCYANIN), "--cofactor", "CU", "--distance", "3.6",
                "--exclude-moieties", "alanine_sidechain", "--mode", "Coord_Network",
                "--compact-html", *extra,
            ],
            cwd=directory,
            env={**os.environ, "MPLBACKEND": "Agg",
                 "MPLCONFIGDIR": str(Path(directory) / "mpl-cache"),
                 "XDG_CACHE_HOME": str(Path(directory) / "cache")},
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return {path.name for path in (Path(directory) / "SSCNA_output").iterdir()}

    def test_default_run_keeps_legacy_output_file_set(self):
        with tempfile.TemporaryDirectory() as directory:
            names = self._run_analyze(directory)
        self.assertEqual(
            names,
            {
                "1ag6.cif_Coord_Breakdown.csv",
                "1ag6.cif_Coord_Breakdown_atoms.csv",
                "1ag6.cif_Coord_Links.csv",
                "1ag6.cif_1_static_element_coloring_mode_both.png",
                "1ag6.cif_1_static_pcs_scs_mode_both.png",
                "1ag6.cif_1_template_coordination_network.html",
            },
        )

    def test_cohesive_viewer_is_opt_in(self):
        with tempfile.TemporaryDirectory() as directory:
            names = self._run_analyze(directory, "--cohesive-viewer")
        self.assertIn("1ag6.cif_coordination_network.html", names)
        self.assertIn("1ag6.cif_Coord_Contacts.csv", names)
        self.assertIn("1ag6.cif_Coordination_Atoms.csv", names)
        self.assertNotIn("1ag6.cif_1_template_coordination_network.html", names)

    def test_compare_subcommand_dispatches_to_pairwise_runner(self):
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory) / "comparison"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "sscna_cli.py"),
                    "compare",
                    "pairwise",
                    "--reference",
                    str(PLASTOCYANIN),
                    "--query",
                    str(PLASTOCYANIN),
                    "--cofactor",
                    "CU",
                    "--output-dir",
                    str(output_dir),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue((output_dir / "pairwise_comparisons.csv").is_file())
            self.assertTrue((output_dir / "pairwise_similarity_heatmap.png").is_file())

    def test_wire_subcommand_writes_cofactor_to_target_tables(self):
        with tempfile.TemporaryDirectory() as directory:
            structure = Path(directory) / "wire_test.pdb"
            structure.write_text(
                "".join(
                    [
                        self._pdb_atom(1, "HETATM", "FE", "HEM", "A", 1500, 0.0, 0.0, 0.0, "FE"),
                        self._pdb_atom(2, "ATOM", "NE1", "TRP", "A", 107, 2.5, 0.0, 0.0, "N"),
                        self._pdb_atom(3, "ATOM", "OH", "TYR", "A", 200, 5.3, 0.0, 0.0, "O"),
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )
            output_dir = Path(directory) / "wire_output"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "sscna_cli.py"),
                    "wire",
                    "--input",
                    str(structure),
                    "--cofactor",
                    "HEM",
                    "--target",
                    "TYR:A:200:OH",
                    "--output-dir",
                    str(output_dir),
                ],
                cwd=directory,
                env={
                    **os.environ,
                    "MPLCONFIGDIR": str(Path(directory) / "mpl-cache"),
                    "XDG_CACHE_HOME": str(Path(directory) / "cache"),
                },
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue((output_dir / "wire_nodes.csv").is_file())
            self.assertTrue((output_dir / "wire_edges.csv").is_file())
            self.assertTrue((output_dir / "wire_paths.csv").is_file())
            self.assertTrue((output_dir / "wire_targets.csv").is_file())
            self.assertTrue((output_dir / "wire_summary.json").is_file())
            self.assertTrue((output_dir / "wire_network.html").is_file())


if __name__ == "__main__":
    unittest.main()
