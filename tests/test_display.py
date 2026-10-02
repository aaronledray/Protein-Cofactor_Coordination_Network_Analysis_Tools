"""Tests for suppressing figure windows in non-interactive runs."""

from pathlib import Path
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from modules import display

ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


class FakePlt:
    def __init__(self, figures=1):
        self.shown, self.closed, self._figures = 0, 0, figures

    def show(self):
        self.shown += 1

    def get_fignums(self):
        return [1] * self._figures

    def gcf(self):
        return "current"

    def close(self, figure):
        self.closed += 1


class FakeFigure:
    shown = 0

    def show(self):
        FakeFigure.shown += 1


class DisplayDecisionTests(unittest.TestCase):
    def tearDown(self):
        display.set_show_figures(None)

    def test_auto_requires_both_streams_to_be_terminals(self):
        for stdin_tty, stdout_tty, expected in ((True, True, True), (True, False, False),
                                                (False, True, False), (False, False, False)):
            with mock.patch.object(display.sys, "stdin") as fake_in, \
                    mock.patch.object(display.sys, "stdout") as fake_out:
                fake_in.isatty.return_value, fake_out.isatty.return_value = stdin_tty, stdout_tty
                self.assertEqual(display.figures_should_show(), expected)

    def test_closed_streams_mean_no_show(self):
        with mock.patch.object(display.sys, "stdin") as fake_in:
            fake_in.isatty.side_effect = ValueError("closed")
            self.assertFalse(display.figures_should_show())

    def test_override_wins_over_auto(self):
        display.set_show_figures(True)
        self.assertTrue(display.figures_should_show())
        display.set_show_figures(False)
        self.assertFalse(display.figures_should_show())

    def test_helpers_show_or_release(self):
        plt, figure = FakePlt(), FakeFigure()
        display.set_show_figures(False)
        display.show_matplotlib(plt)
        display.show_plotly(figure)
        self.assertEqual((plt.shown, plt.closed, FakeFigure.shown), (0, 1, 0))
        display.show_matplotlib(FakePlt(figures=0))  # nothing open: no error
        display.set_show_figures(True)
        display.show_matplotlib(plt)
        display.show_plotly(figure)
        self.assertEqual((plt.shown, FakeFigure.shown), (1, 1))


class LegacyRunDisplayTests(unittest.TestCase):
    def _run(self, directory, *extra):
        marker = Path(directory) / "browser_called"
        launcher = Path(directory) / "browser.py"
        launcher.write_text(
            f"#!{sys.executable}\nimport pathlib, sys, urllib.request\n"
            f"pathlib.Path({str(marker)!r}).write_text(sys.argv[1])\n"
            "urllib.request.urlopen(sys.argv[1]).read()\n"
        )
        launcher.chmod(0o755)
        done = subprocess.run(
            [sys.executable, str(ROOT / "sscna_cli.py"), "analyze", "--input", str(PLASTOCYANIN),
             "--cofactor", "CU", "--mode", "Coord_Network", "--compact-html", *extra],
            cwd=directory, capture_output=True, text=True, timeout=240,
            env={**os.environ, "MPLBACKEND": "Agg", "BROWSER": f"{launcher} %s &",
                 "MPLCONFIGDIR": str(Path(directory) / "mpl")},
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        files = {p.name for p in (Path(directory) / "SSCNA_output").iterdir()}
        return marker.exists(), files

    def test_non_interactive_default_opens_nothing_but_writes_the_same_files(self):
        with tempfile.TemporaryDirectory() as directory:
            opened, files = self._run(directory)
        self.assertFalse(opened)
        self.assertIn("1ag6.cif_1_template_coordination_network.html", files)
        self.assertEqual(len(files), 6)

    def test_no_show_flag_and_show_override(self):
        with tempfile.TemporaryDirectory() as directory:
            opened, _ = self._run(directory, "--no-show")
        self.assertFalse(opened)
        with tempfile.TemporaryDirectory() as directory:
            opened, files = self._run(directory, "--show")
        self.assertTrue(opened)
        self.assertEqual(len(files), 6)


if __name__ == "__main__":
    unittest.main()
