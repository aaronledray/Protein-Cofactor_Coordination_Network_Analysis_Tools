"""Tests for profile-driven conserved-residue viewer output."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from modules.comparison_runner import profile_reference_set
from modules.conservation_viewer import write_profile_conservation_viewer


ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


class ConservationViewerTests(unittest.TestCase):
    def test_profile_viewer_adds_family_conservation_color_mode(self):
        profile_result = profile_reference_set(
            [PLASTOCYANIN],
            "CU",
            query_inputs=[PLASTOCYANIN],
            analysis_options={
                "exclude_moieties": ["alanine_sidechain"],
                "shells": 3,
            },
        )

        with TemporaryDirectory(prefix="sscna-conservation-viewer-") as work_dir:
            output = write_profile_conservation_viewer(
                PLASTOCYANIN,
                "CU",
                profile_result,
                Path(work_dir) / "conservation.html",
                analysis_options={
                    "exclude_moieties": ["alanine_sidechain"],
                    "shells": 3,
                },
                compact_html=True,
            )
            html = output.read_text(encoding="utf-8")

        self.assertIn("Family conservation", html)
        self.assertIn("Family conservation: 100%", html)


if __name__ == "__main__":
    unittest.main()
