"""Regression tests for the pre-refactor coordination-network outputs."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from modules.io_utils import unpack_pdb_file
from modules.moieties import bond_lookup
from modules.structure_processing import (
    generate_coordination_csv_with_moieties,
    identify_coordination_network,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures" / "baseline"


CASES = (
    {
        "name": "1ag6",
        "structure": ROOT / "reference_structures/0_Plastocyanin/1ag6.cif",
        "cofactor": ["CU"],
        "cofactor2": None,
        "combinatorial": False,
        "counts": (1, 5, 2),
    },
    {
        "name": "4ub6",
        "structure": ROOT / "reference_structures/1_OEX/0_Aligned_Reduced/4ub6.pdb",
        "cofactor": ["OEX"],
        "cofactor2": None,
        "combinatorial": False,
        "counts": (10, 21, 16),
    },
    {
        "name": "3u7q",
        "structure": ROOT / "reference_structures/2_Nitrogenase/3u7q_monomer.pdb",
        "cofactor": ["ICS", "CLF"],
        "cofactor2": ["HCA"],
        "combinatorial": True,
        "counts": (47, 38, 62),
    },
    {
        "name": "1a6m",
        "structure": ROOT / "reference_structures/3_Myoglobin/1a6m.pdb",
        "cofactor": ["HEM"],
        "cofactor2": None,
        "combinatorial": False,
        "counts": (43, 17, 15),
    },
)


class BaselineCoordinationOutputTests(unittest.TestCase):
    """Ensure refactors preserve the current CSV output contract."""

    def test_reference_outputs_are_unchanged(self):
        for case in CASES:
            with self.subTest(structure=case["name"]):
                structure, _ = unpack_pdb_file(str(case["structure"]))

                with TemporaryDirectory(prefix=f"sscna-{case['name']}-") as output_dir:
                    with redirect_stdout(StringIO()):
                        cofactor, pcs, scs = identify_coordination_network(
                            structure=structure,
                            cofactor_resname=case["cofactor"],
                            distance_cutoff=3.6,
                            expand_residues=False,
                            combinatorial_mode=case["combinatorial"],
                            combinatorial_cofactor_cutoff=20.0,
                            cofactor_resname2=case["cofactor2"],
                            exclude_moieties=["alanine_sidechain"],
                            output_dir=output_dir,
                            output_prefix="",
                        )
                        generate_coordination_csv_with_moieties(
                            cofactor_sphere=cofactor,
                            pcs_residues=pcs,
                            scs_residues=scs,
                            pdb_file_name=case["structure"].name,
                            structure=structure,
                            bond_lookup=bond_lookup,
                            output_dir=output_dir,
                        )

                    self.assertEqual(
                        (len(cofactor), len(pcs), len(scs)), case["counts"]
                    )

                    actual_breakdown = (
                        Path(output_dir) / f"{case['structure'].name}_Coord_Breakdown.csv"
                    )
                    actual_links = Path(output_dir) / "Coord_Links.csv"
                    expected_dir = FIXTURES / case["name"]

                    self.assertEqual(
                        actual_breakdown.read_bytes(),
                        (expected_dir / "Coord_Breakdown.csv").read_bytes(),
                    )
                    self.assertEqual(
                        actual_links.read_bytes(),
                        (expected_dir / "Coord_Links.csv").read_bytes(),
                    )


if __name__ == "__main__":
    unittest.main()
