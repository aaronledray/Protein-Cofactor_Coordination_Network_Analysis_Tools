"""Tests for assembly/symmetry context and equivalent-site grouping."""

from pathlib import Path
import tempfile
import unittest

import pandas as pd

from modules.assembly_policy import ASSEMBLY_POLICY, equivalent_site_groups, symmetry_context
from modules.ml_export import analyze_labeling_example

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference_structures"

BIOMT_PDB = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
REMARK 350 BIOMOLECULE: 1
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000       10.00000
REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
END
"""


class AssemblyPolicyTests(unittest.TestCase):
    def test_identity_only_files_report_no_unapplied_transforms(self):
        context = symmetry_context(REF / "3_Metals/2zvs.pdb")
        self.assertEqual(context["assembly_policy"], ASSEMBLY_POLICY)
        self.assertEqual(context["n_assemblies"], 3)
        self.assertFalse(context["has_unapplied_transforms"])
        cif = symmetry_context(REF / "0_Plastocyanin/1ag6.cif")
        self.assertEqual(cif["space_group"], "P 31 2 1")
        self.assertFalse(cif["has_unapplied_transforms"])

    def test_non_identity_biomt_is_flagged(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "dimer.pdb"
            path.write_text(BIOMT_PDB)
            context = symmetry_context(path)
        self.assertEqual(context["non_identity_operators"], 1)
        self.assertTrue(context["has_unapplied_transforms"])
        self.assertEqual(context["space_group"], "P 1")

    def test_contract_records_symmetry_context(self):
        _, contract = analyze_labeling_example(REF / "3_Metals/1ca2.pdb", "ZN")
        self.assertEqual(contract["symmetry_context"]["assembly_policy"], ASSEMBLY_POLICY)

    def test_ncs_like_fes_sites_group_at_coordination_level_only(self):
        labels, _ = analyze_labeling_example(REF / "3_Metals/2zvs.pdb", "SF4")
        self.assertEqual(labels["site_id"].nunique(), 6)
        coarse = equivalent_site_groups(labels)
        self.assertEqual(coarse["equivalence_group"].nunique(), 1)
        self.assertTrue((coarse["group_size"] == 6).all())
        strict = equivalent_site_groups(labels, level="full")
        self.assertEqual(strict["equivalence_group"].nunique(), 6)

    def test_grouping_is_deterministic_and_ignores_numbering(self):
        labels, _ = analyze_labeling_example(REF / "3_Metals/2zvs.pdb", "SF4")
        first = equivalent_site_groups(labels, level="full")
        shuffled = labels.sample(frac=1, random_state=1).reset_index(drop=True)
        renumbered = shuffled.assign(residue_number=shuffled["residue_number"].astype(int) + 500,
                                     chain="Z")
        pd.testing.assert_frame_equal(first, equivalent_site_groups(renumbered, level="full"))

    def test_different_chemistry_is_not_grouped_and_bad_level_rejected(self):
        cu, _ = analyze_labeling_example(REF / "0_Plastocyanin/1ag6.cif", "CU")
        zn, _ = analyze_labeling_example(REF / "3_Metals/1ca2.pdb", "ZN")
        groups = equivalent_site_groups(pd.concat([cu, zn], ignore_index=True))
        self.assertEqual(groups["equivalence_group"].nunique(), 2)
        with self.assertRaises(ValueError):
            equivalent_site_groups(cu, level="bogus")
        self.assertTrue(equivalent_site_groups(cu.iloc[0:0]).empty)


if __name__ == "__main__":
    unittest.main()
