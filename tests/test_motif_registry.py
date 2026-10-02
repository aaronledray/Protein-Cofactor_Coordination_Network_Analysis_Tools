"""Validation for the canonical residue-specific motif vocabulary."""

import unittest

from modules.motif_registry import (
    MOTIF_ATOM_MEMBERSHIP,
    motif_atoms,
    motif_for_atom,
)


class MotifRegistryTests(unittest.TestCase):
    def test_histidine_protonation_aliases_share_imidazole_membership(self):
        for residue in ("HIS", "HID", "HIE", "HIP"):
            self.assertEqual(motif_for_atom(residue, "ND1"), "imidazole")
            self.assertEqual(motif_for_atom(residue, "NE2"), "imidazole")
            self.assertIn("ND1", motif_atoms(residue, "imidazole"))
            self.assertIn("NE2", motif_atoms(residue, "imidazole"))

    def test_aspartate_and_glutamate_carboxylates_are_residue_specific(self):
        self.assertEqual(motif_for_atom("ASP", "OD1"), "COO")
        self.assertEqual(motif_for_atom("ASP", "OD2"), "COO")
        self.assertEqual(motif_for_atom("GLU", "OE1"), "COO")
        self.assertEqual(motif_for_atom("GLU", "OE2"), "COO")
        self.assertEqual(motif_atoms("ASP", "carboxylate"), frozenset({"CG", "OD1", "OD2"}))
        self.assertEqual(motif_atoms("GLU", "COO"), frozenset({"CD", "OE1", "OE2"}))
        self.assertNotEqual(motif_for_atom("GLU", "CB"), "COO")

    def test_cysteine_thiol_is_not_the_entire_sidechain(self):
        self.assertEqual(motif_for_atom("CYS", "SG"), "thiol")
        self.assertIn("SG", motif_atoms("CYS", "thiol"))
        self.assertNotEqual(motif_for_atom("CYS", "CB"), "thiol")

    def test_heme_porphyrin_and_propionate_oxygen_are_distinct(self):
        self.assertEqual(motif_for_atom("HEM", "FE"), "heme_fe")
        self.assertEqual(motif_for_atom("HEM", "NA"), "heme_por")
        self.assertEqual(motif_for_atom("HEM", "O1A"), "heme_propionate")
        self.assertIn("NA", motif_atoms("HEM", "heme_porphyrin"))
        self.assertEqual(
            motif_atoms("HEM", "heme_propionate"),
            frozenset({"O1A", "O2A", "O1D", "O2D"}),
        )
        self.assertNotIn("O1A", motif_atoms("HEM", "heme_porphyrin"))

    def test_cluster_vocabularies_cover_reference_cofactor_atom_names(self):
        self.assertEqual(motif_for_atom("SF4", "FE3"), "iron_sulfur_fe")
        self.assertEqual(motif_for_atom("SF4", "S4"), "iron_sulfur_s")
        self.assertEqual(motif_for_atom("SF4", "FE5"), "unknown_motif")
        self.assertEqual(motif_for_atom("BCLF", "FE7"), "iron_sulfur_fe")
        self.assertEqual(motif_for_atom("BCLF", "S4A"), "iron_sulfur_s")
        self.assertEqual(motif_for_atom("ICS", "FE3"), "femoco_fe")
        self.assertEqual(motif_for_atom("ICS", "MO1"), "femoco_mo")
        self.assertEqual(motif_for_atom("ICS", "CX"), "femoco_carbide")
        self.assertEqual(motif_for_atom("HCA", "C7"), "hca_c")
        self.assertEqual(motif_for_atom("HCA", "O7"), "hca_o")
        self.assertEqual(motif_for_atom("OEX", "MN1"), "oxygen_evolving_cluster")
        self.assertIn("SF4", MOTIF_ATOM_MEMBERSHIP)
        self.assertIn("iron_sulfur_fe", MOTIF_ATOM_MEMBERSHIP["SF4"])


if __name__ == "__main__":
    unittest.main()
