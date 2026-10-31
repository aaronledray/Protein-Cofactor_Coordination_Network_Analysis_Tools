"""Offline checks for scripts/plastocyanin_vote_check.py (alignment to the reference numbering) and the opt-in
sulfur flags in the conserved-contact vote."""

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from modules import conserved_interactions as ci

ROOT = Path(__file__).resolve().parents[1]
CIF = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("plastocyanin_vote_check", ROOT / "scripts/plastocyanin_vote_check.py")
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


def residues(sequence, start=1):
    return [((start + n, ""), letter) for n, letter in enumerate(sequence)]


class AlignmentTests(unittest.TestCase):
    def test_a_chain_aligned_to_itself_maps_every_residue_to_its_own_number(self):
        reference = check.chain_residues(CIF, "A")
        self.assertEqual(len(reference), 99)
        mapping = check.reference_map(reference, reference)
        self.assertEqual(mapping, {key: key[0] for key, _ in reference})
        self.assertEqual(reference[83][1], "C")                                        # Cys84

    def test_numbering_offset_and_a_gap_are_absorbed_into_the_reference_numbering(self):
        reference = residues("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV")
        shifted = residues("".join(a for _, a in reference), start=101)
        self.assertEqual(check.reference_map(reference, shifted)[(101 + 10, "")], 11)
        deleted = [r for n, r in enumerate(reference) if n not in (20, 21, 22)]            # a three-residue deletion in the other sequence
        mapping = check.reference_map(reference, deleted)
        self.assertNotIn(21, mapping.values())                                             # the deleted positions have no residue
        self.assertEqual(mapping[deleted[25][0]], deleted[25][0][0])                      # residues after the gap keep the reference number
        inserted = reference[:30] + residues("WWW", start=500) + reference[30:]
        mapping = check.reference_map(reference, inserted)
        self.assertEqual([key for key in mapping if key[0] >= 500], [])                   # inserted residues map to nothing


class VoteSulfurTests(unittest.TestCase):
    def test_conserved_contacts_vote_the_sulfur_flags_when_asked(self):
        rows, structures = [], []
        for n in range(6):
            sid = f"s{n}"
            rows += [(sid, 1, 1, "ASN", "N", "N", 0.0, 0.0, 0.0), (sid, 9, 9, "CYS", "SG", "S", 3.4, 0.0, 0.0)]
            structures.append((sid, "g", f"q{n}"))
        atom_table = pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])
        meta = pd.DataFrame(structures, columns=["structure_id", "group", "sequence_id"])
        plain = ci.conserved_contacts(atom_table, meta)
        self.assertNotIn("polar_s_g", plain.columns)
        voted = ci.conserved_contacts(atom_table, meta, sulfur_acceptor=True)
        pd.testing.assert_frame_equal(voted[list(plain.columns)], plain)                  # the default columns do not move
        row = voted.iloc[0]
        self.assertTrue(row["conserved_polar_s"] and row["conserved_polar_s_bb"])
        self.assertFalse(row["conserved_polar"])                                           # polar needs an N/O pair
        self.assertIn("polar_s_geom", set(ci.rules_of(voted)["name"]))                     # and the legend lists the sulfur rules as computed
        self.assertNotIn("was off", ci.rules_of(voted).set_index("name").loc["polar_s", "rule"])
        self.assertIn("was off", ci.rules_of(plain).set_index("name").loc["polar_s", "rule"])


if __name__ == "__main__":
    unittest.main()
