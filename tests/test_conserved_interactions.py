"""Tests for the conserved-interaction analysis (synthetic geometry, no structure files)."""

import unittest

import numpy as np
import pandas as pd

from modules import conserved_interactions as ci


def atoms_frame(rows):
    return pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])


def meta(rows):
    return pd.DataFrame(rows, columns=["structure_id", "group", "sequence_id"])


class ContactTests(unittest.TestCase):
    def test_hydrophobic_contact_and_separation_rules(self):
        atoms = atoms_frame([
            ("s", 1, 1, "LEU", "CB", "C", 0.0, 0.0, 0.0),
            ("s", 5, 5, "VAL", "CB", "C", 3.6, 0.0, 0.0),
            ("s", 2, 2, "ALA", "CB", "C", 0.0, 3.0, 0.0),     # adjacent to position 1: excluded
            ("s", 9, 9, "ILE", "CB", "C", 0.0, 0.0, 1.5),     # closer than 1.9 A to position 1: covalent candidate
        ])
        c = ci.residue_contacts(atoms)
        pairs = {(r.position_a, r.position_b) for r in c.itertuples()}
        self.assertIn((1, 5), pairs)
        self.assertNotIn((1, 2), pairs)
        self.assertFalse(any(9 in p and 1 in p for p in pairs))
        row = c[(c.position_a == 1) & (c.position_b == 5)].iloc[0]
        self.assertEqual(row["closest"], "sc_sc")
        self.assertTrue(row["hydrophobic"])
        self.assertFalse(row["polar"])

    def test_polar_ionic_and_backbone_flags(self):
        atoms = atoms_frame([
            ("s", 1, 1, "ASP", "OD1", "O", 0.0, 0.0, 0.0),
            ("s", 6, 6, "ARG", "NH1", "N", 2.9, 0.0, 0.0),
            ("s", 20, 20, "GLY", "N", "N", 20.0, 0.0, 0.0),
            ("s", 24, 24, "GLY", "O", "O", 22.9, 0.0, 0.0),
        ])
        c = ci.residue_contacts(atoms).set_index(["position_a", "position_b"])
        side = c.loc[(1, 6)]
        self.assertTrue(side["polar"] and side["polar_sc"] and side["ionic"])
        self.assertFalse(side["polar_bb_bb"])
        back = c.loc[(20, 24)]
        self.assertTrue(back["polar_bb_bb"])
        self.assertEqual(back["closest"], "bb_bb")

    def test_aromatic_ring_atom_contact(self):
        atoms = atoms_frame([
            ("s", 1, 1, "HIS", "NE2", "N", 0.0, 0.0, 0.0),
            ("s", 8, 8, "GLY", "CA", "C", 3.7, 0.0, 0.0),
        ])
        row = ci.residue_contacts(atoms).iloc[0]
        self.assertTrue(row["aromatic"])


class FeatureTests(unittest.TestCase):
    def setUp(self):
        self.atoms = atoms_frame([
            ("s", 1, 1, "HIS", "NE2", "N", 2.2, 0.0, 0.0),
            ("s", 5, 5, "ARG", "NH1", "N", 0.0, 3.0, 0.0),
            ("s", 9, 9, "ASP", "OD1", "O", 0.0, 6.0, 0.0),
            ("s", 12, 12, "LEU", "CD1", "C", 30.0, 0.0, 0.0),
        ])
        self.features = pd.DataFrame([
            ("s", "CuB", "CU", "CU", "metal", 0.0, 0.0, 0.0),
            ("s", "prop", "O1A", "O", "oxygen", 0.0, 0.0, 0.0),
        ], columns=["structure_id", "feature", "atom", "element", "kind", "x", "y", "z"])

    def test_classes_and_absent_residues(self):
        c = ci.cofactor_contacts(self.atoms, self.features).set_index(["feature", "position"])
        self.assertTrue(c.loc[("CuB", 1), "coord"])
        self.assertFalse(c.loc[("prop", 1), "coord"])
        self.assertTrue(c.loc[("prop", 5), "polar"] and c.loc[("prop", 5), "ionic"])
        self.assertNotIn(12, [p for _, p in c.index])            # farther than 8 A: not recorded

    def test_water_bridge_without_direct_contact(self):
        waters = pd.DataFrame([("s", 0.0, 3.0, 2.5)], columns=["structure_id", "x", "y", "z"])
        atoms = atoms_frame([
            ("s", 9, 9, "ASP", "OD1", "O", 0.0, 5.0, 2.5),
            ("s", 5, 5, "ARG", "NH1", "N", 3.0, 0.0, 0.0),
        ])
        waters.loc[0, ["y", "z"]] = [2.0, 1.5]    # 2.5 A from the feature O, 3.2 A from the Asp O
        b = ci.water_bridges(atoms, self.features, waters).set_index(["feature", "position"])
        self.assertTrue(b.loc[("prop", 9), "bridged"])
        self.assertFalse(b.loc[("prop", 9), "direct"])
        self.assertTrue(b.loc[("prop", 9), "indirect"])
        self.assertTrue(b.loc[("prop", 5), "direct"])

    def test_relay_steps_through_side_chains_and_waters(self):
        atoms = atoms_frame([
            ("s", 5, 5, "ARG", "NH1", "N", 0.0, 3.0, 0.0),
            ("s", 5, 5, "ARG", "NH2", "N", 0.0, 3.5, 2.0),
            ("s", 9, 9, "ASP", "OD1", "O", 0.0, 8.0, 3.2),
            ("s", 11, 11, "GLY", "O", "O", 0.0, 9.5, 3.2),       # backbone only: terminal
        ])
        waters = pd.DataFrame([("s", 0.0, 5.0, 3.0)], columns=["structure_id", "x", "y", "z"])
        sources = pd.DataFrame([("s", 0.0, 0.0, 0.0)], columns=["structure_id", "x", "y", "z"])
        steps = ci.relay_steps(atoms, sources, waters, link=3.5, max_steps=6).set_index("position")["steps"]
        self.assertEqual(steps[5], 1)
        self.assertEqual(steps[9], 3)
        dry = ci.relay_steps(atoms, sources, None, link=3.5, max_steps=6).set_index("position")["steps"]
        self.assertNotIn(9, dry.index)


class VoteTests(unittest.TestCase):
    def test_one_vote_per_sequence_and_group_support(self):
        # group A: three sequences, the first with three replicate structures (2 of 3 positive: a positive vote)
        structures = meta([("a1x", "A", "a1"), ("a1y", "A", "a1"), ("a1z", "A", "a1"), ("a2", "A", "a2"), ("a3", "A", "a3"),
                           ("b1", "B", "b1"), ("b2", "B", "b2"), ("b3", "B", "b3")])
        obs = pd.DataFrame({"structure_id": structures.structure_id, "key": "k",
                            "flag": [True, True, False, True, True, True, True, False]})
        v = ci.sequence_vote(obs, structures, ["key"], ["flag"], min_sequences=3, support=0.9)
        row = v.iloc[0]
        self.assertEqual((row["flag_n_A"], row["flag_A"]), (3, 1.0))
        self.assertAlmostEqual(row["flag_B"], 2 / 3, places=2)
        self.assertFalse(row["conserved_flag"])
        v2 = ci.sequence_vote(obs, structures, ["key"], ["flag"], min_sequences=3, support=0.6)
        self.assertTrue(v2.iloc[0]["conserved_flag"])
        v3 = ci.sequence_vote(obs, structures, ["key"], ["flag"], min_sequences=4, support=0.5)
        self.assertFalse(v3.iloc[0]["conserved_flag"])             # too few sequences per group

    def test_conserved_contacts_end_to_end(self):
        rows, struct = [], []
        for g, n in (("A", 4), ("B", 4)):
            for k in range(n):
                sid = f"{g}{k}"
                struct.append((sid, g, f"seq_{sid}"))
                # pair (1,5): always in contact, polar; pair (1,8): in contact in all but one B structure
                rows += [(sid, 1, 1, "ASP", "OD1", "O", 0.0, 0.0, 0.0), (sid, 5, 5, "ARG", "NH1", "N", 2.9, 0.0, 0.0),
                         (sid, 8, 8, "LEU", "CD1", "C", 0.0, 3.8 if not (g == "B" and k == 0) else 9.0, 0.0)]
        atoms, structures = atoms_frame(rows), meta(struct)
        v = ci.conserved_contacts(atoms, structures, min_sequences=4, support=0.9).set_index(["position_a", "position_b"])
        self.assertTrue(v.loc[(1, 5), "conserved_contact"])
        self.assertTrue(v.loc[(1, 5), "conserved_polar"] and v.loc[(1, 5), "conserved_ionic"])
        self.assertFalse(v.loc[(1, 8), "conserved_contact"])
        summary = ci.chemistry_summary(v.reset_index())
        self.assertEqual(int(summary.set_index("class").loc["polar", "all"]), 1)


class LoadTests(unittest.TestCase):
    def test_load_atoms_and_position_map(self):
        from pathlib import Path

        cif = Path(__file__).resolve().parents[1] / "reference_structures/0_Plastocyanin/1ag6.cif"
        atoms = ci.load_atoms(cif, "A")
        self.assertGreater(len(atoms), 500)
        self.assertEqual(set(atoms.columns), {"structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"})
        self.assertFalse(atoms["element"].isin(["H", "D"]).any())
        mapped = ci.load_atoms(cif, "A", {(int(n), ""): n + 1000 for n in atoms["position"].unique()}, structure_id="p")
        self.assertEqual(set(mapped["structure_id"]), {"p"})
        self.assertTrue((mapped["position"] >= 1000).all())
        contacts = ci.residue_contacts(atoms)
        self.assertGreater(len(contacts), 100)

    def test_load_atoms_keeps_the_a_conformer_when_b_has_higher_occupancy(self):
        import tempfile
        from pathlib import Path

        def line(serial, name, alt, res, num, x, occ, element):
            return (f"ATOM  {serial:5d} {name:<4s}{alt}{res:>3s} A{num:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}{occ:6.2f}{10.0:6.2f}          {element:>2s}\n")

        pdb = "".join([
            line(1, "N", " ", "SER", 1, 0.0, 1.0, "N"),
            line(2, "CA", " ", "SER", 1, 1.5, 1.0, "C"),
            line(3, "OG", "A", "SER", 1, 2.0, 0.30, "O"),
            line(4, "OG", "B", "SER", 1, 9.0, 0.70, "O"),
            line(5, "N", " ", "GLY", 2, 3.0, 1.0, "N"),
            line(6, "CA", " ", "GLY", 2, 4.5, 1.0, "C"),
            "END\n",
        ])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "alt.pdb"
            path.write_text(pdb)
            atoms = ci.load_atoms(path, "A")
        og = atoms[(atoms["position"] == 1) & (atoms["atom"] == "OG")]
        self.assertEqual(len(og), 1)
        self.assertAlmostEqual(float(og.iloc[0]["x"]), 2.0)
        self.assertEqual(len(atoms), 6 - 1)


class AuditTests(unittest.TestCase):
    def test_misplaced_anchor_is_flagged_by_sequence(self):
        atoms = atoms_frame([
            ("good", 411, 411, "HIS", "NE2", "N", 2.0, 0.0, 0.0),
            ("bad", 411, 411, "THR", "OG1", "O", 40.0, 0.0, 0.0),
            ("gone", 7, 7, "GLY", "CA", "C", 1.0, 0.0, 0.0),
        ])
        sites = pd.DataFrame([(s, 0.0, 0.0, 0.0) for s in ("good", "bad", "gone")], columns=["structure_id", "x", "y", "z"])
        structures = meta([("good", "A", "s1"), ("bad", "A", "s2"), ("gone", "A", "s3")])
        audit = ci.anchor_audit(atoms, sites, 411, structures, max_distance=8.0).set_index("structure_id")
        self.assertTrue(audit.loc["good", "ok"])
        self.assertFalse(audit.loc["bad", "ok"])
        self.assertFalse(audit.loc["gone", "ok"])
        self.assertEqual(ci.misaligned_sequences(audit.reset_index()), ["s2", "s3"])


if __name__ == "__main__":
    unittest.main()
