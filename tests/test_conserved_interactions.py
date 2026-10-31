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


class SulfurAcceptorTests(unittest.TestCase):
    def contact(self, donor, acceptor, distance, **kwargs):
        """One pair in the table: donor at the origin, acceptor ``distance`` away, far apart in sequence."""
        atoms = atoms_frame([("s", 1, 1, *donor, 0.0, 0.0, 0.0), ("s", 9, 9, *acceptor, distance, 0.0, 0.0), *kwargs.get("extra", [])])
        table = ci.residue_contacts(atoms, sulfur_acceptor=True)
        return table[(table.position_a == 1) & (table.position_b == 9)].iloc[0]

    def test_backbone_nh_to_cys_sulfur_follows_the_window(self):
        for distance, expected in ((3.5, True), (3.7, True), (3.75, False)):
            row = self.contact(("ALA", "N", "N"), ("CYS", "SG", "S"), distance)
            self.assertEqual((row["polar_s"], row["polar_s_bb"]), (expected, expected), distance)
            self.assertFalse(row["polar"])                                  # the N/O rule is untouched
        self.assertFalse(self.contact(("ALA", "N", "N"), ("CYS", "SG", "S"), 3.5)["hydrophobic"])   # an N...S pair is not a C/S pair

    def test_met_window_is_shorter_and_either_residue_order_works(self):
        self.assertTrue(self.contact(("ALA", "N", "N"), ("MET", "SD", "S"), 3.6)["polar_s"])
        self.assertFalse(self.contact(("ALA", "N", "N"), ("MET", "SD", "S"), 3.65)["polar_s"])
        self.assertTrue(self.contact(("CYS", "SG", "S"), ("ALA", "N", "N"), 3.5)["polar_s"])     # acceptor at the lower position

    def test_side_chain_donors_are_capable_but_not_backbone(self):
        row = self.contact(("SER", "OG", "O"), ("CYS", "SG", "S"), 3.4)
        self.assertTrue(row["polar_s"])
        self.assertFalse(row["polar_s_bb"])

    def test_non_donors_and_non_acceptors_are_not_flagged(self):
        self.assertFalse(self.contact(("ALA", "O", "O"), ("CYS", "SG", "S"), 3.4)["polar_s"])      # carbonyl O has no H
        self.assertFalse(self.contact(("PRO", "N", "N"), ("CYS", "SG", "S"), 3.4)["polar_s"])      # Pro N has no H
        self.assertFalse(self.contact(("ALA", "N", "N"), ("LEU", "CD1", "C"), 3.4)["polar_s"])
        self.assertFalse(self.contact(("ALA", "N", "N"), ("CYS", "CB", "C"), 3.4)["polar_s"])      # only the thiol/thiolate S accepts

    def test_side_chain_donors_have_a_shorter_window_than_backbone_n(self):
        cys = ("CYS", "SG", "S")
        for donor, distance, backbone_flag, flag in ((("ALA", "N", "N"), 3.65, True, True), (("SER", "OG", "O"), 3.65, False, False),
                                                       (("SER", "OG", "O"), 3.55, False, True)):
            row = self.contact(donor, cys, distance)
            self.assertEqual((row["polar_s_bb"], row["polar_s"]), (backbone_flag, flag), (donor, distance))

    def residue_chain(self, rows):
        """A donor residue (seq 2, with its preceding C at seq 1) and a Cys at seq 9, from (resname, seq, atom, element, xyz) rows."""
        return atoms_frame([("s", seq, seq, res, atom, el, *xyz) for res, seq, atom, el, xyz in rows])

    def test_geometric_flag_needs_the_hydrogen_to_point_at_the_sulfur(self):
        # Backbone amide placed from C(i-1), N and CA; one sulfur sits along the N-H direction, one across it.
        c_prev, n, ca = np.array([-1.33, 0.0, 0.0]), np.zeros(3), np.array([0.55, 1.3, 0.0])
        h = ci.protonation.backbone_amide_hydrogen(n, c_prev, ca)
        direction = (h - n) / np.linalg.norm(h - n)
        perpendicular = np.cross(direction, [0, 0, 1.0])
        for label, s_xyz, expected in (("along N-H", n + 3.4 * direction, True), ("across N-H", n + 3.4 * perpendicular, False)):
            atoms = self.residue_chain([("ALA", 1, "C", "C", c_prev), ("ALA", 2, "N", "N", n), ("ALA", 2, "CA", "C", ca), ("CYS", 9, "SG", "S", s_xyz)])
            row = ci.residue_contacts(atoms, sulfur_acceptor=True).set_index(["position_a", "position_b"]).loc[(2, 9)]
            self.assertTrue(row["polar_s_bb"], label)                                  # the capable flag only sees the distance
            self.assertEqual(row["polar_s_geom"], expected, label)

    def test_geometric_flag_skips_unplaceable_donors_and_missing_neighbors(self):
        s = ("CYS", 9, "SG", "S", np.array([3.4, 0.0, 0.0]))
        no_previous = self.residue_chain([("ALA", 2, "N", "N", np.zeros(3)), ("ALA", 2, "CA", "C", np.array([-0.5, 1.3, 0.0])), s])
        self.assertFalse(ci.residue_contacts(no_previous, sulfur_acceptor=True)["polar_s_geom"].any())      # no C(i-1) to place the amide H
        hydroxyl = self.residue_chain([("SER", 2, "OG", "O", np.zeros(3)), s])
        row = ci.residue_contacts(hydroxyl, sulfur_acceptor=True).iloc[0]
        self.assertTrue(row["polar_s"] and not row["polar_s_geom"])                      # the OH hydrogen rotates: capable, not testable

    def test_geometric_flag_for_a_planar_side_chain_and_the_metal_exclusion(self):
        # Gln amide in the xy plane: CD at the origin, NE2 at +x, OE1 off to the side; S placed along one N-H.
        cd, ne2, oe1 = np.zeros(3), np.array([1.33, 0.0, 0.0]), np.array([-0.6, 1.1, 0.0])
        h = ci.protonation.amide_nh2_hydrogens(ne2, cd, oe1)[0]
        s_xyz = ne2 + 3.4 * (h - ne2) / np.linalg.norm(h - ne2)
        atoms = self.residue_chain([("GLN", 2, "CD", "C", cd), ("GLN", 2, "NE2", "N", ne2), ("GLN", 2, "OE1", "O", oe1), ("CYS", 9, "SG", "S", s_xyz)])
        row = ci.residue_contacts(atoms, sulfur_acceptor=True).iloc[0]
        self.assertTrue(row["polar_s_geom"] and row["polar_s"] and not row["polar_s_bb"])
        metals = pd.DataFrame([("s", "ZN", *(ne2 + [0.0, 0.0, 2.0]))], columns=["structure_id", "element", "x", "y", "z"])
        row = ci.residue_contacts(atoms, sulfur_acceptor=True, metals=metals).iloc[0]
        self.assertFalse(row["polar_s_geom"] or row["polar_s"])                           # a metal-bound N has no H to donate

    def test_histidine_is_a_capable_donor_unless_a_metal_coordinates_it(self):
        his, cys = ("HIS", "ND1", "N"), ("CYS", "SG", "S")
        atoms = atoms_frame([("s", 1, 1, *his, 0.0, 0.0, 0.0), ("s", 9, 9, *cys, 3.4, 0.0, 0.0)])

        def metals(element, x, sid="s"):
            return pd.DataFrame([(sid, element, x, 0.0, 0.0)], columns=["structure_id", "element", "x", "y", "z"])

        def flagged(table):
            return bool(ci.residue_contacts(atoms, sulfur_acceptor=True, metals=table)["polar_s"].iloc[0])
        self.assertTrue(flagged(None))                                          # no metals given: nothing is excluded
        self.assertTrue(flagged(metals("CU", -3.0)))                            # a metal far from the N does not exclude it
        self.assertFalse(flagged(metals("CU", -2.0)))                           # Cu-N(His) 2.0 A: bond tier, the N is deprotonated
        self.assertTrue(flagged(metals("CU", -2.6)))                            # Cu-N 2.6 A: listed, beyond bond_max 2.4, not coordinated
        self.assertFalse(flagged(metals("ZN", -2.0)))                           # Zn is not in the table: the global 2.8 A applies
        self.assertTrue(flagged(metals("ZN", -3.0)))
        self.assertTrue(flagged(metals("CU", -2.0, sid="other")))               # a metal of another structure is ignored
        with_metal = ci.residue_contacts(atoms, sulfur_acceptor=True, metals=metals("CU", -2.0))
        self.assertFalse(with_metal["polar_s_bb"].iloc[0])

    def test_load_metals_reads_the_ions_of_a_structure(self):
        from pathlib import Path
        cif = Path(__file__).resolve().parents[1] / "reference_structures/0_Plastocyanin/1ag6.cif"
        table = ci.load_metals(cif, structure_id="1ag6")
        self.assertEqual(table["element"].tolist(), ["CU"])
        self.assertEqual(table["structure_id"].tolist(), ["1ag6"])

    def test_a_disulfide_sulfur_is_never_an_acceptor(self):
        partner = [("s", 30, 30, "CYS", "SG", "S", 3.4, 2.05, 0.0)]                                # 2.05 A from the acceptor SG
        row = self.contact(("ALA", "N", "N"), ("CYS", "SG", "S"), 3.4, extra=partner)
        self.assertFalse(row["polar_s"])

    def test_default_output_is_unchanged_and_the_option_needs_chemistry(self):
        atoms = atoms_frame([("s", 1, 1, "ALA", "N", "N", 0.0, 0.0, 0.0), ("s", 9, 9, "CYS", "SG", "S", 3.4, 0.0, 0.0)])
        plain, tiered = ci.residue_contacts(atoms), ci.residue_contacts(atoms, sulfur_acceptor=True)
        self.assertNotIn("polar_s", plain.columns)
        pd.testing.assert_frame_equal(tiered[list(plain.columns)], plain)
        with self.assertRaises(ValueError):
            ci.residue_contacts(atoms, chemistry=False, sulfur_acceptor=True)


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

    def test_load_atoms_keeps_the_first_conformer_when_there_is_no_a(self):
        import tempfile
        from pathlib import Path

        def line(serial, name, alt, num, x, occ, element):
            return (f"ATOM  {serial:5d} {name:<4s}{alt}SER A{num:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}{occ:6.2f}{10.0:6.2f}          {element:>2s}\n")

        pdb = "".join([line(1, "N", " ", 1, 0.0, 1.0, "N"), line(2, "OG", "B", 1, 2.0, 0.4, "O"), line(3, "OG", "C", 1, 9.0, 0.6, "O"), "END\n"])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bc.pdb"
            path.write_text(pdb)
            atoms = ci.load_atoms(path, "A")
        og = atoms[atoms["atom"] == "OG"]
        self.assertEqual(len(og), 1)
        self.assertAlmostEqual(float(og.iloc[0]["x"]), 2.0)

    def test_seq_index_counts_unobserved_and_hetero_residues_as_sequence_steps(self):
        import tempfile
        from pathlib import Path

        def line(serial, name, res, num, x, y=0.0, rec="ATOM  "):
            return f"{rec}{serial:5d} {name:<4s} {res:>3s} A{num:4d}    {x:8.3f}{y:8.3f}{0.0:8.3f}{1.0:6.2f}{10.0:6.2f}          {'C':>2s}\n"

        pdb = "".join([
            line(1, "CB", "LEU", 10, 0.0),
            line(2, "CB", "VAL", 15, 3.8),                       # 11-14 unobserved: 5 apart in sequence
            line(3, "CB", "LEU", 20, 0.0, 30.0),
            line(4, "CB", "MSE", 21, 0.0, 39.0, rec="HETATM"),   # skipped hetero residue between 20 and 22
            line(5, "CB", "VAL", 22, 3.8, 30.0),
            line(6, "CB", "ILE", 30, 0.0, 60.0),
            line(7, "CB", "ALA", 31, 3.8, 60.0),                 # consecutive: neighbors, excluded
            "END\n",
        ])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "gap.pdb"
            path.write_text(pdb)
            atoms = ci.load_atoms(path, "A")
        index = atoms.groupby("position")["seq_index"].first().to_dict()
        self.assertEqual(index[15] - index[10], 5)
        self.assertEqual(index[22] - index[20], 2)
        self.assertEqual(index[31] - index[30], 1)
        pairs = {(r.position_a, r.position_b) for r in ci.residue_contacts(atoms).itertuples()}
        self.assertEqual(pairs, {(10, 15), (20, 22)})


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

    def test_backbone_only_structures_are_flagged(self):
        side = ["CB", "CG", "CD1", "CD2"]
        rows = []
        for sid, atoms_per in (("full", ["N", "CA", "C", "O"] + side), ("trunc", ["N", "CA", "C", "O", "CB"])):
            for pos in range(1, 7):
                rows += [(sid, pos, pos, "LEU", a, a[0], float(pos), 0.0, 0.0) for a in atoms_per]
            rows += [(sid, 7, 7, "GLY", a, a[0], 7.0, 0.0, 0.0) for a in ("N", "CA", "C", "O")]   # Gly and Ala do not count
        rows += [("gly_only", 1, 1, "GLY", "CA", "C", 0.0, 0.0, 0.0)]
        audit = ci.backbone_only_audit(atoms_frame(rows)).set_index("structure_id")
        self.assertFalse(audit.loc["full", "backbone_only"])
        self.assertTrue(audit.loc["trunc", "backbone_only"])
        self.assertEqual(audit.loc["trunc", "n_residues"], 6)
        self.assertTrue(audit.loc["gly_only", "backbone_only"])      # nothing but Gly/Ala is mapped: no side-chain evidence
        half = atoms_frame([("h", p, p, "LEU", a, a[0], float(p), 0.0, 0.0) for p in (1, 2) for a in ["N", "CA", "C", "O", "CB", "CG"]]
                           + [("h", p, p, "LEU", a, a[0], float(p), 0.0, 0.0) for p in (3, 4) for a in ["N", "CA", "C", "O", "CB"]])
        self.assertFalse(ci.backbone_only_audit(half).iloc[0]["backbone_only"])   # exactly half is not more than half


if __name__ == "__main__":
    unittest.main()
