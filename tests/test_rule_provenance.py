"""The rules reported with every table must be the rules the code applies: each documented cutoff is tested by
putting a pair just inside it (flag set) and just outside it (flag clear)."""

import unittest

import numpy as np
import pandas as pd

from modules import conserved_interactions as ci
from modules import protonation

EPS = 0.01


def atoms(*rows):
    return pd.DataFrame([("s", p, p, r, a, e, *xyz) for p, r, a, e, xyz in rows],
                        columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])


def pair(first, second, distance):
    """Two residues far apart in sequence, ``first`` at the origin and ``second`` ``distance`` away along x."""
    (r1, a1, e1), (r2, a2, e2) = first, second
    return atoms((1, r1, a1, e1, (0.0, 0.0, 0.0)), (9, r2, a2, e2, (distance, 0.0, 0.0)))


class ResidueContactRuleTests(unittest.TestCase):
    def flags(self, first, second, distance, **kwargs):
        table = ci.residue_contacts(pair(first, second, distance), cutoff=5.0, **kwargs)
        return None if table.empty else table.iloc[0]

    def cutoff(self, name, **parameters):
        rules = ci.class_rules("residue_contacts", cutoff=5.0, **parameters).set_index("name")
        return float(rules.loc[name, "cutoff_A"])

    def check(self, name, first, second, **kwargs):
        limit = self.cutoff(name, **kwargs)
        inside, outside = self.flags(first, second, limit - EPS, **kwargs), self.flags(first, second, limit + EPS, **kwargs)
        self.assertTrue(inside[name], f"{name} should hold just inside {limit}")
        self.assertFalse(outside[name], f"{name} should not hold just outside {limit}")

    def test_each_default_flag_changes_at_its_documented_cutoff(self):
        self.check("polar", ("ASP", "OD1", "O"), ("ASN", "ND2", "N"))
        self.check("polar_sc", ("ASP", "OD1", "O"), ("ASN", "ND2", "N"))
        self.check("polar_bb_bb", ("GLY", "N", "N"), ("GLY", "O", "O"))
        self.check("hydrophobic", ("LEU", "CB", "C"), ("VAL", "CB", "C"))
        self.check("aromatic", ("PHE", "CG", "C"), ("ALA", "CB", "C"))
        self.check("ionic", ("ASP", "OD1", "O"), ("ARG", "NH1", "N"))

    def test_the_contact_window_is_the_documented_one(self):
        rules = ci.class_rules("residue_contacts", cutoff=4.5, min_distance=2.0).set_index("name")
        self.assertEqual(float(rules.loc["contact", "cutoff_A"]), 4.5)
        self.assertIn("between 2 and 4.5 A", rules.loc["contact", "rule"])
        table = ci.residue_contacts(pair(("LEU", "CB", "C"), ("VAL", "CB", "C"), 4.4), cutoff=4.5)
        self.assertEqual(len(table), 1)
        self.assertEqual(len(ci.residue_contacts(pair(("LEU", "CB", "C"), ("VAL", "CB", "C"), 4.6), cutoff=4.5)), 0)

    def test_a_cutoff_below_a_class_window_caps_the_documented_value(self):
        rules = ci.class_rules("residue_contacts", cutoff=3.2).set_index("name")
        self.assertEqual(float(rules.loc["polar", "cutoff_A"]), 3.2)                 # min(polar window, cutoff)
        self.assertEqual(float(rules.loc["hydrophobic", "cutoff_A"]), 3.2)

    def test_sulfur_rules_change_at_their_documented_windows(self):
        backbone, side = ("ALA", "N", "N"), ("SER", "OG", "O")
        limits = {"polar_s_bb": ci.SULFUR_ACCEPTORS[("CYS", "SG")], "polar_s": ci.SULFUR_SIDE_CHAIN_MAX}
        rules = ci.class_rules("residue_contacts", cutoff=5.0, sulfur_acceptor=True).set_index("name")
        self.assertEqual(float(rules.loc["polar_s_bb", "cutoff_A"]), limits["polar_s_bb"])
        self.assertEqual(float(rules.loc["polar_s", "cutoff_A"]), limits["polar_s"])
        for name, donor in (("polar_s_bb", backbone), ("polar_s", side)):
            limit = float(rules.loc[name, "cutoff_A"])
            inside, outside = (self.flags(donor, ("CYS", "SG", "S"), limit + delta, sulfur_acceptor=True) for delta in (-EPS, EPS))
            self.assertTrue(inside[name], name)
            self.assertFalse(outside[name], name)
        met = ci.SULFUR_ACCEPTORS[("MET", "SD")]
        self.assertTrue(self.flags(backbone, ("MET", "SD", "S"), met - EPS, sulfur_acceptor=True)["polar_s_bb"])
        self.assertFalse(self.flags(backbone, ("MET", "SD", "S"), met + EPS, sulfur_acceptor=True)["polar_s_bb"])
        self.assertIn(f"Met SD {met:g} A", rules.loc["polar_s", "rule"])

    def test_geometric_sulfur_rule_changes_at_its_documented_hydrogen_distance_and_angle(self):
        c_prev, n, ca = np.array([-1.33, 0.0, 0.0]), np.zeros(3), np.array([0.55, 1.3, 0.0])
        h = protonation.backbone_amide_hydrogen(n, c_prev, ca)
        out = (h - n) / np.linalg.norm(h - n)
        perpendicular = np.cross(out, [0.0, 0.0, 1.0])
        rules = ci.class_rules("residue_contacts", sulfur_acceptor=True).set_index("name")
        h_s_limit = float(rules.loc["polar_s_geom", "cutoff_A"])

        def geom(h_s, angle):
            tilt = np.radians(180.0 - angle)                           # angle between the N-H continuation and H->S
            s = h + h_s * (np.cos(tilt) * out + np.sin(tilt) * perpendicular)
            table = atoms((1, "ALA", "C", "C", c_prev), (2, "ALA", "N", "N", n), (2, "ALA", "CA", "C", ca), (9, "CYS", "SG", "S", s))
            result = ci.residue_contacts(table, cutoff=5.0, sulfur_acceptor=True).set_index(["position_a", "position_b"])
            return bool(result.loc[(2, 9), "polar_s_geom"])
        self.assertTrue(geom(h_s_limit - EPS, 180.0))
        self.assertFalse(geom(h_s_limit + EPS, 180.0))
        angle = protonation.SULFUR_ANGLE_MIN
        self.assertTrue(geom(2.5, angle + 1.0))
        self.assertFalse(geom(2.5, angle - 1.0))
        self.assertIn(f"{protonation.SULFUR_ANGLE_MIN:g} degrees", rules.loc["polar_s_geom", "rule"])

    def test_opt_in_rules_are_labeled_and_flagged_as_not_computed_when_off(self):
        off = ci.class_rules("residue_contacts").set_index("name")
        self.assertEqual(off.loc["polar_s", "status"], "opt-in")
        self.assertIn("sulfur_acceptor was off", off.loc["polar_s", "rule"])
        on = ci.class_rules("residue_contacts", sulfur_acceptor=True).set_index("name")
        self.assertNotIn("was off", on.loc["polar_s", "rule"])
        self.assertEqual(off.loc["polar", "status"], "default")


class NeighborRuleTests(unittest.TestCase):
    def test_every_neighbor_class_has_a_documented_cutoff_status_and_selection_note(self):
        rules = ci.class_rules("residue_contacts", neighbor_contacts=True, neighbor_chemistry=("polar", "aromatic")).set_index("name")
        for name in ("polar", "ionic", "hydrophobic", "aromatic"):
            self.assertEqual(rules.loc[f"neighbor_{name}", "status"], "opt-in")
            self.assertEqual(float(rules.loc[f"neighbor_{name}", "cutoff_A"]), 3.5 if name == "polar" else 4.0)
        self.assertNotIn("not selected", rules.loc["neighbor_polar", "rule"])
        self.assertIn("not selected", rules.loc["neighbor_ionic", "rule"])
        self.assertIn("not selected", rules.loc["neighbor_hydrophobic", "rule"])
        self.assertIn("H...A <= 2.5", rules.loc["neighbor_polar", "rule"])
        # the cutoffs themselves are probed with exact geometry in tests/test_neighbor_contacts.py (WindowTests), reading these same rows

    def test_the_legend_says_when_the_option_was_off(self):
        off = ci.class_rules("residue_contacts").set_index("name")
        for name in ("neighbor", "neighbor_polar", "neighbor_aromatic"):
            self.assertIn("neighbor_contacts was off", off.loc[name, "rule"])
            self.assertNotIn("sulfur_acceptor was off", off.loc[name, "rule"])
        self.assertNotIn("was off", ci.class_rules("residue_contacts", neighbor_contacts=True).set_index("name").loc["neighbor_polar", "rule"])
        self.assertNotIn("neighbor_contacts was off", ci.class_rules("residue_contacts").set_index("name").loc["polar_s", "rule"])


class OtherTableRuleTests(unittest.TestCase):
    def features(self, kind, element="CU"):
        return pd.DataFrame([("s", "F", "A1", element, kind, 0.0, 0.0, 0.0)], columns=["structure_id", "feature", "atom", "element", "kind", "x", "y", "z"])

    def test_cofactor_classes_change_at_their_documented_cutoffs(self):
        rules = ci.class_rules("cofactor_contacts").set_index("name")
        cases = (("coord", "metal", ("HIS", "NE2", "N")), ("polar", "oxygen", ("SER", "OG", "O")),
                 ("ionic", "oxygen", ("LYS", "NZ", "N")), ("contact", "oxygen", ("ALA", "CB", "C")))
        for name, kind, residue in cases:
            limit = float(rules.loc[name, "cutoff_A"])
            for delta, expected in ((-EPS, True), (EPS, False)):
                table = atoms((5, residue[0], residue[1], residue[2], (limit + delta, 0.0, 0.0)))
                result = ci.cofactor_contacts(table, self.features(kind), record_distance=8.0)
                self.assertEqual(bool(result[name].iloc[0]), expected, (name, delta))

    def test_cofactor_rows_state_the_record_distance_and_the_tier_rows(self):
        rules = ci.class_rules("cofactor_contacts", record_distance=6.0, coordination_model="by_pair")
        self.assertIn("6 A", rules.set_index("name").loc["(rows)", "rule"])
        self.assertTrue(rules["name"].str.startswith("coord_tier: Cu").any())
        self.assertTrue((rules[rules["name"].str.startswith("coord_tier")]["status"] == "opt-in").all())

    def test_water_bridge_and_relay_cutoffs(self):
        rules = ci.class_rules("water_bridges", bridge=3.2, direct=3.0).set_index("name")
        self.assertEqual((float(rules.loc["bridged", "cutoff_A"]), float(rules.loc["direct", "cutoff_A"])), (3.2, 3.0))
        feature = self.features("oxygen")
        residue = atoms((5, "SER", "OG", "O", (3.0 + EPS, 0.0, 0.0)))
        self.assertFalse(ci.water_bridges(residue, feature, pd.DataFrame(columns=["structure_id", "x", "y", "z"]), bridge=3.2, direct=3.0)["direct"].iloc[0])
        residue = atoms((5, "SER", "OG", "O", (3.0 - EPS, 0.0, 0.0)))
        self.assertTrue(ci.water_bridges(residue, feature, pd.DataFrame(columns=["structure_id", "x", "y", "z"]), bridge=3.2, direct=3.0)["direct"].iloc[0])
        relay = ci.class_rules("relay_steps", link=3.1).set_index("name")
        self.assertEqual(float(relay.loc["steps", "cutoff_A"]), 3.1)

    def test_vote_parameters_appear_in_the_rule(self):
        row = ci.class_rules("conserved_contacts", support=0.8, min_sequences=3, structure_fraction=0.6).set_index("name").loc["conserved_<flag>"]
        for text in ("60%", "3 voting sequences", "80%"):
            self.assertIn(text, row["rule"])


class AttachmentTests(unittest.TestCase):
    def test_results_carry_the_rules_they_were_made_with(self):
        table = ci.residue_contacts(pair(("LEU", "CB", "C"), ("VAL", "CB", "C"), 4.0), cutoff=4.5)
        rules = ci.rules_of(table)
        self.assertEqual(set(rules["table"]), {"residue_contacts"})
        self.assertEqual(float(rules.set_index("name").loc["contact", "cutoff_A"]), 4.5)
        empty = ci.residue_contacts(pair(("LEU", "CB", "C"), ("VAL", "CB", "C"), 9.0))
        self.assertTrue(empty.empty)
        self.assertIn("contact", set(ci.rules_of(empty)["name"]))                  # an empty result still says what it looked for
        features = ci.cofactor_contacts(atoms((5, "HIS", "NE2", "N", (2.0, 0.0, 0.0))), OtherTableRuleTests().features("metal"), coordination_distance=2.5)
        self.assertEqual(float(ci.rules_of(features).set_index("name").loc["coord", "cutoff_A"]), 2.5)
        self.assertTrue(ci.rules_of(pd.DataFrame({"a": [1]})).empty)

    def test_conserved_contacts_carries_the_vote_and_contact_rules(self):
        rows, structures = [], []
        for n in range(6):
            sid = f"s{n}"
            rows += [(sid, 1, 1, "LEU", "CB", "C", 0.0, 0.0, 0.0), (sid, 9, 9, "VAL", "CB", "C", 3.8, 0.0, 0.0)]
            structures.append((sid, "g", f"q{n}"))
        atom_table = pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])
        meta = pd.DataFrame(structures, columns=["structure_id", "group", "sequence_id"])
        result = ci.conserved_contacts(atom_table, meta, support=0.8)
        rules = ci.rules_of(result)
        self.assertIn("conserved_<flag>", set(rules["name"]))
        self.assertIn("80%", rules.set_index("name").loc["conserved_<flag>", "rule"])
        self.assertTrue(result["conserved_contact"].iloc[0])

    def test_table_names_are_checked_and_the_markdown_lists_every_rule(self):
        with self.assertRaises(ValueError):
            ci.class_rules("nope")
        rules = ci.class_rules()
        text = ci.format_rules(rules)
        self.assertEqual(set(rules["table"]), set(ci.RULE_TABLES))
        self.assertEqual(text.count("\n"), len(rules) + 2)
        for name in ("polar_s_geom", "coord_tier", "(rows)", "steps", "bridged"):
            self.assertIn(name, text)


if __name__ == "__main__":
    unittest.main()
