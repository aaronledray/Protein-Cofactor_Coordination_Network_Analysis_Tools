"""vote_coverage and conserved_contacts(coverage=True): the depth of the vote behind each contact."""

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from modules import conserved_interactions as ci

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("leave_one_out", ROOT / "scripts/leave_one_out.py")
loo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(loo)


def votes(key, *groups):
    """Votes of one key (a pair of positions); ``groups`` is (group name, [0/1 votes]) pairs, sequences named <group><index>."""
    return pd.DataFrame([(key[0], key[1], g, f"{g}{n}", bool(v)) for g, vs in groups for n, v in enumerate(vs)],
                        columns=["position_a", "position_b", "group", "sequence_id", "vote"])


class ArithmeticTests(unittest.TestCase):
    def test_spare_negatives_agree_with_the_pass_rule_for_every_group_size_and_support(self):
        for support in (0.9, 0.75):
            for n in range(1, 51):
                for positives in range(n + 1):
                    frame = votes((1, 2), ("g", [1] * positives + [0] * (n - positives)))
                    row = ci.vote_coverage(frame, ["g"], support=support, min_sequences=1).iloc[0]
                    passes = positives / n >= support
                    self.assertEqual(row["spare_negatives_g"] >= 0, passes, (support, n, positives))
                    self.assertEqual(row["failing_groups"] == "", passes, (support, n, positives))


class HandMadeTests(unittest.TestCase):
    def test_a_single_negative_that_holds_a_group_back_is_named(self):
        row = ci.vote_coverage(votes((1, 2), ("g1", [1] * 5), ("g2", [1, 1, 1, 1, 1, 0])), ["g1", "g2"]).iloc[0]
        self.assertEqual((row["conserved_in"], row["failing_groups"], row["held_back_by"]), ("g1", "g2", "g25"))
        self.assertEqual((row["spare_negatives_g1"], row["spare_negatives_g2"]), (0, -1))
        self.assertEqual((row["spare_voters_g1"], row["spare_voters_g2"]), (0, 1))

    def test_a_contact_that_passes_everywhere_reports_its_margins(self):
        row = ci.vote_coverage(votes((1, 2), ("small", [1] * 5), ("big", [1] * 20)), ["small", "big"]).iloc[0]
        self.assertEqual((row["conserved_in"], row["failing_groups"], row["held_back_by"]), ("small+big", "", ""))
        self.assertEqual((row["spare_negatives_small"], row["spare_negatives_big"]), (0, 2))
        self.assertEqual((row["spare_voters_small"], row["spare_voters_big"]), (0, 15))

    def test_a_group_that_is_too_small_is_never_held_back_by_a_sequence(self):
        row = ci.vote_coverage(votes((1, 2), ("g1", [1] * 4), ("g2", [1] * 6)), ["g1", "g2"]).iloc[0]
        self.assertEqual((row["conserved_in"], row["failing_groups"], row["held_back_by"]), ("g2", "g1", ""))
        self.assertEqual(row["spare_voters_g1"], -1)

    def test_dropping_the_negative_of_a_group_at_the_minimum_would_leave_it_too_small_so_nothing_holds_it_back(self):
        row = ci.vote_coverage(votes((1, 2), ("g", [1, 1, 1, 1, 0])), ["g"]).iloc[0]                  # five sequences, one negative: 4 of 5 fails
        self.assertEqual((row["failing_groups"], row["held_back_by"]), ("g", ""))                       # dropping it leaves four voters, below the minimum
        self.assertEqual(row["one_negative_away"], "g4")                                                # but the support would be met: one sequence short
        row = ci.vote_coverage(votes((1, 2), ("g", [1, 1, 1, 1, 1, 0])), ["g"]).iloc[0]                # six sequences: dropping the negative leaves five
        self.assertEqual((row["held_back_by"], row["one_negative_away"]), ("g5", "g5"))                 # strict rescue implies the looser one

    def test_a_small_group_that_is_all_positive_has_no_negative_to_name(self):
        row = ci.vote_coverage(votes((1, 2), ("g1", [1] * 4), ("g2", [1] * 6)), ["g1", "g2"]).iloc[0]
        self.assertEqual((row["failing_groups"], row["held_back_by"], row["one_negative_away"]), ("g1", "", ""))

    def test_two_failing_groups_are_not_held_back_by_one_sequence(self):
        row = ci.vote_coverage(votes((1, 2), ("g1", [1, 1, 1, 1, 1, 0]), ("g2", [1, 1, 1, 1, 1, 0])), ["g1", "g2"]).iloc[0]
        self.assertEqual((row["conserved_in"], row["failing_groups"], row["held_back_by"]), ("", "g1+g2", ""))

    def test_when_several_negatives_could_each_be_dropped_all_are_listed(self):
        # 31 sequences allow 3 negatives; with 4 the group fails, and dropping any one of the four leaves 27 of 30 = 90%
        row = ci.vote_coverage(votes((1, 2), ("g", [1] * 27 + [0] * 4)), ["g"]).iloc[0]
        self.assertEqual(row["failing_groups"], "g")
        self.assertEqual(row["held_back_by"], ";".join(sorted(f"g{n}" for n in range(27, 31))))
        row = ci.vote_coverage(votes((1, 2), ("g", [1] * 26 + [0] * 5)), ["g"]).iloc[0]                 # five negatives: dropping one is not enough
        self.assertEqual(row["held_back_by"], "")

    def test_several_keys_are_kept_apart(self):
        frame = pd.concat([votes((1, 2), ("g", [1] * 6)), votes((3, 4), ("g", [1] * 5 + [0]))], ignore_index=True)
        table = ci.vote_coverage(frame, ["g"]).set_index(["position_a", "position_b"])
        self.assertEqual(table.loc[(1, 2), "conserved_in"], "g")
        self.assertEqual((table.loc[(3, 4), "failing_groups"], table.loc[(3, 4), "held_back_by"]), ("g", "g5"))


class WithTheVoteTests(unittest.TestCase):
    def dataset(self):
        rows, meta = [], []
        for s in range(12):
            sid = f"s{s}"
            meta.append((sid, "g1" if s < 6 else "g2", f"q{s}"))
            for p in range(1, 15):
                rows.append((sid, p, p, "LEU", "CB", "C", 1.9 * p, 3.6 if p % 2 else 0.0, 0.0))
        atoms = pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])
        atoms.loc[(atoms["structure_id"] == "s3") & (atoms["position"] == 5), "y"] += 9.0       # one sequence breaks a few contacts
        return atoms, pd.DataFrame(meta, columns=["structure_id", "group", "sequence_id"])

    def test_the_default_output_is_unchanged_and_the_headline_equals_the_vote(self):
        atoms, structures = self.dataset()
        plain = ci.conserved_contacts(atoms, structures)
        detailed = ci.conserved_contacts(atoms, structures, coverage=True)
        pd.testing.assert_frame_equal(detailed[list(plain.columns)], plain)                        # opt-in: nothing else moves
        self.assertEqual(list(detailed.columns[len(plain.columns):]),
                         ["spare_negatives_g1", "spare_voters_g1", "spare_negatives_g2", "spare_voters_g2", "conserved_in", "failing_groups", "held_back_by",
                          "one_negative_away"])
        self.assertTrue(((detailed["failing_groups"] == "") == detailed["conserved_contact"]).all())
        blocked = detailed[detailed["held_back_by"] != ""]
        self.assertGreater(len(blocked), 0)
        self.assertTrue((blocked["held_back_by"] == "q3").all())                                   # the sequence that breaks them is the one named
        self.assertTrue((blocked["one_negative_away"] == blocked["held_back_by"]).all())            # strict rescue implies the looser one, same sequences

    def test_margins_agree_with_the_leave_one_out_script(self):
        atoms, structures = self.dataset()
        observations = ci.contact_observations(atoms)
        table = ci.sequence_votes(observations, structures, loo.KEYS, "contact")
        margins = loo.margins(table, ["g1", "g2"]).reset_index()
        coverage = ci.vote_coverage(table, ["g1", "g2"])
        coverage = coverage[coverage["failing_groups"] == ""].merge(margins, on=loo.KEYS, suffixes=("", "_script"))
        self.assertGreater(len(coverage), 0)
        for group in ("g1", "g2"):
            self.assertTrue((coverage[f"spare_negatives_{group}"] == coverage[f"spare_negatives_{group}_script"]).all())
            self.assertTrue((coverage[f"spare_voters_{group}"] == coverage[f"spare_voters_{group}_script"]).all())

    def test_the_legend_lists_the_columns_and_says_when_they_were_not_computed(self):
        on = ci.class_rules("conserved_contacts", coverage=True).set_index("name")
        for name in ("conserved_in", "failing_groups", "spare_negatives_<group>", "spare_voters_<group>", "held_back_by", "one_negative_away"):
            self.assertEqual(on.loc[name, "status"], "opt-in")
            self.assertNotIn("coverage was off", on.loc[name, "rule"])
        self.assertIn("90%", on.loc["conserved_in", "rule"])
        off = ci.class_rules("conserved_contacts").set_index("name")
        self.assertIn("coverage was off", off.loc["held_back_by", "rule"])
        atoms, structures = self.dataset()
        self.assertIn("held_back_by", set(ci.rules_of(ci.conserved_contacts(atoms, structures, coverage=True))["name"]))


if __name__ == "__main__":
    unittest.main()
