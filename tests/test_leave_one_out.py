"""Leave-one-sequence-out logic: hand-made vote tables with known answers, and agreement with the real vote."""

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


def table(*groups):
    """One key (1, 2); ``groups`` is (group name, [votes]) pairs, one row per sequence."""
    rows = [(1, 2, name, f"{name}{n}", bool(v)) for name, votes in groups for n, v in enumerate(votes)]
    return pd.DataFrame(rows, columns=["position_a", "position_b", "group", "sequence_id", "vote"])


class RuleTests(unittest.TestCase):
    def test_a_group_of_five_needs_every_vote_and_a_group_of_twenty_can_spare_two(self):
        votes = table(("small", [1] * 5), ("big", [1] * 20))
        n, pos = loo.count_tables(votes, ["small", "big"])
        self.assertTrue(loo.passes(n, pos, 0.9, 5).iloc[0])
        margin = loo.margins(votes, ["small", "big"]).iloc[0]
        self.assertEqual((margin["spare_negatives_small"], margin["spare_negatives_big"], margin["margin"]), (0, 2, 0))
        self.assertEqual((margin["spare_voters_small"], margin["spare_voters_big"]), (0, 15))
        missed = table(("small", [1] * 4 + [0]), ("big", [1] * 20))
        n, pos = loo.count_tables(missed, ["small", "big"])
        self.assertFalse(loo.passes(n, pos, 0.9, 5).iloc[0])

    def test_float_rounding_does_not_shave_a_spare_vote_off_a_group_of_ten(self):
        margin = loo.margins(table(("g", [1] * 10)), ["g"], min_sequences=5).iloc[0]
        self.assertEqual(margin["spare_negatives_g"], 1)                      # floor(10 * (1 - 0.9)) is 1, not 0, whatever 1 - 0.9 rounds to


class LeaveOneOutTests(unittest.TestCase):
    def test_a_sequence_that_is_the_single_blocking_negative_is_found(self):
        # g1: five positive sequences. g2: six sequences with one negative, 5 of 6 is 83%: the contact is blocked, and by that sequence alone.
        votes = table(("g1", [1] * 5), ("g2", [1, 1, 1, 1, 1, 0]))
        n, pos = loo.count_tables(votes, ["g1", "g2"])
        self.assertFalse(loo.passes(n, pos, 0.9, 5).iloc[0])
        result = loo.leave_one_out(votes, ["g1", "g2"]).set_index("sequence_id")
        self.assertEqual(result.loc["g25", ["gained_as_is", "gained_relaxed"]].tolist(), [1, 1])       # without it g2 is five of five
        self.assertEqual(result.loc["g20", ["gained_as_is", "gained_relaxed"]].tolist(), [0, 0])       # leaving out a positive does not help
        self.assertEqual(result.loc["g10", ["gained_as_is", "gained_relaxed"]].tolist(), [0, 0])       # and g1 would fall below five voters

    def test_removing_a_voter_from_a_group_at_the_minimum_loses_the_contact_by_the_size_rule_only(self):
        votes = table(("g1", [1] * 5), ("g2", [1] * 6))
        n, pos = loo.count_tables(votes, ["g1", "g2"])
        self.assertTrue(loo.passes(n, pos, 0.9, 5).iloc[0])
        result = loo.leave_one_out(votes, ["g1", "g2"]).set_index("sequence_id")
        self.assertEqual(result.loc["g10", ["lost_as_is", "lost_relaxed"]].tolist(), [1, 0])           # g1 falls to 4 voters; relaxing the minimum rescues it
        self.assertEqual(result.loc["g20", ["lost_as_is", "lost_relaxed"]].tolist(), [0, 0])           # g2 stays at five voters, all positive
        self.assertTrue((result[["gained_as_is", "gained_relaxed"]] == 0).all().all())

    def test_removing_a_positive_can_lose_a_contact_by_support_and_relaxing_the_minimum_does_not_rescue_it(self):
        votes = table(("g", [1] * 9 + [0]))                                    # ten sequences, one negative: 9 of 10 passes
        n, pos = loo.count_tables(votes, ["g"])
        self.assertTrue(loo.passes(n, pos, 0.9, 5).iloc[0])
        result = loo.leave_one_out(votes, ["g"]).set_index("sequence_id")
        self.assertEqual(result.loc["g9", ["lost_as_is", "gained_as_is"]].tolist(), [0, 0])             # leaving the negative out: 9 of 9, still passes
        self.assertEqual(result.loc["g0", ["lost_as_is", "lost_relaxed"]].tolist(), [1, 1])             # leaving a positive out: 8 of 9 is 89%, lost by support


class AgreementWithTheVoteTests(unittest.TestCase):
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

    def test_the_baseline_equals_conserved_contacts_and_the_per_sequence_votes_add_up(self):
        atoms, structures = self.dataset()
        real = ci.conserved_contacts(atoms, structures)
        observations = ci.contact_observations(atoms)
        votes = ci.sequence_votes(observations, structures, loo.KEYS, "contact")
        n, pos = loo.count_tables(votes, ["g1", "g2"])
        mine = loo.passes(n, pos, 0.9, 5)
        self.assertEqual(set(mine.index[mine]), set(map(tuple, real.loc[real["conserved_contact"], loo.KEYS].to_numpy().tolist())))
        merged = real.set_index(loo.KEYS)
        for group in ("g1", "g2"):
            self.assertTrue((merged.loc[n.index, f"contact_n_{group}"].to_numpy() == n[group].to_numpy()).all())
            np.testing.assert_allclose(merged.loc[n.index, f"contact_{group}"].to_numpy(), (pos[group] / n[group]).round(3).to_numpy())
        # leaving out the sequence that breaks contacts frees them
        result = loo.leave_one_out(votes, ["g1", "g2"]).set_index("sequence_id")
        self.assertGreater(result.loc["q3", "gained_relaxed"], 0)
        self.assertEqual(result["gained_relaxed"].drop("q3").sum(), 0)


if __name__ == "__main__":
    unittest.main()
