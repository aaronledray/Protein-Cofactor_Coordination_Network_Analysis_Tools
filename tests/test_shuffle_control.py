"""The shuffle control must shuffle only what it says, and a real conserved set must collapse under it."""

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from modules import conserved_interactions as ci

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("shuffle_control", ROOT / "scripts/shuffle_control.py")
shuffle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(shuffle)


def dataset(n_structures=8, n_positions=14):
    """Residues on a zig-zag so that positions i and i+3 touch (CB atoms 3.8 A apart) in every structure, in two groups of sequences."""
    rows, meta = [], []
    for s in range(n_structures):
        sid = f"s{s}"
        meta.append((sid, "g1" if s < n_structures // 2 else "g2", f"q{s}"))
        for p in range(1, n_positions + 1):
            x = 1.9 * p
            y = 3.6 if p % 2 else 0.0
            rows.append((sid, p, p, "LEU", "CB", "C", x, y, 0.0))
    atoms = pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])
    return atoms, pd.DataFrame(meta, columns=["structure_id", "group", "sequence_id"])


class PermutationTests(unittest.TestCase):
    def setUp(self):
        self.atoms, self.structures = dataset()

    def test_full_permutation_keeps_geometry_and_sequence_order_and_permutes_each_structure(self):
        out = shuffle.permute_labels(self.atoms, np.random.default_rng(1))
        for column in ("structure_id", "seq_index", "resname", "atom", "element", "x", "y", "z"):
            pd.testing.assert_series_equal(out[column], self.atoms[column])
        for sid, group in out.groupby("structure_id"):
            self.assertEqual(sorted(group["position"]), sorted(self.atoms.loc[self.atoms["structure_id"] == sid, "position"]))
        self.assertGreater((out["position"] != self.atoms["position"]).mean(), 0.7)                    # most labels moved
        again = shuffle.permute_labels(self.atoms, np.random.default_rng(1))
        pd.testing.assert_frame_equal(out, again)                                                      # same seed, same shuffle
        other = shuffle.permute_labels(self.atoms, np.random.default_rng(2))
        self.assertFalse(out["position"].equals(other["position"]))

    def test_each_structure_gets_its_own_permutation(self):
        out = shuffle.permute_labels(self.atoms, np.random.default_rng(3))
        maps = [tuple(g["position"]) for _, g in out.groupby("structure_id")]
        self.assertGreater(len(set(maps)), 1)

    def test_window_permutation_stays_inside_blocks(self):
        out = shuffle.permute_labels(self.atoms, np.random.default_rng(4), window=3)
        for sid, group in out.groupby("structure_id"):
            before = self.atoms.loc[self.atoms["structure_id"] == sid, "position"].to_numpy()
            after = group["position"].to_numpy()
            for start in range(0, 14, 3):
                self.assertEqual(sorted(after[start:start + 3]), sorted(before[start:start + 3]))
        self.assertTrue((out["position"] != self.atoms["position"]).any())

    def test_only_limits_the_shuffle_to_the_named_structures(self):
        out = shuffle.permute_labels(self.atoms, np.random.default_rng(5), only={"s0", "s1"})
        rest = self.atoms["structure_id"].isin(["s0", "s1"]) == False  # noqa: E712
        pd.testing.assert_series_equal(out.loc[rest, "position"], self.atoms.loc[rest, "position"])
        self.assertTrue((out.loc[~rest, "position"] != self.atoms.loc[~rest, "position"]).any())

    def test_scrambling_sequences_changes_exactly_the_chosen_ones_with_one_permutation_each(self):
        out, chosen = shuffle.scramble_sequences(self.atoms, self.structures, 2, np.random.default_rng(6))
        self.assertEqual(len(chosen), 2)
        moved = set(self.structures.loc[self.structures["sequence_id"].isin(chosen), "structure_id"])
        changed = set(out.loc[out["position"] != self.atoms["position"], "structure_id"])
        self.assertTrue(changed <= moved and changed)
        untouched = ~out["structure_id"].isin(moved)
        pd.testing.assert_series_equal(out.loc[untouched, "position"], self.atoms.loc[untouched, "position"])
        with self.assertRaises(ValueError):
            shuffle.scramble_sequences(self.atoms, self.structures, 99, np.random.default_rng(1))


class CollapseTests(unittest.TestCase):
    def test_a_real_conserved_set_collapses_under_the_full_shuffle_and_survives_a_single_misaligned_sequence(self):
        atoms, structures = dataset(n_structures=12)
        groups = ["g1", "g2"]
        real, _ = shuffle.conserved_set(atoms, structures, groups, min_sequences=5)
        self.assertGreater(len(real), 5)
        for seed in range(5):
            shuffled = shuffle.permute_labels(atoms, np.random.default_rng(seed))
            found, _ = shuffle.conserved_set(shuffled, structures, groups, min_sequences=5)
            self.assertEqual(len(found & real), 0)
            self.assertLessEqual(len(found), 1)
        # twelve sequences in two groups of six: the 90% rule needs 6 of 6 per group, so one scrambled sequence keeps only the pairs it
        # happens to preserve by chance: a subset of the real set, and a small one
        one, _ = shuffle.scramble_sequences(atoms, structures, 1, np.random.default_rng(0))
        survivors = shuffle.conserved_set(one, structures, groups, min_sequences=5)[0]
        self.assertTrue(survivors <= real)
        self.assertLess(len(survivors), len(real) / 2)
        # with support at 80% one misaligned sequence in six is tolerated
        tolerant, _ = shuffle.conserved_set(one, structures, groups, min_sequences=5, support=0.8)
        self.assertGreater(len(tolerant & real), 0)

    def test_the_plan_names_every_scheme_and_replicate(self):
        tasks = shuffle.plan(2, [3], [1, 2])
        self.assertEqual([t[0] for t in tasks], ["full", "full", "local-3", "local-3", "sequences-1", "sequences-1", "sequences-2", "sequences-2"])


if __name__ == "__main__":
    unittest.main()
