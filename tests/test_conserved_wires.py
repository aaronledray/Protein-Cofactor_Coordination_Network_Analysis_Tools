"""Tests for water-aware wires across structures (synthetic geometry)."""

import unittest

import pandas as pd

from modules import conserved_wires as cw


def atoms_frame(rows):
    return pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])


def wframe(rows):
    return pd.DataFrame(rows, columns=["structure_id", "x", "y", "z"])


def toy(sid="s"):
    """Source at the origin; residue 5 (side chain) 3 A away; water w0 bridging 5 to residue 9; water w1 giving a second route to 9;
    residue 11 reached only through its backbone, with a side chain that would reach residue 13."""
    atoms = atoms_frame([
        (sid, 5, 5, "ARG", "NH1", "N", 3.0, 0.0, 0.0),
        (sid, 9, 9, "ASP", "OD1", "O", 3.0, 6.0, 0.0),
        (sid, 11, 11, "GLY", "O", "O", -3.0, 0.0, 0.0),
        (sid, 11, 11, "GLY", "N", "N", -5.0, 3.0, 0.0),
        (sid, 13, 13, "SER", "OG", "O", -5.0, 6.4, 0.0),
    ])
    waters = wframe([(sid, 3.0, 3.0, 0.0), (sid, 0.0, 3.0, 0.0), (sid, 1.5, 4.5, 0.0)])
    sources = wframe([(sid, 0.0, 0.0, 0.0)])
    return atoms, sources, waters


class WireTests(unittest.TestCase):
    def test_signatures_steps_and_terminal_backbone(self):
        atoms, sources, waters = toy()
        adj = cw.relay_graph(atoms, sources, waters)
        found = cw.shortest_signatures(adj, [5, 9, 11, 13], max_steps=6)
        self.assertEqual(found[5][0], 1)
        self.assertEqual(found[5][1], [(5,)])
        self.assertEqual(found[9][0], 3)                   # src -> w -> w -> 9 or src -> 5 -> w -> 9 (both three links)
        self.assertIn((5, "w", 9), found[9][1])
        self.assertEqual(found[11][0], 1)                  # reached through the backbone oxygen
        self.assertNotIn(13, found)                        # backbone-only entry does not relay onward

    def test_disjoint_paths_count_independent_routes(self):
        atoms, sources, waters = toy()
        adj = cw.relay_graph(atoms, sources, waters)
        self.assertGreaterEqual(cw.disjoint_path_count(adj, 9), 2)
        only_one = wframe([("s", 3.0, 3.0, 0.0)])
        adj1 = cw.relay_graph(atoms, sources, only_one)
        self.assertEqual(cw.disjoint_path_count(adj1, 9), 1)
        self.assertEqual(cw.disjoint_path_count(adj1, 99), 0)

    def test_vote_across_sequences(self):
        frames, srcs, wats, meta = [], [], [], []
        for g, n in (("A", 3), ("B", 3)):
            for k in range(n):
                sid = f"{g}{k}"
                a, s, w = toy(sid)
                if g == "B" and k == 0:
                    w = w.iloc[:0]                          # no waters modeled in one structure: wire to 9 lost
                frames.append(a), srcs.append(s), wats.append(w)
                meta.append((sid, g, f"seq_{sid}"))
        atoms, sources, waters = pd.concat(frames), pd.concat(srcs), pd.concat(wats)
        structures = pd.DataFrame(meta, columns=["structure_id", "group", "sequence_id"])
        reach, signatures, paths = cw.conserved_wires(atoms, structures, sources, waters, [5, 9], min_sequences=3, support=0.9, max_steps=4)
        r = reach.set_index("target")
        self.assertTrue(r.loc[5, "conserved_within_1"])
        self.assertFalse(r.loc[9, "conserved_within_3"])    # group B is at 2/3
        self.assertAlmostEqual(r.loc[9, "within_3_B"], 2 / 3, places=2)
        loose = cw.conserved_wires(atoms, structures, sources, waters, [5, 9], min_sequences=3, support=0.6, max_steps=4)[0].set_index("target")
        self.assertTrue(loose.loc[9, "conserved_within_3"])
        self.assertIn("signature", signatures.columns)
        self.assertTrue(((paths[paths.target == 5]).disjoint_paths >= 1).all())


if __name__ == "__main__":
    unittest.main()
