"""The opt-in neighbor-contact class (residues one apart in sequence, judged per atom pair).

Real proteins are recomputed independently from the coordinates (a different bond-counting method and explicit donor and
acceptor lists); synthetic chains with exact geometry test each window and exclusion.
"""

import unittest
from pathlib import Path

import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import neighbor_geometry as geometry  # noqa: E402

from modules import conserved_interactions as ci  # noqa: E402
from modules import protonation  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
STRUCTURES = {"1ag6": (ROOT / "reference_structures/0_Plastocyanin/1ag6.cif", "A"),
              "1a6m": (ROOT / "reference_structures/3_Myoglobin/1a6m.pdb", "A"),
              "1ca2": (ROOT / "reference_structures/3_Metals/1ca2.pdb", "A")}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
DONORS = {("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("ASN", "ND2"), ("GLN", "NE2"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
          ("LYS", "NZ"), ("TRP", "NE1"), ("HIS", "ND1"), ("HIS", "NE2")}
RINGS = {"PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}, "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
         "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}, "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"}}
POSITIVE = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2")}
NEGATIVE = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}


def residues_of(structure_path, chain_id):
    from Bio.PDB import MMCIFParser, PDBParser

    parser = MMCIFParser(QUIET=True) if str(structure_path).endswith(".cif") else PDBParser(QUIET=True)
    model = parser.get_structure("x", str(structure_path))[0]
    return [r for r in model[chain_id] if r.id[0] == " "]


def expected_neighbor_rows(residues, classes):
    """{(first, second): {class, ...}} by brute force: one bond graph over both residues and explicit definitions of each class."""
    expected = {}
    beyond = lambda n: n not in BACKBONE and n != "CB"
    for first, second in zip(residues, residues[1:]):
        if second.id[1] - first.id[1] != 1 or np.linalg.norm(first["C"].coord - second["N"].coord) >= 1.6:
            continue
        atoms = list(first) + list(second)
        coords = np.array([a.coord for a in atoms])
        bonds = np.where(np.linalg.norm(coords[:, None] - coords[None], axis=-1) <= 1.9, 1, 99)
        np.fill_diagonal(bonds, 0)
        for k in range(len(atoms)):                                      # shortest bond counts (Floyd-Warshall on a tiny graph)
            bonds = np.minimum(bonds, bonds[:, [k]] + bonds[[k], :])
        for x, a in enumerate(first):
            for y, b in enumerate(second, start=len(first)):
                if bonds[x, y] < 4:
                    continue
                ida, idb, d = a.get_id(), b.get_id(), float(np.linalg.norm(a.coord - b.coord))
                ra, rb = first.resname, second.resname
                hits = {
                    "hydrophobic": a.element in ("C", "S") and b.element in ("C", "S") and beyond(ida) and beyond(idb) and d <= 4.0,
                    "aromatic": ((ida in RINGS.get(ra, ()) and beyond(idb)) or (idb in RINGS.get(rb, ()) and beyond(ida))) and d <= 4.0,
                    "ionic": (((ra, ida) in POSITIVE and (rb, idb) in NEGATIVE) or ((ra, ida) in NEGATIVE and (rb, idb) in POSITIVE)) and d <= 4.0,
                }
                for name in classes:
                    if hits[name]:
                        expected.setdefault((first.id[1], second.id[1]), set()).add(name)
    return expected


class RealProteinTests(unittest.TestCase):
    def test_hydrophobic_aromatic_and_ionic_classes_match_an_independent_recomputation(self):
        classes = ("ionic", "hydrophobic", "aromatic")
        for label, (path, chain) in STRUCTURES.items():
            residues = residues_of(path, chain)
            expected = expected_neighbor_rows(residues, classes)
            atoms = ci.load_atoms(path, chain, structure_id=label)
            table = ci.residue_contacts(atoms, neighbor_contacts=True, neighbor_chemistry=classes)
            found = {(r.position_a, r.position_b): {c for c in classes if getattr(r, c)} for r in table[table["neighbor"]].itertuples()}
            self.assertEqual(found, expected, label)
            self.assertTrue(len(expected) > 0 or label == "1ag6", label)                    # the rows are not trivially empty

    def test_polar_neighbor_contacts_are_hydrogen_bonds_and_follow_the_stated_geometry(self):
        for label, (path, chain) in STRUCTURES.items():
            atoms = ci.load_atoms(path, chain, structure_id=label)
            default = ci.residue_contacts(atoms)
            table = ci.residue_contacts(atoms, neighbor_contacts=True, neighbor_chemistry=("polar",))
            pd.testing.assert_frame_equal(table[~table["neighbor"]][list(default.columns)].reset_index(drop=True), default)    # nothing else moves
            rows = table[table["neighbor"]]
            self.assertTrue((rows["polar"] & rows["polar_sc"] & ~rows["polar_bb_bb"] & ~rows["hydrophobic"] & ~rows["aromatic"] & ~rows["ionic"]).all(), label)
            self.assertTrue((rows["min_distance"] <= 3.5).all(), label)
            self.assertTrue((rows["closest"] != "bb_bb").all(), label)
            # both ways: every row has a donor-acceptor pair with the stated geometry, and every adjacent pair that has one is a row
            residues = {r.id[1]: r for r in residues_of(path, chain)}
            for row in rows.itertuples():
                self.assertTrue(self.has_hydrogen_bond(residues, row.position_a, row.position_b), (label, row.position_a, row.position_b))
            adjacent = [(n, n + 1) for n in residues if n + 1 in residues and np.linalg.norm(residues[n]["C"].coord - residues[n + 1]["N"].coord) < 1.6]
            expected = {pair for pair in adjacent if self.has_hydrogen_bond(residues, *pair)}
            self.assertEqual({(r.position_a, r.position_b) for r in rows.itertuples()}, expected, label)
        proline = ci.load_atoms(STRUCTURES["1a6m"][0], "A", structure_id="1a6m")
        rows = ci.residue_contacts(proline, neighbor_contacts=True, neighbor_chemistry=("polar",))
        by_second = proline.drop_duplicates("position").set_index("position")["resname"]
        self.assertNotIn("PRO", set(by_second.loc[rows.loc[rows["neighbor"], "position_b"]]))                               # a Pro N has no hydrogen to donate

    def has_hydrogen_bond(self, residues, first, second):
        low, high = residues[first], residues[second]
        for donor_residue, acceptor_residue in ((low, high), (high, low)):
            before = residues.get(donor_residue.id[1] - 1)
            for donor in donor_residue:
                key = (donor_residue.resname, donor.get_id())
                if donor.get_id() == "N" and donor_residue.resname != "PRO":
                    if before is None or "C" not in before or np.linalg.norm(before["C"].coord - donor.coord) >= 1.6:
                        continue
                    hydrogens = [protonation.backbone_amide_hydrogen(donor.coord, before["C"].coord, donor_residue["CA"].coord)]
                elif key in protonation.RING_DONORS:
                    a, b = protonation.RING_DONORS[key]
                    hydrogens = [protonation.ring_nh_hydrogen(donor.coord, donor_residue[a].coord, donor_residue[b].coord)] if a in donor_residue and b in donor_residue else []
                elif key in protonation.AMIDE_DONORS:
                    carbon, reference = protonation.AMIDE_DONORS[key]
                    hydrogens = protonation.amide_nh2_hydrogens(donor.coord, donor_residue[carbon].coord, donor_residue[reference].coord) if carbon in donor_residue and reference in donor_residue else []
                elif key in DONORS:
                    hydrogens = None                                                                     # hydroxyl or Lys NZ: distance only
                else:
                    continue
                for acceptor in acceptor_residue:
                    if not (acceptor.element == "O" or (acceptor_residue.resname == "HIS" and acceptor.get_id() in ("ND1", "NE2"))):
                        continue
                    if donor.get_id() in BACKBONE and acceptor.get_id() in BACKBONE:
                        continue                                                                         # needs a side-chain atom
                    if np.linalg.norm(donor.coord - acceptor.coord) > 3.5:
                        continue
                    if hydrogens is None or not hydrogens:
                        if hydrogens is None:
                            return True
                        continue
                    for h in hydrogens:
                        h_a, angle = protonation.hydrogen_bond_geometry(h, donor.coord, acceptor.coord)
                        if h_a <= 2.5 and angle >= 120.0:
                            return True
        return False


class WindowTests(unittest.TestCase):
    """Each class changes at its documented distance (probed 0.01 A either side) with exact geometry; nothing else moves."""

    def rows(self, atoms, classes=ci.NEIGHBOR_CHEMISTRY_ALL, **kwargs):
        table = ci.residue_contacts(atoms, neighbor_contacts=True, neighbor_chemistry=classes, **kwargs)
        return table[table["neighbor"]]

    def window(self, name):
        rules = ci.class_rules("residue_contacts", neighbor_contacts=True, neighbor_chemistry=ci.NEIGHBOR_CHEMISTRY_ALL).set_index("name")
        return float(rules.loc[f"neighbor_{name}", "cutoff_A"])

    def test_the_polar_window_for_a_hydroxyl_donor(self):
        limit = self.window("polar")
        self.assertEqual(limit, 3.5)
        inside = self.rows(geometry.hydroxyl_to_carbonyl(limit - 0.01))
        self.assertEqual(len(inside), 1)
        row = inside.iloc[0]
        self.assertTrue(row["polar"] and row["polar_sc"] and not row["polar_bb_bb"] and row["neighbor"])
        self.assertFalse(row["hydrophobic"] or row["aromatic"] or row["ionic"])
        self.assertAlmostEqual(row["min_distance"], limit - 0.01, places=2)
        self.assertTrue(self.rows(geometry.hydroxyl_to_carbonyl(limit + 0.01)).empty)

    def test_the_hydrophobic_window_and_c_beta_exclusion(self):
        limit = self.window("hydrophobic")
        self.assertEqual(limit, 4.0)
        row = self.rows(geometry.leucine_valine(limit - 0.01)).iloc[0]
        self.assertTrue(row["hydrophobic"])
        self.assertFalse(row["polar"] or row["aromatic"] or row["ionic"])
        self.assertTrue(self.rows(geometry.leucine_valine(limit + 0.01)).empty)
        # the two C-beta atoms are 2.8 A apart in this chain and the backbone is everywhere close: none of that is hydrophobic chemistry
        self.assertTrue(self.rows(geometry.hydroxyl_to_carbonyl(5.0)).empty)

    def test_the_aromatic_window_needs_a_ring_atom_against_a_side_chain_atom_beyond_c_beta(self):
        limit = self.window("aromatic")
        row = self.rows(geometry.phenylalanine_serine(limit - 0.01)).iloc[0]
        self.assertTrue(row["aromatic"])
        self.assertFalse(row["hydrophobic"] or row["polar"])                      # Phe CG is a carbon, Ser OG an oxygen
        self.assertTrue(self.rows(geometry.phenylalanine_serine(limit + 0.01)).empty)

    def test_only_the_selected_classes_are_reported(self):
        near = geometry.leucine_valine(3.8)
        self.assertEqual(len(self.rows(near, ("hydrophobic",))), 1)
        self.assertTrue(self.rows(near, ("polar", "ionic")).empty)                # the default classes do not include it
        self.assertTrue(self.rows(near, ("aromatic",)).empty)
        mixed = self.rows(geometry.phenylalanine_serine(3.8), ("aromatic", "polar"))
        self.assertTrue(mixed.iloc[0]["aromatic"] and not mixed.iloc[0]["polar"])
        for bad in ((), ("polar", "nope")):
            with self.assertRaises(ValueError):
                ci.residue_contacts(near, neighbor_contacts=True, neighbor_chemistry=bad)

    def test_a_metal_coordinated_donor_is_not_a_donor(self):
        atoms = geometry.hydroxyl_to_carbonyl(3.3)
        self.assertEqual(len(self.rows(atoms)), 1)
        og = atoms[atoms["atom"] == "OG"][["x", "y", "z"]].to_numpy(float)[0]
        metals = pd.DataFrame([("s", "ZN", *(og + [0.0, 0.0, 2.0]))], columns=["structure_id", "element", "x", "y", "z"])
        self.assertTrue(self.rows(atoms, metals=metals).empty)                    # Zn 2.0 A from OG: the hydroxyl is deprotonated and bound
        far = pd.DataFrame([("s", "ZN", *(og + [0.0, 0.0, 4.0]))], columns=["structure_id", "element", "x", "y", "z"])
        self.assertEqual(len(self.rows(atoms, metals=far)), 1)

    def test_the_ionic_window_and_the_donor_acceptor_rules_at_function_level(self):
        names = np.array(["NZ", "OD1", "N", "O"])
        res = np.array(["LYS", "ASP", "PRO", "GLY"])
        elem = np.array(["N", "O", "N", "O"])
        bb = np.array([False, False, True, True])
        xyz = np.zeros((4, 3))
        nohydrogen = np.full((4, 2, 3), np.nan)
        none = np.zeros(4, bool)

        def flags(a, b, d, excluded=none):
            return ci._neighbor_flags(elem, res, names, bb, np.array([a]), np.array([b]), np.array([d]), 5.0, xyz, nohydrogen, excluded)
        self.assertTrue(flags(0, 1, 3.99)["ionic"][0])
        self.assertFalse(flags(0, 1, 4.01)["ionic"][0])
        self.assertTrue(flags(0, 1, 3.0)["polar"][0])                              # a Lys NZ (donor) to an Asp OD1 (acceptor), no hydrogen to test
        self.assertFalse(flags(0, 1, 3.0, excluded=np.array([True, False, False, False]))["polar"][0])
        self.assertFalse(flags(2, 3, 3.0)["polar"][0])                             # a Pro N has no hydrogen, and backbone-backbone never counts
        self.assertFalse(flags(1, 3, 3.0)["polar"][0])                             # two acceptors

    def test_the_hydrogen_must_point_at_the_acceptor_when_it_can_be_placed(self):
        names, res = np.array(["N", "OD1"]), np.array(["GLU", "ASP"])
        elem, bb = np.array(["N", "O"]), np.array([True, False])
        xyz = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        pointing, away = np.array([[[1.0, 0.0, 0.0], [np.nan] * 3], [[np.nan] * 3] * 2]), np.array([[[-1.0, 0.0, 0.0], [np.nan] * 3], [[np.nan] * 3] * 2])
        none = np.zeros(2, bool)

        def polar(hydrogens, d):
            return bool(ci._neighbor_flags(elem, res, names, bb, np.array([0]), np.array([1]), np.array([d]), 5.0, xyz, hydrogens, none)["polar"][0])
        self.assertTrue(polar(pointing, 3.0))                                     # H...O 2.0 A, angle 180
        self.assertFalse(polar(away, 3.0))                                        # same N...O, but the N-H points away: a close pair is not a hydrogen bond
        self.assertFalse(polar(np.array([[[0.0, 1.0, 0.0], [np.nan] * 3], [[np.nan] * 3] * 2]), 3.0))   # perpendicular


class VoteTests(unittest.TestCase):
    def test_the_vote_carries_neighbor_contacts_when_asked(self):
        frames, structures = [], []
        for n in range(6):
            frame = geometry.leucine_valine(3.8)
            frame["structure_id"] = f"s{n}"
            frames.append(frame)
            structures.append((f"s{n}", "g", f"q{n}"))
        atoms = pd.concat(frames, ignore_index=True)
        meta = pd.DataFrame(structures, columns=["structure_id", "group", "sequence_id"])
        self.assertTrue(ci.conserved_contacts(atoms, meta).empty)
        self.assertTrue(ci.conserved_contacts(atoms, meta, neighbor_contacts=True).empty)                  # the default classes do not include hydrophobic
        voted = ci.conserved_contacts(atoms, meta, neighbor_contacts=True, neighbor_chemistry=("hydrophobic",))
        self.assertEqual(len(voted), 1)
        self.assertTrue(voted["conserved_contact"].iloc[0] and voted["conserved_hydrophobic"].iloc[0])
        rules = ci.rules_of(voted).set_index("name")
        self.assertIn("neighbor_hydrophobic", rules.index)
        self.assertNotIn("not selected", rules.loc["neighbor_hydrophobic", "rule"])
        self.assertIn("not selected", rules.loc["neighbor_polar", "rule"])


if __name__ == "__main__":
    unittest.main()
