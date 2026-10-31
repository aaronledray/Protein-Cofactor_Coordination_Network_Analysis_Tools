"""Known-answer check on spinach plastocyanin (1AG6, 1.6 A X-ray, one Cu, 99 residues, chain A).

The expected values are computed here from the coordinates with Biopython and numpy alone, so they do not
depend on the code under test. They agree with the type 1 copper site in the literature: ligands His37 ND1,
Cys84 SG, His87 ND1 and a weak axial Met92 SD (Cu-S(Met) 2.7 to 2.9 A across plastocyanins), and the
Asn38 backbone NH donating to the Cys84 thiolate. See docs/scope_checks/01_plastocyanin.md.
"""

import itertools
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from modules import conserved_interactions as ci
from modules.coordination_api import analyze_structure

CIF = Path(__file__).resolve().parents[1] / "reference_structures/0_Plastocyanin/1ag6.cif"
LIGANDS = {(37, "HIS", "ND1"), (84, "CYS", "SG"), (87, "HIS", "ND1"), (92, "MET", "SD")}


def _chain():
    from Bio.PDB import MMCIFParser

    return MMCIFParser(QUIET=True).get_structure("pc", str(CIF))[0]["A"]


class PlastocyaninKnownAnswerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chain = _chain()
        cls.residues = [r for r in cls.chain if r.id[0] == " "]
        cls.cu = next(a for r in cls.chain for a in r if a.element == "CU")
        cls.atoms = ci.load_atoms(CIF, "A", structure_id="1ag6")
        cls.features = pd.DataFrame([("1ag6", "Cu", "CU", "CU", "metal", *map(float, cls.cu.coord))],
                                    columns=["structure_id", "feature", "atom", "element", "kind", "x", "y", "z"])

    def distance(self, residue, atom):
        return float(np.linalg.norm(self.chain[residue][atom].coord - self.cu.coord))

    def test_geometry_has_four_ligand_atoms_within_three_angstrom(self):
        near = {(r.id[1], r.resname, a.get_id()) for r in self.residues for a in r
                if a.element != "C" and np.linalg.norm(a.coord - self.cu.coord) <= 3.0}
        self.assertEqual(near, LIGANDS)
        self.assertLess(self.distance(37, "ND1"), 2.1)
        self.assertLess(self.distance(84, "SG"), 2.3)
        self.assertGreater(self.distance(92, "SD"), 2.8)          # weak axial ligand, just beyond the default 2.8 A

    def test_cofactor_contacts_recover_the_ligands_and_depend_on_the_cutoff(self):
        by_cutoff = {}
        for cutoff in (2.8, 3.0):
            table = ci.cofactor_contacts(self.atoms, self.features, coordination_distance=cutoff)
            by_cutoff[cutoff] = set(table.loc[table["coord"], "position"])
        self.assertEqual(by_cutoff[2.8], {37, 84, 87})
        self.assertEqual(by_cutoff[3.0], {37, 84, 87, 92})
        table = ci.cofactor_contacts(self.atoms, self.features)
        self.assertEqual(set(table.loc[table["contact"], "position"]), {36, 37, 84, 87, 92})   # Pro36 O at 3.47 A is a contact, not a ligand

    def test_by_pair_model_tags_three_bonds_and_one_weak_axial_ligand(self):
        plain = ci.cofactor_contacts(self.atoms, self.features)
        tiered = ci.cofactor_contacts(self.atoms, self.features, coordination_model="by_pair")
        pd.testing.assert_frame_equal(tiered[list(plain.columns)], plain)         # opt-in: the default columns do not move
        self.assertEqual(list(plain.columns), ["structure_id", "feature", "position", "resname", "min_distance", *ci.FEATURE_CLASSES])
        tier = tiered.set_index("position")
        self.assertEqual({p: tier.loc[p, "coord_tier"] for p in (37, 84, 87, 92)}, {37: "bond", 84: "bond", 87: "bond", 92: "weak"})
        self.assertEqual(tier.loc[92, "coord_donor_class"], "Met S")
        self.assertEqual(tier.loc[92, "coord_rule"], "Cu-Met S: bond <= 2.7, weak <= 3.3")
        self.assertEqual(tier.loc[36, "coord_tier"], "none")                      # Pro36 O at 3.47 A: a contact, not a ligand
        self.assertEqual(set(tier.index[tier["coord_tier"].isin(["bond", "weak"])]), {37, 84, 87, 92})
        with self.assertRaises(ValueError):
            ci.cofactor_contacts(self.atoms, self.features, coordination_model="other")

    def test_shell_analysis_finds_primary_and_secondary_shell(self):
        result = analyze_structure(CIF, "CU", direct_coordination=True)
        shells = result["residues"].groupby("shell")["residue_number"].apply(set).to_dict()
        self.assertEqual(shells["PCS"], {36, 37, 84, 87, 92})
        self.assertEqual(shells["SCS"], {38, 85})                  # the two backbone NH groups near the thiolate
        contacts = result["contacts"]
        direct = contacts[contacts["direct_coordination"]]
        self.assertEqual(set(direct["dst_resnum"]), {37, 84, 87})   # the 2.6 A default leaves Met92 (2.88 A) as a shell contact
        reported = {(r.dst_resnum, r.dst_atom): r.distance_A for r in contacts[contacts["src_shell"] == "Cofactor"].itertuples()}
        for residue, atom in ((37, "ND1"), (84, "SG"), (87, "ND1"), (92, "SD")):
            self.assertAlmostEqual(reported[(residue, atom)], self.distance(residue, atom), places=3)

    def test_analyze_structure_contacts_can_be_tiered_without_changing_them(self):
        from modules import coordination_table as ct

        contacts = analyze_structure(CIF, "CU", direct_coordination=True)["contacts"]
        tiered = ct.annotate_contacts(contacts)
        pd.testing.assert_frame_equal(tiered[list(contacts.columns)], contacts)
        rows = tiered[tiered["coord_tier"].isin(["bond", "weak"])]
        self.assertEqual({(r.dst_resnum, r.dst_atom, r.coord_tier) for r in rows.itertuples()},
                         {(37, "ND1", "bond"), (84, "SG", "bond"), (87, "ND1", "bond"), (92, "SD", "weak")})
        shell_rows = tiered[tiered["link_type"] != "cofactor->pcs"]
        self.assertTrue((shell_rows["coord_tier"] == "").all())
        self.assertEqual(len(ct.annotate_contacts(contacts.iloc[0:0])), 0)

    def test_neighbor_contacts_are_none_for_plastocyanin_once_a_hydrogen_bond_must_point_at_its_acceptor(self):
        # Asp42 OD1 is 2.98 A from the N of Glu43, but that N-H points away (H...O 3.01 A, 79 degrees): a close pair, not a hydrogen bond.
        typed = ci.residue_contacts(self.atoms, neighbor_contacts=True, neighbor_chemistry=ci.NEIGHBOR_CHEMISTRY_ALL,
                                    metals=ci.load_metals(CIF, structure_id="1ag6"))
        self.assertEqual(int(typed["neighbor"].sum()), 0)
        self.assertEqual(len(typed), len(ci.residue_contacts(self.atoms)))
        loose = {(r.position_a, r.position_b) for r in ci.residue_contacts(self.atoms, min_separation=1).itertuples()}
        self.assertNotIn((42, 43), loose)                                  # min_separation=1 cannot show it either: the peptide bond is the closest pair

    def test_residue_contacts_equal_a_brute_force_recomputation(self):
        coords = {r.id[1]: np.array([a.coord for a in r]) for r in self.residues}
        expected = {}
        for a, b in itertools.combinations(sorted(coords), 2):
            if b - a < 2:
                continue
            closest = np.linalg.norm(coords[a][:, None] - coords[b][None], axis=-1).min()
            if 1.9 <= closest <= 4.0:
                expected[(a, b)] = closest
        found = {(r.position_a, r.position_b): r.min_distance for r in ci.residue_contacts(self.atoms).itertuples()}
        self.assertEqual(set(found), set(expected))
        for key, value in expected.items():
            self.assertAlmostEqual(found[key], value, places=6)

    def test_second_shell_hydrogen_bonds(self):
        self.assertLess(float(np.linalg.norm(self.chain[38]["N"].coord - self.chain[84]["SG"].coord)), 3.6)   # Asn38 NH...S thiolate
        contacts = ci.residue_contacts(self.atoms).set_index(["position_a", "position_b"])
        clasp = contacts.loc[(38, 85)]                               # Asn38 OD1 ... Ser85 N
        self.assertTrue(clasp["polar_sc"])
        self.assertLess(clasp["min_distance"], 3.0)
        self.assertNotIn((84, 85), contacts.index)                   # sequence neighbors are excluded
        # The polar flag needs an N/O pair, so the NH...S hydrogen bond to the thiolate is not typed polar by default.
        self.assertFalse(contacts.loc[(38, 84), "polar"])

    def test_sulfur_acceptor_option_types_the_thiolate_hydrogen_bond(self):
        plain = ci.residue_contacts(self.atoms)
        typed = ci.residue_contacts(self.atoms, sulfur_acceptor=True)
        pd.testing.assert_frame_equal(typed[list(plain.columns)], plain)         # opt-in: existing columns do not move
        typed = typed.set_index(["position_a", "position_b"])
        self.assertTrue(typed.loc[(38, 84), "polar_s"] and typed.loc[(38, 84), "polar_s_bb"])      # Asn38 N...SG 3.45 A
        self.assertNotIn((84, 85), typed.index)                                    # Ser85 N...SG (3.58 A) is still hidden by the neighbor rule (TODO D)
        expected = set()                                                          # brute force from the coordinates (1AG6 has no disulfide)
        for sulfur in self.residues:
            window = {"CYS": 3.7, "MET": 3.6}.get(sulfur.resname)
            atom = sulfur["SG"] if sulfur.resname == "CYS" else sulfur["SD"] if sulfur.resname == "MET" else None
            for other in self.residues:
                if window and abs(other.id[1] - sulfur.id[1]) >= 2 and other.resname != "PRO" and \
                        float(np.linalg.norm(other["N"].coord - atom.coord)) <= window:
                    expected.add(tuple(sorted((sulfur.id[1], other.id[1]))))
        self.assertEqual(set(typed.index[typed["polar_s_bb"]]), expected)
        self.assertIn((38, 84), expected)
        self.assertTrue(typed.loc[(37, 92), "polar_s"])                             # without the metal, His37 ND1 (3.40 A from Met92 SD) reads as a capable donor
        self.assertFalse(typed.loc[(84, 87), "polar_s"])                            # His87 ND1 is 3.603 A from Cys84 SG: beyond the 3.6 A side-chain window
        with_metals = ci.residue_contacts(self.atoms, sulfur_acceptor=True, metals=ci.load_metals(CIF, structure_id="1ag6"))
        with_metals = with_metals.set_index(["position_a", "position_b"])
        self.assertTrue(typed.loc[(38, 84), "polar_s_geom"] and with_metals.loc[(38, 84), "polar_s_geom"])   # the Asn38 N-H points at SG
        self.assertGreater(typed["polar_s_geom"].sum(), 0)
        self.assertTrue(with_metals.loc[(38, 84), "polar_s_bb"])                    # the Asn38 amide is untouched
        self.assertFalse(with_metals.loc[(37, 92), "polar_s"])                      # Met92 SD...His37 ND1 (3.40 A): ND1 is bonded to the copper
        self.assertFalse(with_metals.loc[(37, 92), "polar_s_geom"])
        self.assertTrue(set(typed.index[typed["polar_s_bb"]]) <= set(typed.index[typed["polar_s"]]))


if __name__ == "__main__":
    unittest.main()
