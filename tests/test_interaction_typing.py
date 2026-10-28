import unittest

from modules.interaction_typing import (
    angle_degrees,
    geometry_rule_accepts,
    interaction_type_eligibility,
    normalize_arpeggio_contacts,
    sequence_balanced_type_summary,
    summarize_residue_pairs,
    _endpoint,
)


def atom(comp, chain, number, name, comp_type="polypeptide(L)"):
    return {
        "label_comp_id": comp, "auth_asym_id": chain, "auth_seq_id": number,
        "pdbx_PDB_ins_code": "?", "auth_atom_id": name,
        "label_comp_type": comp_type,
    }


class InteractionTypingTests(unittest.TestCase):
    def test_atom_pair_emits_all_types_and_disulfide_subtype(self):
        left, right = atom("CYS", "A", 12, "SG"), atom("CYS", "A", 46, "SG")
        raw = [{
            "type": "atom-atom", "bgn": left, "end": right, "distance": 2.04,
            "contact": ["covalent", "hydrophobic"],
        }]
        rows = normalize_arpeggio_contacts(raw, structure_id="fixture")
        self.assertEqual({row["interaction_type"] for row in rows}, {"disulfide", "hydrophobic"})
        self.assertTrue(all(row["analysis_layer"] == "direct_protein_protein" for row in rows))
        self.assertTrue(all(row["src_resnum"] == "12" for row in rows))
        self.assertIn("criterion", rows[0])
        self.assertIn("distance_A", rows[0])

    def test_water_and_cofactor_layers_are_independent(self):
        protein = atom("HIS", "A", 10, "NE2")
        water = atom("HOH", "W", 3, "O", "water")
        cofactor = atom("HEM", "A", 500, "FE", "non-polymer")
        rows = normalize_arpeggio_contacts([
            {"type": "atom-atom", "bgn": protein, "end": water, "distance": 2.9, "contact": ["hbond"]},
            {"type": "atom-atom", "bgn": protein, "end": cofactor, "distance": 2.1, "contact": ["metal_complex"]},
        ])
        self.assertEqual({row["analysis_layer"] for row in rows}, {
            "water_mediated_or_water_endpoint", "protein_cofactor_or_ligand",
        })

    def test_ring_plane_geometry_categories_include_offset_stacking(self):
        rows = normalize_arpeggio_contacts([{
            "type": "plane-plane", "bgn": atom("PHE", "A", 10, "CG,CD1,CD2,CE1,CE2,CZ"),
            "end": atom("TYR", "A", 30, "CG,CD1,CD2,CE1,CE2,CZ"),
            "distance": 5.2, "contact": ["OF"],
        }])
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["interaction_type"], "aromatic_stacking")

    def test_pair_summary_keeps_multiple_types_and_baseline_only_edges(self):
        evidence = normalize_arpeggio_contacts([{
            "type": "atom-atom", "bgn": atom("ASP", "A", 1, "OD1"),
            "end": atom("ARG", "A", 2, "NH1"), "distance": 3.0,
            "contact": ["ionic", "hbond"],
        }], structure_id="x")
        baseline = [{
            "structure_id": "x", "src_resname": "GLY", "src_chain": "A", "src_resnum": 3,
            "dst_resname": "TYR", "dst_chain": "A", "dst_resnum": 4,
            "baseline_contact": True,
        }]
        result = summarize_residue_pairs(evidence, baseline_pairs=baseline)
        self.assertIn("hydrogen_bond;ionic", {r["interaction_types"] for r in result})
        self.assertIn("proximity_only_unclassified", {r["outcome"] for r in result})

    def test_sequence_balancing_uses_assessable_type_denominators_and_half_replicates(self):
        rows = [
            {"edge_id": "e1", "group": "A", "sequence_hash": "s1", "eligible": True,
             "eligible_by_type": {"ionic": True, "disulfide": False}, "interaction_types": "ionic"},
            {"edge_id": "e1", "group": "A", "sequence_hash": "s1", "eligible": True,
             "eligible_by_type": {"ionic": True, "disulfide": False}, "interaction_types": ""},
            {"edge_id": "e1", "group": "A", "sequence_hash": "s2", "eligible": True,
             "eligible_by_type": {"ionic": False, "disulfide": True}, "interaction_types": "disulfide"},
        ]
        result = sequence_balanced_type_summary(rows)
        ionic = next(row for row in result if row["interaction_type"] == "ionic")
        disulfide = next(row for row in result if row["interaction_type"] == "disulfide")
        self.assertEqual((ionic["eligible_sequences"], ionic["positive_sequences"]), (1, 1))
        self.assertEqual((disulfide["eligible_sequences"], disulfide["positive_sequences"]), (1, 1))
        self.assertEqual(ionic["eligible_structures"], 2)
        self.assertEqual(disulfide["eligible_structures"], 1)

    def test_angle_helper_covers_valid_and_near_miss_geometry(self):
        self.assertAlmostEqual(angle_degrees((1, 0, 0), (0, 0, 0), (0, 1, 0)), 90.0)
        self.assertAlmostEqual(angle_degrees((1, 0, 0), (0, 0, 0), (-1, 0, 0)), 180.0)
        with self.assertRaises(ValueError):
            angle_degrees((0, 0, 0), (0, 0, 0), (1, 0, 0))

    def test_each_implemented_type_has_a_positive_and_distance_or_geometry_near_miss(self):
        cases = {
            "hydrogen_bond": ({"distance_A": 3.0, "donor_h_acceptor_angle_deg": 120}, {"distance_A": 3.0, "donor_h_acceptor_angle_deg": 89}),
            "weak_hydrogen_bond": ({"distance_A": 3.4, "donor_h_acceptor_angle_deg": 140, "substituent_angle_deg": 90}, {"distance_A": 3.4, "donor_h_acceptor_angle_deg": 129, "substituent_angle_deg": 90}),
            "ionic": ({"distance_A": 3.8}, {"distance_A": 4.01}),
            "hydrophobic": ({"distance_A": 4.4}, {"distance_A": 4.51}),
            "aromatic_atom_contact": ({"distance_A": 3.9}, {"distance_A": 4.01}),
            "aromatic_stacking": ({"distance_A": 5.5, "plane_angle_deg": 20, "centroid_vector_angle_deg": 30}, {"distance_A": 5.5, "plane_angle_deg": 20, "centroid_vector_angle_deg": 91}),
            "cation_pi": ({"distance_A": 4.0, "centroid_vector_angle_deg": 20}, {"distance_A": 4.0, "centroid_vector_angle_deg": 31}),
            "donor_pi": ({"distance_A": 4.0, "centroid_vector_angle_deg": 20}, {"distance_A": 4.51, "centroid_vector_angle_deg": 20}),
            "sulfur_pi": ({"distance_A": 5.9}, {"distance_A": 6.01}),
            "metal_coordination": ({"distance_A": 2.7}, {"distance_A": 2.81}),
            "covalent": ({"distance_A": 1.8, "covalent_bond": True}, {"distance_A": 1.8, "covalent_bond": False}),
            "disulfide": ({"distance_A": 2.0, "covalent_bond": True, "cysteine_sg_pair": True}, {"distance_A": 2.0, "covalent_bond": True, "cysteine_sg_pair": False}),
        }
        for kind, (positive, near_miss) in cases.items():
            with self.subTest(kind=kind):
                self.assertTrue(geometry_rule_accepts(kind, chemistry_supported=True, **positive))
                self.assertFalse(geometry_rule_accepts(kind, chemistry_supported=True, **near_miss))
                self.assertFalse(geometry_rule_accepts(kind, chemistry_supported=False, **positive))

    def test_atom_mapping_is_independent_of_input_atom_order(self):
        first, second = atom("ASP", "A", 11, "OD1"), atom("ARG", "A", 12, "NH1")
        contacts = [
            {"type": "atom-atom", "bgn": first, "end": second, "distance": 3.0, "contact": ["ionic"]},
            {"type": "atom-atom", "bgn": second, "end": first, "distance": 3.0, "contact": ["ionic"]},
        ]
        rows = normalize_arpeggio_contacts(contacts)
        self.assertEqual(rows[0]["interaction_type"], rows[1]["interaction_type"])
        self.assertEqual({rows[0]["src_atom"], rows[0]["dst_atom"]}, {"OD1", "NH1"})

    def test_altloc_ambiguity_and_missing_atom_are_explicit(self):
        endpoint = _endpoint(atom("HIS", "A", 10, "NE2"), {
            ("A", "10", "", "HIS", "NE2"): {"element": "N", "altloc": "ambiguous:A|B"},
        })
        self.assertEqual(endpoint["altloc"], "ambiguous:A|B")
        rows = normalize_arpeggio_contacts([{
            "type": "atom-atom", "bgn": atom("HIS", "A", 10, "NE2"),
            "end": atom("ASP", "A", 20, "OD1"), "distance": 3.1, "contact": [],
        }])
        self.assertEqual(rows, [])

    def test_eligibility_denominators_require_observed_type_atoms(self):
        phe = {"resname": "PHE", "atom_names": {"N", "CA", "C", "O", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"}}
        lys = {"resname": "LYS", "atom_names": {"N", "CA", "C", "O", "NZ", "CB"}}
        complete = interaction_type_eligibility(phe, lys)
        self.assertTrue(complete["cation_pi"])
        self.assertTrue(complete["donor_pi"])
        self.assertFalse(complete["metal_coordination"])
        missing_ring_atom = dict(phe, atom_names=phe["atom_names"] - {"CZ"})
        incomplete = interaction_type_eligibility(missing_ring_atom, lys)
        self.assertFalse(incomplete["aromatic_stacking"])
        self.assertFalse(incomplete["cation_pi"])


if __name__ == "__main__":
    unittest.main()
