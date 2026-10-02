"""Tests for hypothetical substrate-point seeds and chains."""

from pathlib import Path
import unittest

import numpy as np
import pandas as pd

from modules.coordination_api import analyze_structure
from modules.substrate_seeds import (
    analyze_substrate_seed,
    resolve_substrate_point,
    substrate_chains,
)

ROOT = Path(__file__).resolve().parents[1]
MYOGLOBIN = ROOT / "reference_structures/3_Myoglobin/1a6m.pdb"
FERREDOXIN = ROOT / "reference_structures/3_Metals/2zvs.pdb"
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


def distal_point(tables):
    """Fe + 2 Å along the proximal-His -> Fe axis (the ligand-binding side)."""
    atoms = tables["atoms"]
    fe = atoms[(atoms["shell"] == "Cofactor") & (atoms["atom_name"] == "FE")][["x", "y", "z"]].to_numpy(float)[0]
    ne2 = atoms[(atoms["residue_name"] == "HIS") & (atoms["residue_number"] == 93)
                & (atoms["atom_name"] == "NE2")][["x", "y", "z"]].to_numpy(float)[0]
    return fe + 2.0 * (fe - ne2) / np.linalg.norm(fe - ne2)


class SubstrateSeedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables = analyze_structure(
            MYOGLOBIN, "HEM", shells=3, site_mode="per-site", direct_coordination=True
        )
        cls.point = distal_point(cls.tables)

    def test_distal_point_chains_run_through_the_distal_histidine(self):
        result = substrate_chains(self.tables, self.point)
        chains = result["chains"]
        self.assertEqual(list(chains["rank"]), sorted(chains["rank"]))
        self.assertTrue((chains["terminal_atom"] == "FE").all())
        distal = chains[chains["entry_residue"] == "HIS64"].iloc[0]
        self.assertEqual(distal["shells_traversed"], "SCS>PCS>Cofactor")
        self.assertEqual(distal["n_hops"], 3)
        steps = result["steps"]
        path = steps[steps["chain_id"] == distal["chain_id"]]
        self.assertEqual(path["step"].tolist(), [0, 1, 2, 3])
        self.assertEqual(path["node_type"].tolist(), ["substrate_point", "network", "network", "cofactor"])
        self.assertEqual(path["edge_from_previous"].iloc[1:].tolist(),
                         ["point_contact", "shell_contact", "direct_coordination"])
        # hops and total length are consistent with the per-edge distances
        self.assertAlmostEqual(path["edge_distance_A"].sum(), distal["total_length_A"], places=2)

    def test_relative_and_explicit_points_agree(self):
        site = self.tables["atoms"]
        offset = tuple(self.point - resolve_substrate_point(site, {"cofactor_atom": "FE", "offset": (0, 0, 0)}))
        relative = substrate_chains(self.tables, {"cofactor_atom": "FE", "offset": offset})
        explicit = substrate_chains(self.tables, self.point)
        pd.testing.assert_frame_equal(relative["chains"], explicit["chains"])

    def test_contacts_are_sorted_and_within_cutoff(self):
        contacts = substrate_chains(self.tables, self.point, contact_cutoff=3.0)["contacts"]
        self.assertTrue((contacts["distance_A"] <= 3.0).all())
        self.assertEqual(list(contacts["distance_A"]), sorted(contacts["distance_A"]))

    def test_far_point_has_no_contacts_or_chains(self):
        result = substrate_chains(self.tables, (500.0, 500.0, 500.0))
        self.assertEqual(len(result["seed"]), 1)
        self.assertTrue(result["contacts"].empty and result["chains"].empty and result["steps"].empty)

    def test_max_chains_and_determinism(self):
        first = substrate_chains(self.tables, self.point, max_chains=2)
        self.assertEqual(len(first["chains"]), 2)
        shuffled = {k: v.sample(frac=1, random_state=4).reset_index(drop=True) for k, v in self.tables.items()}
        again = substrate_chains(shuffled, self.point, max_chains=2)
        pd.testing.assert_frame_equal(first["chains"], again["chains"])

    def test_invalid_inputs_are_rejected(self):
        for point in ((1.0, 2.0), (float("nan"), 0.0, 0.0), {"cofactor_atom": "FE"},
                      {"cofactor_atom": "NOPE", "offset": (0, 0, 1)}):
            with self.assertRaises(ValueError):
                substrate_chains(self.tables, point)
        with self.assertRaises(ValueError):
            substrate_chains(self.tables, self.point, contact_cutoff=0)
        with self.assertRaises(ValueError):
            substrate_chains(self.tables, self.point, max_chains=0)

    def test_each_cofactor_site_resolves_its_own_relative_point(self):
        result = analyze_substrate_seed(
            FERREDOXIN, "SF4", {"cofactor_atom": "FE1", "offset": (0.0, 0.0, 2.5)}
        )
        seeds = result["seed"]
        self.assertEqual(len(seeds), 6)
        self.assertEqual(len(seeds[["x", "y", "z"]].drop_duplicates()), 6)
        self.assertTrue((result["chains"]["terminal_atom"].str.startswith(("FE", "S"))).all())

    def test_wrapper_returns_analysis_tables_plus_seed_outputs(self):
        result = analyze_substrate_seed(PLASTOCYANIN, "CU", {"cofactor_atom": "CU", "offset": (0, 0, 2.2)})
        for key in ("residues", "atoms", "links", "contacts", "seed", "substrate_contacts", "chains", "steps"):
            self.assertIn(key, result)
        self.assertGreater(len(result["chains"]), 0)


if __name__ == "__main__":
    unittest.main()
