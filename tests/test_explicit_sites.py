"""Tests for caller-defined cofactor sites."""

from pathlib import Path
import unittest

from modules.coordination_api import analyze_structure

ROOT = Path(__file__).resolve().parents[1]
FERREDOXIN = ROOT / "reference_structures/3_Metals/2zvs.pdb"
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"


def sf4_keys():
    tables = analyze_structure(FERREDOXIN, "SF4", site_mode="per-site")
    cofactors = tables["atoms"][tables["atoms"]["shell"] == "Cofactor"]
    return sorted({(r.residue_name, r.residue_number, r.chain, r.insertion_code, r.hetero_flag)
                   for r in cofactors.itertuples()})


class ExplicitSiteTests(unittest.TestCase):
    def test_named_sites_group_chosen_residues_in_given_order(self):
        keys = sf4_keys()
        self.assertEqual(len(keys), 6)
        sites = {"pair": [keys[0], keys[1]], "single": [keys[2]]}
        tables = analyze_structure(FERREDOXIN, "SF4", site_mode="per-site", explicit_sites=sites, shells=3)
        atoms = tables["atoms"]
        self.assertEqual(list(dict.fromkeys(atoms["site_id"])), ["pair", "single"])
        cofactor = atoms[atoms["shell"] == "Cofactor"]
        self.assertEqual(cofactor[cofactor["site_id"] == "pair"]["residue_number"].nunique(), 2)
        self.assertEqual(cofactor[cofactor["site_id"] == "single"]["residue_number"].nunique(), 1)

    def test_single_residue_site_matches_default_per_site_analysis(self):
        keys = sf4_keys()
        default = analyze_structure(FERREDOXIN, "SF4", site_mode="per-site", shells=3)
        explicit = analyze_structure(FERREDOXIN, "SF4", site_mode="per-site", shells=3,
                                     explicit_sites={"site_1": [keys[0]]})
        first = default["atoms"][default["atoms"]["site_id"] == "site_1"].reset_index(drop=True)
        self.assertTrue(first.equals(explicit["atoms"].reset_index(drop=True)))

    def test_invalid_use_is_rejected(self):
        keys = sf4_keys()
        with self.assertRaises(ValueError):
            analyze_structure(FERREDOXIN, "SF4", site_mode="union", explicit_sites={"a": [keys[0]]})
        with self.assertRaises(ValueError):
            analyze_structure(FERREDOXIN, "SF4", site_mode="per-site", combinatorial=True,
                              explicit_sites={"a": [keys[0]]})
        with self.assertRaises(ValueError):
            analyze_structure(FERREDOXIN, "SF4", site_mode="per-site", explicit_sites={"a": [("SF4", 999, "Z", "", "H_SF4")]})
        with self.assertRaises(ValueError):
            analyze_structure(PLASTOCYANIN, "CU", site_mode="per-site", explicit_sites={"a": []})

    def test_default_behavior_is_unchanged_without_the_option(self):
        tables = analyze_structure(FERREDOXIN, "SF4", site_mode="per-site")
        self.assertEqual(tables["atoms"]["site_id"].nunique(), 6)


if __name__ == "__main__":
    unittest.main()
