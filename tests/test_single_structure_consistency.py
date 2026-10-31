"""Scope check 4: the CNI functions and ``sscna analyze`` must agree on one structure wherever their rules overlap.

Three reference structures with different cofactors: spinach plastocyanin (Cu), sperm whale myoglobin (heme), human carbonic
anhydrase II (Zn). The differences that remain are pinned as known and explained, so a change in either side is noticed.
"""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from modules import analyze_cli
from modules import conserved_interactions as ci
from modules import coordination_table
from modules.coordination_api import analyze_structure

ROOT = Path(__file__).resolve().parents[1]
CASES = {"1ag6": (ROOT / "reference_structures/0_Plastocyanin/1ag6.cif", "CU", "A", ("CU",)),
         "1a6m": (ROOT / "reference_structures/3_Myoglobin/1a6m.pdb", "HEM", "A", ("FE",)),
         "1ca2": (ROOT / "reference_structures/3_Metals/1ca2.pdb", "ZN", "A", ("ZN",))}
STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def cofactor_coordinates(path, resname):
    """Heavy atoms of the cofactor residue(s) as (name, element, xyz) by an independent parse."""
    from Bio.PDB import MMCIFParser, PDBParser

    parser = MMCIFParser(QUIET=True) if str(path).endswith(".cif") else PDBParser(QUIET=True)
    model = next(parser.get_structure("x", str(path)).get_models())
    return [(a.get_id(), a.element.upper(), np.array(a.coord, float)) for r in model.get_residues() if r.resname == resname
            for a in r if a.element.upper() not in ("H", "D")]


def features_for(sid, atoms):
    rows = [(sid, "cofactor", name, element, "metal" if element in ("CU", "FE", "ZN") else "oxygen", *xyz) for name, element, xyz in atoms]
    return pd.DataFrame(rows, columns=["structure_id", "feature", "atom", "element", "kind", "x", "y", "z"])


class CommandLineTests(unittest.TestCase):
    def test_the_command_writes_exactly_the_tables_the_api_returns(self):
        path, cofactor, _, _ = CASES["1ag6"]
        result = analyze_structure(path, cofactor, shells=3, direct_coordination=True, first_model_only=True)
        with tempfile.TemporaryDirectory() as out:
            self.assertEqual(analyze_cli.main(["--input", str(path), "--cofactor", cofactor, "--no-plots", "--output-dir", out,
                                               "--direct-coordination", "--first-model"]), 0)
            for table, suffix in (("residues", "Coordination_Residues"), ("atoms", "Coordination_Atoms"), ("links", "Coord_Links"), ("contacts", "Coord_Contacts")):
                written = pd.read_csv(Path(out) / f"{path.name}_{suffix}.csv")
                expected = result[table]
                self.assertEqual(list(written.columns), list(expected.columns), table)
                self.assertEqual(len(written), len(expected), table)
                for column in expected.columns:
                    numeric = pd.to_numeric(expected[column], errors="coerce")
                    if expected[column].dtype.kind == "f" or (expected[column].dtype == object and numeric.notna().all() and column.endswith("_A")):
                        # the links table keeps distances as three-decimal strings ('2.010', a pinned legacy format); compare the values
                        np.testing.assert_allclose(pd.to_numeric(written[column]).to_numpy(float), numeric.to_numpy(float), atol=1e-9, err_msg=f"{table}.{column}")
                    else:
                        self.assertEqual(written[column].fillna("").astype(str).tolist(), expected[column].fillna("").astype(str).tolist(), f"{table}.{column}")

    def test_the_command_defaults_differ_from_the_api_in_the_documented_way(self):
        self.assertEqual((analyze_cli.DEFAULT_SHELLS, analyze_cli.DEFAULT_DIRECT_COORDINATION_CUTOFF), (3, 2.6))
        path, cofactor, _, _ = CASES["1ag6"]
        shells = analyze_structure(path, cofactor, first_model_only=True)["residues"]["shell"].unique()
        self.assertNotIn("TCS", shells)                                      # the API default is two shells, the command's is three


class OverlapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = {}
        for sid, (path, cofactor, chain, metals) in CASES.items():
            analysis = analyze_structure(path, cofactor, shells=3, direct_coordination=True, first_model_only=True)
            atoms = ci.load_atoms(path, chain, structure_id=sid)
            coordinates = cofactor_coordinates(path, cofactor)
            cls.data[sid] = {"analysis": analysis, "atoms": atoms, "features": features_for(sid, coordinates), "cofactor": coordinates, "metals": metals}

    def protein(self, residues):
        return residues[residues["residue_name"].isin(STANDARD)]

    def test_the_primary_shell_is_every_residue_with_an_n_o_or_s_atom_within_3_6_angstrom_of_any_cofactor_atom(self):
        for sid, d in self.data.items():
            residues = d["analysis"]["residues"]
            shell = {(int(r.residue_number), r.residue_name) for r in residues[residues["shell"] == "PCS"].itertuples()}
            _, cofactor_name, chain, _ = CASES[sid]
            from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
            parser = MMCIFParser(QUIET=True) if str(CASES[sid][0]).endswith(".cif") else PDBParser(QUIET=True)
            model = next(parser.get_structure("x", str(CASES[sid][0])).get_models())
            heavy = [a for a in model.get_atoms() if a.element.upper() not in ("H", "D")]
            search = NeighborSearch(heavy)
            expected = set()
            for _, _, xyz in d["cofactor"]:
                for atom in search.search(xyz, 3.6):
                    residue = atom.get_parent()
                    if residue.resname != cofactor_name and atom.element.upper() in ("N", "O", "S"):
                        expected.add((int(residue.id[1]), residue.resname))
            self.assertEqual(shell, expected, sid)

    def test_cnis_cofactor_contacts_cover_the_protein_primary_shell_with_matching_distances(self):
        for sid, d in self.data.items():
            contacts = ci.cofactor_contacts(d["atoms"], d["features"])
            recorded = contacts.set_index("position")
            analysis = d["analysis"]["contacts"]
            first = analysis[(analysis["src_shell"] == "Cofactor") & analysis["dst_resname"].isin(STANDARD)]
            self.assertTrue(len(first) > 0, sid)
            for residue, group in first.groupby("dst_resnum"):
                self.assertIn(residue, recorded.index, sid)
                self.assertTrue(bool(recorded.loc[residue, "contact"]), (sid, residue))                      # within 3.6 A implies within CNI's 4.0 A contact
                self.assertLessEqual(float(recorded.loc[residue, "min_distance"]), float(group["distance_A"].min()) + 1e-6, (sid, residue))   # CNI's minimum includes carbons
            xyz = d["atoms"].set_index(["position", "atom"])[["x", "y", "z"]]
            by_name = {name: coord for name, _, coord in d["cofactor"]}
            for row in first.itertuples():                                                                     # the loader's coordinates give the analysis's distance
                distance = float(np.linalg.norm(xyz.loc[(row.dst_resnum, row.dst_atom)].to_numpy(float) - by_name[row.src_atom]))
                self.assertAlmostEqual(distance, row.distance_A, places=5, msg=(sid, row.dst_resnum, row.dst_atom))

    def test_direct_coordination_and_coord_agree_at_the_same_cutoff_for_protein_residues(self):
        for sid, d in self.data.items():
            for cutoff in (2.6, 2.8):
                analysis = analyze_structure(CASES[sid][0], CASES[sid][1], shells=3, direct_coordination=True, direct_coordination_cutoff=cutoff, first_model_only=True)["contacts"]
                direct = {int(r.dst_resnum) for r in analysis[analysis["direct_coordination"] & analysis["dst_resname"].isin(STANDARD)].itertuples()}
                metal_features = d["features"][d["features"]["element"].isin(d["metals"])]
                coord = ci.cofactor_contacts(d["atoms"], metal_features, coordination_distance=cutoff)
                self.assertEqual({int(p) for p in coord.loc[coord["coord"], "position"]}, direct, (sid, cutoff))

    def test_the_tiers_agree_between_the_two_routes(self):
        d = self.data["1ag6"]
        route_one = coordination_table.annotate_contacts(d["analysis"]["contacts"])
        first = route_one[(route_one["src_shell"] == "Cofactor") & route_one["dst_resname"].isin(STANDARD)]
        metal = d["features"][d["features"]["element"] == "CU"]
        route_two = ci.cofactor_contacts(d["atoms"], metal, coordination_model="by_pair").set_index("position")
        for residue, group in first.groupby("dst_resnum"):
            tiers = {t for t in group["coord_tier"] if t in ("bond", "weak")}
            expected = route_two.loc[residue, "coord_tier"]
            self.assertEqual(tiers, {expected} if expected in ("bond", "weak") else set(), residue)

    def test_the_known_differences_are_what_we_say_they_are(self):
        # 1. waters and ligands: analyze_structure counts them (the zinc-bound water of carbonic anhydrase is a direct ligand); the loader keeps standard residues only.
        zinc = self.data["1ca2"]["analysis"]["contacts"]
        direct = zinc[zinc["direct_coordination"]]
        self.assertEqual(set(direct["dst_resname"]), {"HIS", "HOH"})
        self.assertEqual(set(self.data["1ca2"]["atoms"]["resname"]) - STANDARD, set())
        # 2. the oxygen ligand of myoglobin is a direct ligand too, and not a protein residue
        iron = self.data["1a6m"]["analysis"]["contacts"]
        self.assertIn("OXY", set(iron.loc[iron["direct_coordination"], "dst_resname"]))
        # 3. the two cutoffs: the command's 2.6 A and CNI's 2.8 A; the axial Met92 of plastocyanin (2.88 A) is in neither, and only the tiers name it
        plastocyanin = self.data["1ag6"]
        at_default = ci.cofactor_contacts(plastocyanin["atoms"], plastocyanin["features"])
        self.assertNotIn(92, set(at_default.loc[at_default["coord"], "position"]))
        tiered = ci.cofactor_contacts(plastocyanin["atoms"], plastocyanin["features"], coordination_model="by_pair").set_index("position")
        self.assertEqual(tiered.loc[92, "coord_tier"], "weak")
        # 4. shell membership is by moiety atoms for non-metal cofactors: the heme primary shell includes residues touching only its methyl carbons
        heme = self.data["1a6m"]["analysis"]["contacts"]
        carbon_sourced = heme[(heme["src_shell"] == "Cofactor") & (heme["src_element"] == "C")]
        self.assertGreater(len(carbon_sourced), 0)


if __name__ == "__main__":
    unittest.main()
