"""load_metal_ligands: the non-protein atoms bound to metals, reported per structure (descriptive, never voted)."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from modules import conserved_interactions as ci

ROOT = Path(__file__).resolve().parents[1]
CASES = {"1ag6": ROOT / "reference_structures/0_Plastocyanin/1ag6.cif", "1a6m": ROOT / "reference_structures/3_Myoglobin/1a6m.pdb",
         "1ca2": ROOT / "reference_structures/3_Metals/1ca2.pdb"}
METALS = {"CU", "FE", "ZN", "MG", "MN", "CA", "NA", "K", "NI", "CO", "CD", "HG"}


def line(serial, name, alt, resname, chain, number, xyz, element, record="HETATM", occupancy=1.0):
    return (f"{record:<6s}{serial:5d} {name:<4s}{alt}{resname:>3s} {chain}{number:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"{occupancy:6.2f}{10.0:6.2f}          {element:>2s}")


def brute_force(path):
    """{(metal chain, resname, number, atom, ligand resname, number, atom)}: every hetero atom within 2.8 A of a metal atom of another residue."""
    from Bio.PDB import MMCIFParser, PDBParser

    parser = MMCIFParser(QUIET=True) if str(path).endswith(".cif") else PDBParser(QUIET=True)
    model = next(parser.get_structure("x", str(path)).get_models())
    atoms = [(chain.id, residue, atom) for chain in model for residue in chain if residue.id[0] != " " for atom in residue if atom.element.upper() not in ("H", "D")]
    found = set()
    for chain, residue, atom in atoms:
        if atom.element.upper() not in METALS:
            continue
        for other_chain, other, partner in atoms:
            if other is residue:
                continue
            if float(np.linalg.norm(atom.coord - partner.coord)) <= 2.8:
                found.add((chain, residue.resname.strip(), int(residue.id[1]), atom.get_id(), other.resname.strip(), int(other.id[1]), partner.get_id()))
    return found


class ReferenceStructureTests(unittest.TestCase):
    def test_the_rows_equal_an_independent_recomputation(self):
        for sid, path in CASES.items():
            table = ci.load_metal_ligands(path, structure_id=sid)
            found = {(r.metal_chain, r.metal_resname, r.metal_number, r.metal_atom, r.ligand_resname, r.ligand_number, r.ligand_atom) for r in table.itertuples()}
            self.assertEqual(found, brute_force(path), sid)
            self.assertTrue((table["structure_id"] == sid).all())
            self.assertTrue((table["distance"] <= 2.8).all())

    def test_the_known_ligands(self):
        zinc = ci.load_metal_ligands(CASES["1ca2"])
        self.assertEqual(len(zinc), 1)                                          # the three histidines are protein, so only the water is listed
        row = zinc.iloc[0]
        self.assertEqual((row["metal_resname"], row["ligand_resname"], row["ligand_number"], bool(row["water"])), ("ZN", "HOH", 263, True))
        self.assertAlmostEqual(row["distance"], 2.087, places=3)
        iron = ci.load_metal_ligands(CASES["1a6m"])
        self.assertEqual(set(iron["ligand_resname"]), {"OXY"})                  # the porphyrin nitrogens are the heme's own residue
        self.assertEqual(sorted(iron["ligand_atom"]), ["O1", "O2"])
        self.assertAlmostEqual(float(iron.loc[iron["ligand_atom"] == "O1", "distance"].iloc[0]), 1.806, places=3)
        self.assertFalse(iron["water"].any())
        copper = ci.load_metal_ligands(CASES["1ag6"])
        self.assertTrue(copper.empty)                                           # plastocyanin copper has only protein ligands
        self.assertEqual(list(copper.columns), list(zinc.columns))

    def test_the_reach_is_a_parameter_and_the_legend_says_so(self):
        self.assertTrue(ci.load_metal_ligands(CASES["1ca2"], max_distance=2.0).empty)
        self.assertEqual(len(ci.load_metal_ligands(CASES["1ca2"], max_distance=2.1)), 1)
        rules = ci.rules_of(ci.load_metal_ligands(CASES["1ca2"], max_distance=3.1)).set_index("name")
        self.assertEqual(float(rules.loc["ligand", "cutoff_A"]), 3.1)
        self.assertEqual(rules.loc["ligand", "status"], "opt-in")
        self.assertIn("never voted", rules.loc["ligand", "rule"])


class SyntheticRuleTests(unittest.TestCase):
    def table(self, lines):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "s.pdb"
            path.write_text("\n".join(lines) + "\nEND\n")
            return ci.load_metal_ligands(path, structure_id="s")

    def test_each_rule(self):
        lines = [
            line(1, "ZN", " ", "ZN", "A", 101, (0, 0, 0), "ZN"),
            line(2, "O", "A", "HOH", "A", 201, (2.0, 0, 0), "O", occupancy=0.4),             # first conformer, 2.0 A
            line(3, "O", "B", "HOH", "A", 201, (3.5, 0, 0), "O", occupancy=0.6),             # higher occupancy but beyond reach: ignored
            line(4, "O1", "B", "PER", "A", 202, (0, 2.3, 0), "O"),                           # only a B conformer: kept
            line(5, "OQ1", " ", "KCX", "A", 17, (0, 0, 2.1), "O"),                           # a modified amino acid: a hetero residue, so a ligand
            line(6, "H1", " ", "HOH", "A", 203, (0, 0, -1.0), "H"),                          # hydrogens are not loaded
            line(7, "ND1", " ", "HIS", "A", 40, (-2.0, 0, 0), "N", record="ATOM"),            # a protein atom: never a ligand
            line(8, "FE", " ", "HEM", "A", 301, (10, 0, 0), "FE"),
            line(9, "NA", " ", "HEM", "A", 301, (10, 2.0, 0), "N"),                           # the metal's own residue: not a ligand
            line(10, "O", " ", "HOH", "B", 204, (10, 0, 2.1), "O"),                          # a water on the iron, in another chain
            line(11, "ZN", " ", "ZN", "A", 102, (0, 0, 2.7), "ZN"),                          # a second ion 2.7 A from the first: a ligand atom (its own residue differs)
        ]
        table = self.table(lines)
        zinc = table[(table["metal_resname"] == "ZN") & (table["metal_number"] == 101)].set_index("ligand_atom")
        self.assertAlmostEqual(zinc.loc["O", "distance"], 2.0, places=3)                    # HOH 201: the A conformer
        self.assertAlmostEqual(zinc.loc["O1", "distance"], 2.3, places=3)                   # PER 202: the only conformer
        self.assertEqual(zinc.loc["OQ1", "ligand_resname"], "KCX")
        self.assertEqual(set(zinc["ligand_resname"]), {"HOH", "PER", "KCX", "ZN"})          # no HIS, no hydrogen
        self.assertEqual(bool(zinc.loc["O", "water"]), True)
        self.assertEqual(bool(zinc.loc["O1", "water"]), False)
        iron = table[table["metal_resname"] == "HEM"]
        self.assertEqual(list(iron["ligand_resname"]), ["HOH"])
        self.assertEqual(list(iron["ligand_chain"]), ["B"])
        self.assertEqual(table.equals(table.sort_values(["metal_chain", "metal_number", "metal_atom", "distance"], ignore_index=True)), True)

    def test_the_default_reach_is_the_documented_2_8_angstrom(self):
        rules = ci.rules_of(self.table([line(1, "ZN", " ", "ZN", "A", 101, (0, 0, 0), "ZN")])).set_index("name")
        self.assertEqual(float(rules.loc["ligand", "cutoff_A"]), ci.METAL_LIGAND_DISTANCE)          # the legend follows the default
        limit = ci.METAL_LIGAND_DISTANCE
        self.assertEqual(limit, 2.8)
        inside = self.table([line(1, "ZN", " ", "ZN", "A", 101, (0, 0, 0), "ZN"), line(2, "O", " ", "HOH", "A", 201, (limit - 0.01, 0, 0), "O")])
        outside = self.table([line(1, "ZN", " ", "ZN", "A", 101, (0, 0, 0), "ZN"), line(2, "O", " ", "HOH", "A", 201, (limit + 0.01, 0, 0), "O")])
        self.assertEqual((len(inside), len(outside)), (1, 0))

    def test_a_structure_without_metals_or_ligands_gives_an_empty_table_with_the_columns(self):
        for lines in ([line(1, "ND1", " ", "HIS", "A", 40, (0, 0, 0), "N", record="ATOM")], [line(1, "ZN", " ", "ZN", "A", 101, (0, 0, 0), "ZN")]):
            table = self.table(lines)
            self.assertTrue(table.empty)
            self.assertIn("distance", table.columns)


if __name__ == "__main__":
    unittest.main()
