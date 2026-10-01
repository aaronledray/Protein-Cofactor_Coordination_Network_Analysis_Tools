"""Tests for insertion-code and model-selection handling."""

from contextlib import redirect_stdout
from io import StringIO
import unittest

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from modules.structure_processing import identify_coordination_network


def atom(name, coordinate, element, serial):
    return Atom(
        name,
        np.asarray(coordinate, dtype=float),
        1.0,
        1.0,
        " ",
        f"{name:>4}",
        serial,
        element=element,
    )


def structure_with_insertions():
    structure = Structure("insertions")
    model = Model(0)
    chain = Chain("A")
    cofactor = Residue(("H_CU", 200, " "), "CU", " ")
    cofactor.add(atom("CU", (0.0, 0.0, 0.0), "CU", 1))
    chain.add(cofactor)

    residue_52 = Residue((" ", 52, " "), "HIS", " ")
    residue_52.add(atom("ND1", (2.0, 0.0, 0.0), "N", 2))
    chain.add(residue_52)
    residue_52a = Residue((" ", 52, "A"), "HIS", " ")
    residue_52a.add(atom("ND1", (2.5, 0.0, 0.0), "N", 3))
    chain.add(residue_52a)
    model.add(chain)
    structure.add(model)
    return structure


def structure_with_two_models():
    structure = Structure("models")
    for model_id, offset in ((0, 0.0), (1, 10.0)):
        model = Model(model_id)
        chain = Chain("A")
        cofactor = Residue(("H_CU", 200, " "), "CU", " ")
        cofactor.add(atom("CU", (offset, 0.0, 0.0), "CU", model_id * 3 + 1))
        chain.add(cofactor)
        residue = Residue((" ", 52, " "), "HIS", " ")
        residue.add(atom("ND1", (offset + 2.0, 0.0, 0.0), "N", model_id * 3 + 2))
        chain.add(residue)
        model.add(chain)
        structure.add(model)
    return structure


class StructureIdentityTests(unittest.TestCase):
    def test_insertion_codes_keep_same_numbered_residues_distinct(self):
        with redirect_stdout(StringIO()):
            cofactor, pcs, scs = identify_coordination_network(
                structure_with_insertions(),
                ["CU"],
                distance_cutoff=3.6,
                expand_residues=False,
                combinatorial_mode=False,
                write_coord_links=False,
            )

        self.assertEqual(len(cofactor), 1)
        self.assertEqual(len(pcs), 2)
        self.assertEqual(
            {atom["insertion_code"] for atom in pcs},
            {"", "A"},
        )

    def test_first_model_option_limits_model_iteration(self):
        with redirect_stdout(StringIO()):
            all_models = identify_coordination_network(
                structure_with_two_models(),
                ["CU"],
                distance_cutoff=3.6,
                expand_residues=False,
                combinatorial_mode=False,
                write_coord_links=False,
            )
            first_model = identify_coordination_network(
                structure_with_two_models(),
                ["CU"],
                distance_cutoff=3.6,
                expand_residues=False,
                combinatorial_mode=False,
                write_coord_links=False,
                first_model_only=True,
            )

        self.assertEqual(len(all_models[0]), 2)
        self.assertEqual(len(first_model[0]), 1)
        self.assertEqual(len(first_model[1]), 1)


if __name__ == "__main__":
    unittest.main()
