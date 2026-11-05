import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from modules import conserve_cli, conserved_interactions as ci, state_stratification

ROOT = Path(__file__).resolve().parents[1]

ATOM_COLUMNS = ["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"]


def fixture():
    structures, atoms = [], []
    for seq in range(6):
        for state in ("ox", "red"):
            for replicate in range(2):
                sid = f"q{seq}_{state}_{replicate}"
                structures.append((sid, "A", f"q{seq}", state))
                has_b = not (state == "red" and seq == 5)
                contact = state == "ox" or (state == "red" and (seq < 3 or (seq == 3 and replicate == 0)))
                atoms.append((sid, 10, 1, "ALA", "CA", "C", 0.0, 0.0, 0.0))
                if has_b:
                    atoms.append((sid, 20, 5, "VAL", "CA", "C", 3.0 if contact else 5.0, 0.0, 0.0))
    for seq in range(6):
        sid = f"q{seq}_apo_0"
        structures.append((sid, "A", f"q{seq}", "apo"))
        atoms.append((sid, 10, 1, "ALA", "CA", "C", 0.0, 0.0, 0.0))
    return (pd.DataFrame(atoms, columns=ATOM_COLUMNS),
            pd.DataFrame(structures, columns=["structure_id", "group", "sequence_id", "state"]))


class StateStratificationTests(unittest.TestCase):
    def test_votes_are_independent_and_missing_positions_are_not_negative(self):
        atoms, structures = fixture()
        result = state_stratification.state_stratified_contacts(atoms, structures, support=0.9, min_sequences=5).set_index("state")
        structure_votes = state_stratification.state_structure_votes(atoms, structures).set_index(
            ["position_a", "position_b", "state", "sequence_id"])
        ox, red, apo = result.loc["ox"], result.loc["red"], result.loc["apo"]
        self.assertEqual((ox["contact_n_A"], ox["contact_positive_A"], ox["contact_A"], ox["passes_contact_A"], ox["conserved_contact"]), (6, 6, 1.0, True, True))
        self.assertEqual((red["contact_n_A"], red["contact_positive_A"], red["contact_A"], red["passes_contact_A"], red["conserved_contact"]), (5, 4, 0.8, False, False))
        self.assertEqual((apo["contact_n_A"], apo["conserved_contact"]), (0, False))
        self.assertTrue(pd.isna(apo["contact_A"]))
        self.assertEqual(red["one_negative_away"], "q4")
        self.assertEqual(apo["spare_voters_A"], -5)
        self.assertEqual(set(result.index), {"apo", "ox", "red"})
        self.assertEqual(tuple(structure_votes.loc[(10, 20, "red", "q3"), ["eligible_structures", "positive_structures", "structure_fraction", "sequence_vote"]]), (2, 1, 0.5, True))
        self.assertEqual(tuple(structure_votes.loc[(10, 20, "red", "q5"), ["eligible_structures", "positive_structures"]]), (0, 0))
        self.assertTrue(pd.isna(structure_votes.loc[(10, 20, "red", "q5"), "sequence_vote"]))
        self.assertEqual(tuple(structure_votes.loc[(10, 20, "apo", "q0"), ["eligible_structures", "positive_structures"]]), (0, 0))

    def test_each_state_matches_the_existing_vote_on_its_subset(self):
        atoms, structures = fixture()
        result = state_stratification.state_stratified_contacts(atoms, structures, support=0.9, min_sequences=5)
        ox_ids = set(structures.loc[structures["state"] == "ox", "structure_id"])
        direct = ci.conserved_contacts(atoms[atoms["structure_id"].isin(ox_ids)],
                                       structures[structures["state"] == "ox"].drop(columns="state"),
                                       support=0.9, min_sequences=5)
        state_ox = result[result["state"] == "ox"]
        self.assertEqual(state_ox["conserved_contact"].tolist(), direct["conserved_contact"].tolist())
        self.assertEqual(state_ox["contact_n_A"].tolist(), direct["contact_n_A"].tolist())

    def test_input_contract(self):
        atoms, structures = fixture()
        with self.assertRaisesRegex(ValueError, "at least two"):
            state_stratification.state_stratified_contacts(atoms, structures.assign(state="ox"))
        with self.assertRaisesRegex(ValueError, "non-empty state"):
            state_stratification.state_stratified_contacts(atoms, structures.assign(state=""))

    def test_cli_writes_the_opt_in_table_and_keeps_pooled_table(self):
        structure_file = ROOT / "reference_structures" / "0_Plastocyanin" / "1ag6.cif"
        if not structure_file.is_file():
            self.skipTest("plastocyanin structure fixture is unavailable")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = []
            for i in range(5):
                for state in ("ox", "red"):
                    rows.append((f"s{i}_{state}", str(structure_file), f"q{i}", "all", "A", state))
            manifest = root / "manifest.csv"
            pd.DataFrame(rows, columns=["structure_id", "path", "sequence_id", "group", "chain", "state"]).to_csv(manifest, index=False)
            code = conserve_cli.main(["--manifest", str(manifest), "--output-dir", str(root / "log"), "--stratify-state"])
            self.assertEqual(code, 0)
            log = root / "log"
            pooled = pd.read_csv(log / "conserved_contacts.csv")
            by_state = pd.read_csv(log / "state_contacts.csv")
            structure_votes = pd.read_csv(log / "state_structure_votes.csv")
            self.assertEqual(set(by_state["state"]), {"ox", "red"})
            self.assertEqual(len(by_state), 2 * len(pooled))
            self.assertEqual(len(structure_votes), 10 * len(pooled))
            self.assertTrue((structure_votes["eligible_structures"] == 1).all())
            self.assertTrue((structure_votes["positive_structures"] == 1).all())
            metadata = json.loads((log / "manifest.json").read_text())
            self.assertTrue(metadata["parameters"]["stratify_state"])
            self.assertIn("state_contacts", metadata["tables"])
            self.assertIn("do not test for a state effect", (log / "SUMMARY.md").read_text())
            self.assertEqual(set(pd.read_csv(log / "structures.csv")["state"]), {"ox", "red"})


if __name__ == "__main__":
    unittest.main()
