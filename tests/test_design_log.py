"""The log folder: tables, manifest with provenance and legends, summary. It describes; it carries no directives."""

import json
import re
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from modules import conserved_interactions as ci
from modules import design_log

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = re.compile(r"important|critical|essential|functional|stabiliz|drives|predict|must be retained|should be kept", re.IGNORECASE)


def dataset():
    rows, meta = [], []
    for s in range(12):
        sid = f"s{s}"
        meta.append((sid, "g1" if s < 6 else "g2", f"q{s}"))
        for p in range(1, 15):
            rows.append((sid, p, p, "LEU", "CB", "C", 1.9 * p, 3.6 if p % 2 else 0.0, 0.0))
    atoms = pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])
    atoms.loc[(atoms["structure_id"] == "s3") & (atoms["position"] == 5), "y"] += 9.0
    return atoms, pd.DataFrame(meta, columns=["structure_id", "group", "sequence_id"])


class LogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        atoms, cls.structures = dataset()
        cls.contacts = ci.conserved_contacts(atoms, cls.structures, coverage=True)
        cls.ligands = ci.load_metal_ligands(ROOT / "reference_structures/3_Metals/1ca2.pdb", structure_id="1ca2")
        cls.parameters = {"support": 0.9, "min_sequences": 5, "cutoff": 4.0, "coverage": True, "nan": float("nan"), "numpy": np.int64(3)}

    def write(self, directory, **kwargs):
        return design_log.write_log(directory, {"conserved_contacts": self.contacts, "metal_ligands": self.ligands}, structures=self.structures,
                                    title="Test log", parameters=self.parameters, created="2026-11-02", **kwargs)

    def test_the_folder_has_the_tables_the_manifest_and_the_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = self.write(Path(tmp) / "log")
            self.assertEqual(sorted(p.name for p in out.iterdir()), ["SUMMARY.md", "conserved_contacts.csv", "manifest.json", "metal_ligands.csv", "structures.csv"])
            manifest = json.loads((out / "manifest.json").read_text())
            self.assertEqual((manifest["log_format"], manifest["log_version"], manifest["created"]), ("cni-log", 1, "2026-11-02"))
            self.assertEqual((manifest["dataset"]["structures"], manifest["dataset"]["sequences"]), (12, 12))
            self.assertEqual(manifest["dataset"]["groups"]["g1"], {"sequences": 6, "structures": 6, "tolerated_negatives": 0, "meets_minimum": True})
            self.assertEqual(manifest["tables"]["conserved_contacts"]["rows"], len(self.contacts))
            self.assertEqual(manifest["tables"]["metal_ligands"]["rows"], len(self.ligands))
            self.assertEqual(manifest["parameters"]["nan"], None)                  # NaN is not valid JSON
            self.assertEqual(manifest["parameters"]["numpy"], 3)
            self.assertEqual(len(pd.read_csv(out / "structures.csv")), 12)

    def test_every_column_is_listed_and_the_legend_travels_with_the_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = json.loads((self.write(Path(tmp) / "log") / "manifest.json").read_text())
        info = manifest["tables"]["conserved_contacts"]
        self.assertEqual([c["name"] for c in info["columns"]], list(self.contacts.columns))
        by_name = {c["name"]: c for c in info["columns"]}
        self.assertIn("rule", by_name["conserved_in"])                                  # exact legend name
        self.assertIn("rule", by_name["spare_negatives_g1"])                            # a <group> placeholder in the legend matches any group
        self.assertIn("held_back_by", {r["name"] for r in info["rules"]})
        self.assertTrue(any(r["name"] == "ligand" for r in manifest["tables"]["metal_ligands"]["rules"]))
        self.assertTrue(info["what"])

    def test_the_csv_round_trips(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = self.write(Path(tmp) / "log")
            back = pd.read_csv(out / "conserved_contacts.csv", keep_default_na=False)
        self.assertEqual(len(back), len(self.contacts))
        self.assertEqual(back["conserved_contact"].astype(bool).tolist(), self.contacts["conserved_contact"].tolist())
        self.assertEqual(back["conserved_in"].astype(str).tolist(), self.contacts["conserved_in"].fillna("").astype(str).tolist())

    def test_writing_twice_gives_identical_files_and_an_existing_log_is_protected(self):
        with tempfile.TemporaryDirectory() as tmp:
            first, second = self.write(Path(tmp) / "a"), self.write(Path(tmp) / "b")
            for name in ("manifest.json", "SUMMARY.md", "conserved_contacts.csv", "metal_ligands.csv", "structures.csv"):
                self.assertEqual((first / name).read_bytes(), (second / name).read_bytes(), name)
            with self.assertRaises(FileExistsError):
                self.write(first)
            (first / "mine.txt").write_text("keep")
            self.write(first, overwrite=True)
            self.assertEqual((first / "mine.txt").read_text(), "keep")              # other files are never touched

    def test_bad_inputs_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            for bad in ("has space", "manifest", "structures", "a/b"):
                with self.assertRaises(ValueError):
                    design_log.write_log(Path(tmp) / "x", {bad: self.ligands}, structures=self.structures, title="t")
            with self.assertRaises(ValueError):
                design_log.write_log(Path(tmp) / "x", {}, structures=self.structures[["structure_id"]], title="t")

    def test_the_summary_states_the_meaning_and_the_depth_of_each_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            text = (self.write(Path(tmp) / "log") / "SUMMARY.md").read_text()
            small = (design_log.write_log(Path(tmp) / "small", {}, structures=self.structures, title="t", parameters={"min_sequences": 7}) / "SUMMARY.md").read_text()
        self.assertIn(design_log.LOG_MEANING, text)
        self.assertIn("Group g1 has 6 sequence(s) and tolerates no negative sequence", text)
        self.assertIn("Group g1 has fewer than the minimum of 7 voting sequences", small)
        self.assertIn("## Rules in force: conserved_contacts", text)
        self.assertIn("HOH", text)                                                      # the metal-ligand section lists the zinc water

    def test_nothing_outside_the_meaning_statement_claims_importance_or_prediction(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = self.write(Path(tmp) / "log")
            summary = (out / "SUMMARY.md").read_text().replace(design_log.LOG_MEANING, "")
            manifest = json.loads((out / "manifest.json").read_text())
            manifest.pop("meaning")
            self.assertIsNone(FORBIDDEN.search(summary), FORBIDDEN.search(summary))
            self.assertIsNone(FORBIDDEN.search(json.dumps(manifest)), FORBIDDEN.search(json.dumps(manifest)))
            self.assertEqual([p.name for p in out.iterdir() if "directive" in p.name.lower()], [])      # the log writes no directives
            self.assertNotIn("directive", json.dumps(manifest).lower())


if __name__ == "__main__":
    unittest.main()
