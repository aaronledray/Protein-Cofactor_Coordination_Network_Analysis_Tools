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


def specific_dataset():
    """Like dataset(), but the contacts exist only in group g2: in g1 positions i and i+2 are moved 20 A apart (the offset alternates every two positions)."""
    atoms, structures = dataset()
    g1 = atoms["structure_id"].isin(set(structures.loc[structures["group"] == "g1", "structure_id"]))
    atoms.loc[g1, "y"] += 20.0 * ((atoms.loc[g1, "position"] // 2) % 2)
    return atoms, structures


class DirectiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        atoms, cls.structures = dataset()
        cls.contacts = ci.conserved_contacts(atoms, cls.structures, coverage=True)
        cls.directives = design_log.contact_directives(cls.contacts, cls.structures)

    def test_tiers_and_evidence_follow_the_vote(self):
        d = self.directives.set_index("subject")
        self.assertEqual(self.directives["tier"].value_counts().to_dict(), {"conserved_in_all": 10, "conserved_in_subset": 2})
        row = d.loc["contact 1-3"]
        self.assertEqual((row["tier"], row["conserved_in"], row["failing_groups"], row["voters_g1"], row["positives_g1"], row["spare_negatives_g1"]),
                         ("conserved_in_all", "g1+g2", "", 6, 6, 0))
        self.assertEqual(row["chemistry"], "sc_sc;hydrophobic")
        row = d.loc["contact 3-5"]                                       # the sequence that breaks it is in g1: 5 of 6 fails there, g2 passes
        self.assertEqual((row["tier"], row["conserved_in"], row["failing_groups"], row["positives_g1"], row["spare_negatives_g1"]), ("conserved_in_subset", "g2", "g1", 5, -1))
        self.assertEqual(set(self.directives["claim_level"]), {"descriptive"})
        self.assertEqual(self.directives["directive_id"].iloc[0], "retain_contact:1-3")

    def test_a_group_specific_contact_needs_the_others_to_show_it_in_under_half_their_sequences(self):
        atoms, structures = specific_dataset()
        contacts = ci.conserved_contacts(atoms, structures, coverage=True)
        d = design_log.contact_directives(contacts, structures)
        self.assertEqual(d["tier"].value_counts().to_dict(), {"family_specific": 12})
        self.assertEqual(set(d["conserved_in"]), {"g2"})
        self.assertEqual(set(d["positives_g1"]), {0})
        # a looser threshold on "mostly absent" turns the same rows into subset directives only when the others show them often enough
        d = design_log.contact_directives(contacts, structures, specific_below=0.0)
        self.assertEqual(d["tier"].value_counts().to_dict(), {"conserved_in_subset": 12})

    def hand_made(self, g1, g2):
        """A contacts table with one pair; g1 and g2 are (voting sequences, positives)."""
        table = pd.DataFrame([{"position_a": 1, "position_b": 3, "contact_n_g1": g1[0], "contact_g1": round(g1[1] / g1[0], 3) if g1[0] else np.nan,
                               "contact_n_g2": g2[0], "contact_g2": round(g2[1] / g2[0], 3) if g2[0] else np.nan}])
        structures = pd.DataFrame([("s0", "g1", "q0"), ("s1", "g2", "q1")], columns=["structure_id", "group", "sequence_id"])
        return design_log.contact_directives(table, structures)

    def test_the_boundaries_of_the_tiers(self):
        self.assertEqual(self.hand_made((10, 9), (10, 10)).iloc[0]["tier"], "conserved_in_all")                        # 90% is inclusive
        self.assertEqual(self.hand_made((10, 8), (10, 10)).iloc[0]["tier"], "conserved_in_subset")                     # 80%: g1 fails but is not mostly absent
        self.assertEqual(self.hand_made((10, 5), (10, 10)).iloc[0]["tier"], "conserved_in_subset")                     # exactly half is not "fewer than half"
        self.assertEqual(self.hand_made((10, 4), (10, 10)).iloc[0]["tier"], "family_specific")
        self.assertEqual(self.hand_made((10, 0), (10, 10)).iloc[0]["tier"], "family_specific")
        self.assertEqual(len(self.hand_made((4, 4), (4, 4))), 0)                                                       # fewer than five voters: conserved nowhere, no directive
        self.assertEqual(self.hand_made((4, 4), (10, 10)).iloc[0]["tier"], "conserved_in_subset")                      # g1 has voters, but too few to count, and all positive
        self.assertEqual(self.hand_made((0, 0), (10, 10)).iloc[0]["tier"], "conserved_in_subset")                      # no voters in g1: absence is not shown
        self.assertEqual(len(self.hand_made((10, 5), (10, 5))), 0)

    def test_a_group_without_voters_is_named_in_the_statement_and_not_called_absent(self):
        row = self.hand_made((0, 0), (10, 10)).iloc[0]
        self.assertEqual(row["tier"], "conserved_in_subset")
        self.assertIn("g1 0/0", row["statement"])
        self.assertTrue(row["statement"].endswith("; no voting sequences in g1)."))
        self.assertNotIn("no voting sequences", self.hand_made((10, 10), (10, 10)).iloc[0]["statement"])             # only when it matters
        self.assertNotIn("no voting sequences", self.hand_made((10, 4), (10, 10)).iloc[0]["statement"])             # g1 votes (and is mostly absent): no clause

    def test_the_statements_with_a_silent_group_still_match_the_template(self):
        pattern = re.compile(r"^[^:]+: (conserved in all groups|conserved in some groups only|conserved in one group and mostly absent from the others) "
                             r"\(conserved in [\w+ ]+; sequences with it / sequences voting: [^)]*?(; no voting sequences in [\w+ ]+)?\)\.$")
        for g1, g2 in (((0, 0), (10, 10)), ((10, 4), (10, 10)), ((10, 10), (10, 10))):
            self.assertRegex(self.hand_made(g1, g2).iloc[0]["statement"], pattern)

    def test_a_group_with_no_voting_sequences_does_not_make_a_contact_family_specific(self):
        atoms, structures = dataset()
        only = structures[structures["group"] == "g2"]
        contacts = ci.conserved_contacts(atoms[atoms["structure_id"].isin(only["structure_id"])], structures, coverage=True, groups=["g1", "g2"])
        d = design_log.contact_directives(contacts, structures)
        self.assertEqual(set(d["tier"]), {"conserved_in_subset"})        # g1 has no voters: absent from the evidence, not shown to be absent

    def test_every_directive_restates_a_table_row(self):
        table = self.contacts.set_index(["position_a", "position_b"])
        for row in self.directives.itertuples():
            source = table.loc[(row.position_a, row.position_b)]
            for g in ("g1", "g2"):
                self.assertEqual(int(source[f"contact_n_{g}"]), getattr(row, f"voters_{g}"))
                self.assertEqual(int(round(source[f"contact_{g}"] * source[f"contact_n_{g}"])), getattr(row, f"positives_{g}"))
                self.assertEqual(int(source[f"spare_negatives_{g}"]), getattr(row, f"spare_negatives_{g}"))        # the coverage column agrees
            self.assertEqual(source["conserved_in"], row.conserved_in)
            self.assertEqual(source["failing_groups"], row.failing_groups)
            self.assertEqual(row.tier == "conserved_in_all", bool(source["conserved_contact"]))                       # the headline is the first tier
        self.assertEqual(len(self.directives), int((self.contacts["conserved_in"] != "").sum()))                       # one per pair conserved somewhere

    def test_statements_are_generated_from_the_templates_only(self):
        pattern = re.compile(r"^[^:]+: (conserved in all groups|conserved in some groups only|conserved in one group and mostly absent from the others) "
                             r"\(conserved in [\w+ ]+; sequences with it / sequences voting: [^)]*?(; no voting sequences in [\w+ ]+)?\)\.$")
        for statement in self.directives["statement"]:
            self.assertRegex(statement, pattern)
            self.assertIsNone(FORBIDDEN.search(statement))

    def test_cofactor_directives_cover_each_class(self):
        rows, meta, features = [], [], []
        for s in range(12):
            sid = f"s{s}"
            group = "g1" if s < 6 else "g2"
            meta.append((sid, group, f"q{s}"))
            rows.append((sid, 1, 1, "HIS", "ND1", "N", 2.0, 0.0, 0.0))                       # a bonded His in every structure
            rows.append((sid, 2, 2, "SER", "OG", "O", 3.9 if group == "g2" else 8.0, 0.0, 0.0))        # a contact only in g2
            features.append((sid, "Cu", "CU", "CU", "metal", 0.0, 0.0, 0.0))
        atoms = pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])
        structures = pd.DataFrame(meta, columns=["structure_id", "group", "sequence_id"])
        features = pd.DataFrame(features, columns=["structure_id", "feature", "atom", "element", "kind", "x", "y", "z"])
        cof = ci.cofactor_contacts(atoms, features, coordination_model="by_pair", record_distance=9.0)
        cof["bond"] = cof["coord_tier"] == "bond"
        vote = ci.sequence_vote(cof, structures, ["feature", "position"], ["bond", "contact"])
        d = design_log.cofactor_directives(vote, structures, ["bond", "contact"]).set_index("directive_id")
        self.assertEqual(d.loc["retain_cofactor_contact:Cu:1:bond", "tier"], "conserved_in_all")
        self.assertEqual(d.loc["retain_cofactor_contact:Cu:1:contact", "tier"], "conserved_in_all")
        self.assertEqual(d.loc["retain_cofactor_contact:Cu:2:contact", "tier"], "family_specific")                    # in g2 only; g1 has voters and none positive
        self.assertNotIn("retain_cofactor_contact:Cu:2:bond", d.index)                                                # not conserved anywhere
        self.assertEqual(set(d["cofactor_class"]), {"bond", "contact"})
        self.assertEqual(design_log.cofactor_directives(vote.iloc[0:0], structures, ["bond"]).shape[0], 0)

    def test_the_log_carries_the_policy_and_the_directives_and_still_says_nothing_about_importance(self):
        atoms, structures = specific_dataset()
        contacts = ci.conserved_contacts(atoms, structures, coverage=True)
        directives = design_log.contact_directives(contacts, structures)
        with tempfile.TemporaryDirectory() as tmp:
            out = design_log.write_log(Path(tmp) / "log", {"conserved_contacts": contacts}, structures=structures, title="t", created="2026-11-02", directives=directives)
            again = design_log.write_log(Path(tmp) / "again", {"conserved_contacts": contacts}, structures=structures, title="t", created="2026-11-02", directives=directives)
            for name in ("manifest.json", "SUMMARY.md", "directives.csv"):
                self.assertEqual((out / name).read_bytes(), (again / name).read_bytes(), name)
            manifest = json.loads((out / "manifest.json").read_text())
            self.assertEqual(manifest["directive_policy"], design_log.DIRECTIVE_POLICY)
            self.assertEqual(manifest["directive_tiers"], {"family_specific": 12})
            self.assertEqual(manifest["tables"]["directives"]["rows"], 12)
            summary = (out / "SUMMARY.md").read_text()
            self.assertIn("## Directives", summary)
            self.assertIn("| family_specific | 12 |", summary)
            body = summary.replace(design_log.LOG_MEANING, "")
            self.assertIsNone(FORBIDDEN.search(body), FORBIDDEN.search(body))
            text = json.dumps({k: v for k, v in manifest.items() if k != "meaning"})
            self.assertIsNone(FORBIDDEN.search(text), FORBIDDEN.search(text))
            self.assertIsNone(FORBIDDEN.search((out / "directives.csv").read_text()))

    def test_directives_that_break_the_policy_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            for bad in (self.directives.assign(tier="must_keep"), self.directives.assign(claim_level="predictive"), self.directives.drop(columns=["statement"])):
                with self.assertRaises(ValueError):
                    design_log.write_log(Path(tmp) / "x", {"conserved_contacts": self.contacts}, structures=self.structures, title="t", directives=bad)

    def test_the_policy_text_itself_claims_nothing_about_importance(self):
        self.assertIsNone(FORBIDDEN.search(design_log.DIRECTIVE_POLICY))
        self.assertIn("makes no claim", design_log.DIRECTIVE_POLICY)


if __name__ == "__main__":
    unittest.main()
