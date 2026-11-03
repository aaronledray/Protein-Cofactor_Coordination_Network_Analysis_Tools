"""Single-protein mode: contacts across the structures of one sequence (evidence is counts and distances, never the vote)."""

import contextlib
import io
import json
import re
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from modules import conserve_cli, design_log, single_protein
from modules import conserved_interactions as ci

ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"
FORBIDDEN = re.compile(r"important|critical|essential|functional|stabiliz|drives|predict|must be retained|should be kept|conserved|conservation", re.IGNORECASE)
COLUMNS = ["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"]

# closest distance of each pair in each of four structures; the pairs share residue 1 and lie along different axes, so only they can touch
DISTANCES = {3: [3.6, 3.6, 3.6, 3.6],       # in contact in all four
             5: [3.6, 3.6, 4.4, 5.0],       # in two, median exactly at the cutoff
             7: [3.0, 9.0, 9.0, 9.0]}       # in one, the others far away
AXES = {3: (1, 0, 0), 5: (0, 1, 0), 7: (0, 0, 1)}


def fixture(missing=None):
    rows, meta = [], []
    for s in range(4):
        sid = f"s{s}"
        meta.append((sid, "all", "q", "ox" if s < 2 else "red"))
        rows.append((sid, 1, 1, "LEU", "CB", "C", 0.0, 0.0, 0.0))
        for position, distances in DISTANCES.items():
            if missing == (sid, position):
                continue
            rows.append((sid, position, position, "LEU", "CB", "C", *(distances[s] * a for a in AXES[position])))
    return pd.DataFrame(rows, columns=COLUMNS), pd.DataFrame(meta, columns=["structure_id", "group", "sequence_id", "state"])


class PersistenceTests(unittest.TestCase):
    def setUp(self):
        self.atoms, self.structures = fixture()
        self.table = single_protein.contact_persistence(self.atoms, self.structures).set_index("position_b")

    def test_counts_and_distances_per_pair(self):
        t = self.table
        self.assertEqual(list(t.index), [3, 5, 7])
        self.assertEqual(t["structures_present"].tolist(), [4, 4, 4])
        self.assertEqual(t["structures_in_contact"].tolist(), [4, 2, 1])
        self.assertEqual(t["fraction_in_contact"].tolist(), [1.0, 0.5, 0.25])
        self.assertEqual(t["closest_min"].tolist(), [3.6, 3.6, 3.0])
        self.assertEqual(t["closest_median"].tolist(), [3.6, 4.0, 6.0])
        self.assertEqual(t["closest_max"].tolist(), [3.6, 5.0, 6.0])                  # 9 A is recorded as the limit, cutoff + 2

    def test_near_cutoff_uses_the_median_and_the_margin(self):
        self.assertEqual(self.table["near_cutoff"].tolist(), [False, True, False])
        wide = single_protein.contact_persistence(self.atoms, self.structures, near_margin=0.5).set_index("position_b")
        self.assertEqual(wide["near_cutoff"].tolist(), [True, True, False])           # 3.6 is 0.4 from 4.0
        edge = single_protein.contact_persistence(self.atoms, self.structures, near_margin=0.4).set_index("position_b")
        self.assertTrue(edge.loc[3, "near_cutoff"])                                  # the margin is inclusive
        none = single_protein.contact_persistence(self.atoms, self.structures, near_margin=0.0).set_index("position_b")
        self.assertEqual(none["near_cutoff"].tolist(), [False, True, False])

    def test_a_structure_that_lacks_a_position_is_not_a_negative(self):
        atoms, structures = fixture(missing=("s3", 7))
        t = single_protein.contact_persistence(atoms, structures).set_index("position_b")
        self.assertEqual((int(t.loc[7, "structures_present"]), int(t.loc[7, "structures_in_contact"])), (3, 1))
        self.assertEqual(int(t.loc[3, "structures_present"]), 4)
        self.assertEqual(float(t.loc[7, "fraction_in_contact"]), round(1 / 3, 3))

    def test_chemistry_is_what_every_contact_structure_shows(self):
        rows = []
        for sid, distance in (("a", 3.0), ("b", 3.8)):
            rows += [(sid, 1, 1, "ASN", "ND2", "N", 0.0, 0.0, 0.0), (sid, 3, 3, "ASP", "OD1", "O", distance, 0.0, 0.0)]
        structures = pd.DataFrame({"structure_id": ["a", "b"], "group": "all", "sequence_id": "q"})
        table = single_protein.contact_persistence(pd.DataFrame(rows, columns=COLUMNS), structures)
        self.assertEqual(table["chemistry"].tolist(), ["sc_sc"])                      # a polar pair in a, a plain contact in b
        atoms = pd.DataFrame(rows, columns=COLUMNS)
        only_a = single_protein.contact_persistence(atoms[atoms["structure_id"] == "a"], structures[structures["structure_id"] == "a"])
        self.assertIn("polar", only_a["chemistry"].iloc[0])
        with self.assertRaisesRegex(ValueError, "same structures"):
            single_protein.contact_persistence(atoms, structures[structures["structure_id"] == "a"])

    def test_a_pair_never_in_contact_is_not_listed(self):
        atoms, structures = fixture()
        atoms.loc[(atoms["position"] == 7), ["x", "y", "z"]] = [0.0, 0.0, 9.0]
        self.assertNotIn(7, single_protein.contact_persistence(atoms, structures)["position_b"].tolist())

    def test_states_are_counted_separately(self):
        t = self.table
        self.assertEqual((t["present_ox"].tolist(), t["in_contact_ox"].tolist()), ([2, 2, 2], [2, 2, 1]))
        self.assertEqual((t["present_red"].tolist(), t["in_contact_red"].tolist()), ([2, 2, 2], [2, 0, 0]))
        one_state = single_protein.contact_persistence(self.atoms, self.structures.assign(state="ox"))
        self.assertFalse([c for c in one_state.columns if c.endswith("_ox")])           # one state is no stratification

    def test_one_sequence_is_required_and_state_labels_are_checked(self):
        with self.assertRaisesRegex(ValueError, "exactly one sequence_id"):
            single_protein.contact_persistence(self.atoms, self.structures.assign(sequence_id=["a", "b", "a", "b"]))
        with self.assertRaisesRegex(ValueError, "state labels"):
            single_protein.contact_persistence(self.atoms, self.structures.assign(state="ox state"))

    def test_the_legend_travels_and_names_the_cutoff(self):
        rules = ci.rules_of(self.table.reset_index())
        text = " ".join(rules["rule"])
        self.assertIn("within 0.3 A of the 4 A cutoff", text)
        self.assertIn("recorded as 6", text)
        self.assertEqual(set(rules["table"]), {"residue_contacts", "contact_persistence"})

    def test_real_structure_equals_the_contact_rule(self):
        copies = [f"c{i}" for i in range(3)]
        atoms = pd.concat([ci.load_atoms(str(PLASTOCYANIN), "A", structure_id=s) for s in copies], ignore_index=True)
        structures = pd.DataFrame({"structure_id": copies, "group": "all", "sequence_id": "q"})
        table = single_protein.contact_persistence(atoms, structures)
        direct = ci.residue_contacts(atoms[atoms["structure_id"] == "c0"])
        self.assertEqual(len(table), len(direct))
        self.assertTrue((table["structures_in_contact"] == 3).all() and (table["structures_present"] == 3).all())
        merged = table.merge(direct, on=["position_a", "position_b"])
        np.testing.assert_allclose(merged["closest_median"], merged["min_distance"], atol=1e-3)          # one distance in every copy
        self.assertEqual(sorted(merged["chemistry"].str.contains("polar")), sorted(direct["polar"]))


class CofactorPersistenceTests(unittest.TestCase):
    def test_counts_per_class_and_distances(self):
        rows, meta, features = [], [], []
        for s in range(4):
            sid = f"s{s}"
            meta.append((sid, "all", "q"))
            rows.append((sid, 1, 1, "HIS", "ND1", "N", 2.0, 0.0, 0.0))
            rows.append((sid, 2, 2, "SER", "OG", "O", 3.0 if s < 2 else 7.0, 0.0, 0.0))
            features.append((sid, "Cu", "CU", "CU", "metal", 0.0, 0.0, 0.0))
        atoms = pd.DataFrame(rows, columns=COLUMNS)
        structures = pd.DataFrame(meta, columns=["structure_id", "group", "sequence_id"])
        features = pd.DataFrame(features, columns=["structure_id", "feature", "atom", "element", "kind", "x", "y", "z"])
        cofactor = ci.cofactor_contacts(atoms, features, coordination_model="by_pair", record_distance=9.0)
        cofactor["bond"] = cofactor["coord_tier"] == "bond"
        table = single_protein.cofactor_persistence(cofactor, structures, ["bond", "contact"]).set_index("position")
        self.assertEqual((table.loc[1, "bond_structures"], table.loc[1, "structures_recorded"]), (4, 4))
        self.assertEqual((table.loc[2, "contact_structures"], table.loc[2, "bond_structures"]), (2, 0))
        self.assertEqual((table.loc[2, "distance_min"], table.loc[2, "distance_max"]), (3.0, 7.0))
        self.assertIn("structures_recorded", " ".join(ci.rules_of(table.reset_index())["name"]))
        self.assertNotIn(3, table.index)


class DirectiveTests(unittest.TestCase):
    def test_tier_boundaries(self):
        tier = design_log.observed_tier
        cases = {(1, 1): "observed_in_one_structure", (1, 8): "observed_in_one_structure", (2, 2): "observed_in_all_structures",
                 (8, 8): "observed_in_all_structures", (6, 8): "observed_in_some_structures", (2, 3): "observed_in_some_structures"}
        for (seen, present), expected in cases.items():
            self.assertEqual(tier(seen, present), expected, (seen, present))

    def test_directives_restate_the_table_and_use_the_template(self):
        atoms, structures = fixture()
        table = single_protein.contact_persistence(atoms, structures)
        directives = design_log.observed_contact_directives(table, structures).set_index("position_b")
        self.assertEqual(directives["tier"].tolist(), ["observed_in_all_structures", "observed_in_some_structures", "observed_in_one_structure"])
        self.assertEqual(directives["claim_level"].unique().tolist(), ["observed"])
        self.assertEqual(directives["structures_in_contact"].tolist(), [4, 2, 1])
        self.assertEqual(directives.loc[5, "near_cutoff"], True)
        pattern = re.compile(r"^contact 1-\d: observed in (all|some) structures|one structure \(in contact in \d of \d structures with both positions present; "
                             r"closest distance median \d\.\d\d A, range \d\.\d\d to \d\.\d\d A(; within the near-cutoff margin of the cutoff)?\)\.$")
        for statement in directives["statement"]:
            self.assertRegex(statement, pattern)
            self.assertIsNone(FORBIDDEN.search(statement), statement)
        self.assertIn("within the near-cutoff margin", directives.loc[5, "statement"])
        self.assertNotIn("near-cutoff", directives.loc[3, "statement"])

    def test_one_structure_says_nothing_about_persistence(self):
        atoms, structures = fixture()
        one, one_structures = atoms[atoms["structure_id"] == "s0"], structures[structures["structure_id"] == "s0"]
        directives = design_log.observed_contact_directives(single_protein.contact_persistence(one, one_structures), one_structures)
        self.assertEqual(set(directives["tier"]), {"observed_in_one_structure"})
        for statement in directives["statement"]:
            self.assertIn("nothing is shown about persistence", statement)

    def test_the_wording_of_the_mode_claims_nothing_about_importance(self):
        self.assertIsNone(FORBIDDEN.search(design_log.SINGLE_DIRECTIVE_POLICY.replace("evidence of conservation", "")), design_log.SINGLE_DIRECTIVE_POLICY)
        self.assertIn("not evidence of conservation", design_log.SINGLE_MEANING)

    def test_writer_modes_accept_only_their_own_directives(self):
        atoms, structures = fixture()
        table = single_protein.contact_persistence(atoms, structures)
        observed = design_log.observed_contact_directives(table, structures)
        with tempfile.TemporaryDirectory() as tmp:
            out = design_log.write_log(Path(tmp) / "ok", {"contact_persistence": table}, structures=structures, title="t", directives=observed, mode="single_protein", created="2026-11-03")
            manifest = json.loads((out / "manifest.json").read_text())
            summary = (out / "SUMMARY.md").read_text()
            self.assertEqual((manifest["mode"], manifest["directive_policy"]), ("single_protein", design_log.SINGLE_DIRECTIVE_POLICY))
            self.assertEqual(manifest["dataset"]["states"], {"ox": 2, "red": 2})
            self.assertEqual(manifest["meaning"], design_log.SINGLE_MEANING)
            self.assertNotIn("groups", manifest["dataset"])
            self.assertIn("| observed_in_one_structure | 1 |", summary)
            body = summary.replace(design_log.SINGLE_MEANING, "").replace(design_log.SINGLE_DIRECTIVE_POLICY, "")
            self.assertIsNone(FORBIDDEN.search(body), FORBIDDEN.search(body))
            self.assertIn("state", pd.read_csv(out / "structures.csv").columns)
            with self.assertRaises(ValueError):                                          # observed directives in the vote mode
                design_log.write_log(Path(tmp) / "a", {"contact_persistence": table}, structures=structures, title="t", directives=observed)
            vote = observed.assign(tier="conserved_in_all", claim_level="descriptive")
            with self.assertRaises(ValueError):                                          # vote directives in the single-protein mode
                design_log.write_log(Path(tmp) / "b", {"contact_persistence": table}, structures=structures, title="t", directives=vote, mode="single_protein")
            with self.assertRaises(ValueError):                                          # several sequences
                design_log.write_log(Path(tmp) / "c", {}, structures=structures.assign(sequence_id=["a", "b", "a", "b"]), title="t", mode="single_protein")
            with self.assertRaises(ValueError):
                design_log.write_log(Path(tmp) / "d", {}, structures=structures, title="t", mode="other")


def run(argv):
    err = io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        code = conserve_cli.main(argv)
    return code, err.getvalue()


class CommandTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def manifest(self, name="m.csv", n=4, **columns):
        frame = pd.DataFrame({"structure_id": [f"s{i}" for i in range(n)], "path": [str(PLASTOCYANIN)] * n, "sequence_id": "q", **columns})
        frame.to_csv(self.dir / name, index=False)
        return str(self.dir / name)

    def test_log_equals_the_api_and_has_no_vote(self):
        code, err = run(["--manifest", self.manifest(state=["ox", "ox", "red", "red"]), "--single-protein", "--metals", "CU", "--coordination-tiers",
                         "--directives", "--output-dir", str(self.dir / "log")])
        self.assertEqual(code, 0, err)
        log = self.dir / "log"
        self.assertEqual({p.stem for p in log.glob("*.csv")}, {"structures", "contact_persistence", "cofactor_persistence", "audit_backbone_only", "directives"})
        table = pd.read_csv(log / "contact_persistence.csv")
        atoms = pd.concat([ci.load_atoms(str(PLASTOCYANIN), "A", structure_id=f"s{i}") for i in range(4)], ignore_index=True)
        structures = pd.DataFrame({"structure_id": [f"s{i}" for i in range(4)], "group": "all", "sequence_id": "q", "state": ["ox", "ox", "red", "red"]})
        direct = single_protein.contact_persistence(atoms, structures)
        self.assertEqual(len(table), len(direct))
        self.assertEqual(table["structures_in_contact"].tolist(), direct["structures_in_contact"].tolist())
        cofactor = pd.read_csv(log / "cofactor_persistence.csv").set_index("position")
        self.assertEqual({int(p) for p in cofactor.index[cofactor["bond_structures"] == 4]}, {37, 84, 87})
        manifest = json.loads((log / "manifest.json").read_text())
        self.assertEqual((manifest["mode"], manifest["title"]), ("single_protein", "Single-protein contact log"))
        self.assertNotIn("support", manifest["parameters"])
        directives = pd.read_csv(log / "directives.csv")
        self.assertEqual(set(directives["claim_level"]), {"observed"})
        self.assertEqual(set(directives["directive"]), {"observed_contact", "observed_cofactor_contact"})
        raw = (log / "directives.csv").read_text().splitlines()
        contact_lines = [line.split(",") for line in raw[1:] if line.startswith("observed_contact:")]
        self.assertGreater(len(contact_lines), 100)
        self.assertTrue(all(re.fullmatch(r"\d+", fields[3]) and re.fullmatch(r"\d+", fields[4]) for fields in contact_lines))      # positions are not written as 1.0

    def test_one_structure_is_allowed_and_noted(self):
        code, err = run(["--manifest", self.manifest(n=1), "--single-protein", "--output-dir", str(self.dir / "one")])
        self.assertEqual(code, 0, err)
        manifest = json.loads((self.dir / "one/manifest.json").read_text())
        self.assertTrue(any("nothing can be said about persistence" in n for n in manifest["notes"]))
        self.assertIn("nothing can be said", (self.dir / "one/SUMMARY.md").read_text())

    def test_refusals(self):
        two = self.manifest("two.csv", sequence_id=["a", "a", "b", "b"])
        code, err = run(["--manifest", two, "--single-protein", "--output-dir", str(self.dir / "o")])
        self.assertEqual((code, "exactly one sequence_id" in err), (2, True))
        base = self.manifest()
        for flags in (["--support", "0.8"], ["--min-sequences", "2"], ["--sulfur-acceptor"], ["--neighbor-contacts"], ["--metals", "CU", "--anchor", "84"]):
            code, err = run(["--manifest", base, "--single-protein", "--output-dir", str(self.dir / "o"), *flags])
            self.assertEqual((code, "do not apply with --single-protein" in err), (2, True), flags)
        code, err = run(["--manifest", self.manifest("bad.csv", state=["a b"] * 4), "--single-protein", "--output-dir", str(self.dir / "o")])
        self.assertEqual((code, "state labels" in err), (2, True))
        self.assertFalse((self.dir / "o").exists())

    def test_default_mode_is_unchanged_and_ignores_a_state_column(self):
        code, err = run(["--manifest", self.manifest(sequence_id=["a", "b", "c", "d"], state="ox"), "--output-dir", str(self.dir / "vote"), "--min-sequences", "1"])
        self.assertEqual(code, 0, err)
        manifest = json.loads((self.dir / "vote/manifest.json").read_text())
        self.assertEqual(manifest["mode"], "vote")
        self.assertIn("groups", manifest["dataset"])
        self.assertTrue((self.dir / "vote/conserved_contacts.csv").exists())


if __name__ == "__main__":
    unittest.main()
