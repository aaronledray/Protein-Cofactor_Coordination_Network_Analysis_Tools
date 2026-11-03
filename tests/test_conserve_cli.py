"""``sscna conserve``: a thin layer, so its log must equal what the API functions return for the same inputs."""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

import sscna_cli
from modules import conserve_cli, design_log
from modules import conserved_interactions as ci

ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"
SEQUENCES = [f"seq{i}" for i in range(5)]


def run(argv):
    err = io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        code = conserve_cli.main(argv)
    return code, err.getvalue()


class ConserveTests(unittest.TestCase):
    """Five copies of one structure as five sequences: every contact that exists is positive in every sequence."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.dir = Path(cls.tmp.name)
        cls.manifest = cls.dir / "manifest.csv"
        pd.DataFrame({"structure_id": [f"s{i}" for i in range(5)], "path": [str(PLASTOCYANIN)] * 5, "sequence_id": SEQUENCES}).to_csv(cls.manifest, index=False)
        cls.out = cls.dir / "log"
        code, err = run(["--manifest", str(cls.manifest), "--output-dir", str(cls.out), "--metals", "CU", "--coordination-tiers",
                         "--metal-ligands", "--directives", "--anchor", "84", "--title", "t"])
        assert code == 0, err
        cls.log = {p.stem: pd.read_csv(p) for p in cls.out.glob("*.csv")}
        cls.manifest_json = json.loads((cls.out / "manifest.json").read_text())

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_files_and_manifest(self):
        self.assertEqual(set(self.log), {"structures", "conserved_contacts", "audit_backbone_only", "cofactor_contacts", "audit_anchor", "metal_ligands", "directives"})
        self.assertTrue((self.out / "SUMMARY.md").exists())
        parameters = self.manifest_json["parameters"]
        self.assertEqual((parameters["support"], parameters["min_sequences"], parameters["cutoff"]), (0.9, 5, 4.0))
        self.assertEqual(self.manifest_json["dataset"]["groups"]["all"]["sequences"], 5)

    def test_contacts_equal_the_api(self):
        atoms = pd.concat([ci.load_atoms(str(PLASTOCYANIN), "A", structure_id=f"s{i}") for i in range(5)], ignore_index=True)
        structures = pd.DataFrame({"structure_id": [f"s{i}" for i in range(5)], "group": "all", "sequence_id": SEQUENCES})
        direct = ci.conserved_contacts(atoms, structures, coverage=True)
        mine = self.log["conserved_contacts"]
        keys = ["position_a", "position_b"]
        self.assertEqual(len(mine), len(direct))
        self.assertEqual(int(mine["conserved_contact"].sum()), int(direct["conserved_contact"].sum()))
        pd.testing.assert_frame_equal(mine.sort_values(keys, ignore_index=True)[keys + ["conserved_contact"]],
                                      direct.sort_values(keys, ignore_index=True)[keys + ["conserved_contact"]], check_dtype=False)
        self.assertGreater(int(mine["conserved_contact"].sum()), 100)

    def test_cofactor_tiers_and_directives(self):
        cofactor = self.log["cofactor_contacts"].set_index("position")
        self.assertEqual({int(p) for p in cofactor.index[cofactor["conserved_bond"]]}, {37, 84, 87})     # His37, Cys84, His87 (reference numbering of 1ag6)
        directives = self.log["directives"]
        self.assertEqual(set(directives["tier"]), {"conserved_in_all"})
        self.assertEqual(set(directives["directive"]), {"conserved_contact", "conserved_cofactor_contact"})
        self.assertEqual(int((directives["directive"] == "conserved_contact").sum()), int(self.log["conserved_contacts"]["conserved_contact"].sum()))

    def test_audits_and_ligands(self):
        self.assertTrue(self.log["audit_anchor"]["ok"].all())
        self.assertFalse(self.log["audit_backbone_only"]["backbone_only"].any())
        self.assertEqual(list(self.log["metal_ligands"].columns[:2]), ["structure_id", "metal_chain"])      # the copper of 1ag6 has no hetero ligand: an empty table

    def test_metal_ligands_equal_the_api(self):
        myoglobin = ROOT / "reference_structures/3_Myoglobin/1a6m.pdb"
        pd.DataFrame({"structure_id": ["m1", "m2"], "path": [str(myoglobin)] * 2, "sequence_id": ["x", "y"]}).to_csv(self.dir / "myo.csv", index=False)
        code, err = run(["--manifest", str(self.dir / "myo.csv"), "--output-dir", str(self.dir / "myo_log"), "--metal-ligands", "--min-sequences", "1"])
        self.assertEqual(code, 0, err)
        mine = pd.read_csv(self.dir / "myo_log/metal_ligands.csv")
        direct = pd.concat([ci.load_metal_ligands(str(myoglobin), structure_id=s) for s in ("m1", "m2")], ignore_index=True)
        self.assertGreater(len(direct), 0)
        self.assertEqual(len(mine), len(direct))
        self.assertEqual(sorted(mine["ligand_resname"]), sorted(direct["ligand_resname"]))

    def test_defaults_are_off(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, err = run(["--manifest", str(self.manifest), "--output-dir", tmp + "/plain"])
            self.assertEqual(code, 0, err)
            names = {p.stem for p in Path(tmp, "plain").glob("*.csv")}
            self.assertEqual(names, {"structures", "conserved_contacts", "audit_backbone_only"})
            contacts = pd.read_csv(Path(tmp, "plain/conserved_contacts.csv"))
            self.assertFalse([c for c in contacts.columns if any(c.startswith(f"{flag}_") for flag in ci.SULFUR_FLAGS)])      # no sulfur-acceptor flags
            manifest = json.loads(Path(tmp, "plain/manifest.json").read_text())
            self.assertNotIn("directive_policy", manifest)

    def test_sulfur_and_neighbor_options_equal_the_api(self):
        anhydrase = ROOT / "reference_structures/3_Metals/1ca2.pdb"                       # plastocyanin has no neighbor polar contact to add
        ids = [f"c{i}" for i in range(5)]
        pd.DataFrame({"structure_id": ids, "path": [str(anhydrase)] * 5, "sequence_id": SEQUENCES}).to_csv(self.dir / "ca.csv", index=False)
        atoms = pd.concat([ci.load_atoms(str(anhydrase), "A", structure_id=s) for s in ids], ignore_index=True)
        structures = pd.DataFrame({"structure_id": ids, "group": "all", "sequence_id": SEQUENCES})
        metals = pd.concat([ci.load_metals(str(anhydrase), structure_id=s) for s in ids], ignore_index=True)
        with tempfile.TemporaryDirectory() as tmp:
            for flags, kwargs in ((["--sulfur-acceptor", "--neighbor-contacts"], {"sulfur_acceptor": True, "neighbor_contacts": True, "metals": metals}),
                                  ([], {})):
                code, err = run(["--manifest", str(self.dir / "ca.csv"), "--output-dir", tmp + "/log", "--overwrite", *flags])
                self.assertEqual(code, 0, err)
                mine = pd.read_csv(Path(tmp, "log/conserved_contacts.csv"))
                direct = ci.conserved_contacts(atoms, structures, coverage=True, **kwargs)
                self.assertEqual((len(mine), int(mine["conserved_contact"].sum())), (len(direct), int(direct["conserved_contact"].sum())), flags)
                if flags:
                    with_options = len(mine)
                    self.assertIn("polar_s_all", mine.columns)
        self.assertGreater(with_options, len(mine))                                       # the options reach the analysis: neighbor pairs add candidates

    def test_position_map_relabels_and_drops(self):
        residues = ci.load_atoms(str(PLASTOCYANIN), "A").drop_duplicates("position")["position"].tolist()
        kept = [n for n in residues if n <= 50]
        rows = [{"structure_id": f"s{i}", "residue_number": n, "insertion_code": "", "position": f"P{n}"} for i in range(5) for n in kept]
        with tempfile.TemporaryDirectory() as tmp:
            mapping = Path(tmp) / "map.csv"
            pd.DataFrame(rows).to_csv(mapping, index=False)
            code, err = run(["--manifest", str(self.manifest), "--position-map", str(mapping), "--output-dir", tmp + "/log"])
            self.assertEqual(code, 0, err)
            contacts = pd.read_csv(Path(tmp, "log/conserved_contacts.csv"))
        self.assertTrue(contacts["position_a"].astype(str).str.startswith("P").all())
        self.assertTrue(contacts["position_b"].astype(str).str.startswith("P").all())
        self.assertGreater(len(contacts), 0)


class InputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def manifest(self, **columns):
        frame = pd.DataFrame({"structure_id": ["a", "b"], "path": [str(PLASTOCYANIN)] * 2, "sequence_id": ["x", "y"], **columns})
        path = self.dir / "m.csv"
        frame.to_csv(path, index=False)
        return path

    def test_manifest_problems_are_reported_not_raised(self):
        pd.DataFrame({"structure_id": ["a"], "path": [str(PLASTOCYANIN)]}).to_csv(self.dir / "m0.csv", index=False)
        code, err = run(["--manifest", str(self.dir / "m0.csv"), "--output-dir", str(self.dir / "o")])
        self.assertEqual(code, 2)
        self.assertIn("lacks column", err)
        pd.DataFrame({"structure_id": ["a", "a"], "path": [str(PLASTOCYANIN)] * 2, "sequence_id": ["x", "y"]}).to_csv(self.dir / "m1.csv", index=False)
        code, err = run(["--manifest", str(self.dir / "m1.csv"), "--output-dir", str(self.dir / "o")])
        self.assertEqual((code, "duplicate" in err), (2, True))
        pd.DataFrame({"structure_id": ["a"], "path": ["nowhere.cif"], "sequence_id": ["x"]}).to_csv(self.dir / "m2.csv", index=False)
        code, err = run(["--manifest", str(self.dir / "m2.csv"), "--output-dir", str(self.dir / "o")])
        self.assertEqual((code, "not found" in err), (2, True))
        self.assertFalse((self.dir / "o").exists())

    def test_relative_paths_are_taken_from_the_manifest_folder(self):
        (self.dir / "sub").mkdir()
        (self.dir / "sub/copy.cif").write_bytes(PLASTOCYANIN.read_bytes())
        pd.DataFrame({"structure_id": ["a"], "path": ["sub/copy.cif"], "sequence_id": ["x"]}).to_csv(self.dir / "m.csv", index=False)
        table = conserve_cli.read_manifest(self.dir / "m.csv")
        self.assertEqual(table["path"].iloc[0], str(self.dir / "sub/copy.cif"))
        self.assertEqual((table["group"].iloc[0], table["chain"].iloc[0]), ("all", "A"))

    def test_option_dependencies_and_overwrite(self):
        manifest = self.manifest()
        for flags, text in ((["--coordination-tiers"], "--metals"), (["--anchor", "84"], "--metals")):
            code, err = run(["--manifest", str(manifest), "--output-dir", str(self.dir / "o"), *flags])
            self.assertEqual(code, 2)
            self.assertIn(text, err)
        args = ["--manifest", str(manifest), "--output-dir", str(self.dir / "log"), "--min-sequences", "1"]
        self.assertEqual(run(args)[0], 0)
        self.assertGreater(int(pd.read_csv(self.dir / "log/conserved_contacts.csv")["conserved_contact"].sum()), 100)      # --min-sequences reaches the vote
        code, err = run(args)
        self.assertEqual((code, "overwrite" in err), (2, True))
        self.assertEqual(run(args + ["--overwrite"])[0], 0)

    def test_position_map_problems(self):
        manifest = self.manifest()
        pd.DataFrame({"structure_id": ["a", "a"], "residue_number": [1, 2], "position": [5, 5]}).to_csv(self.dir / "clash.csv", index=False)
        code, err = run(["--manifest", str(manifest), "--position-map", str(self.dir / "clash.csv"), "--output-dir", str(self.dir / "o")])
        self.assertEqual((code, "two residues" in err), (2, True))
        pd.DataFrame({"structure_id": ["a"], "residue_number": [1], "position": [1]}).to_csv(self.dir / "partial.csv", index=False)
        code, err = run(["--manifest", str(manifest), "--position-map", str(self.dir / "partial.csv"), "--output-dir", str(self.dir / "o")])
        self.assertEqual((code, "no rows for structure b" in err), (2, True))

    def test_small_group_is_noted_and_nothing_is_conserved(self):
        manifest = self.manifest()
        code, err = run(["--manifest", str(manifest), "--output-dir", str(self.dir / "log"), "--directives"])
        self.assertEqual(code, 0, err)
        notes = json.loads((self.dir / "log/manifest.json").read_text())["notes"]
        self.assertTrue(any("fewer than 5 sequences" in n for n in notes))
        self.assertTrue(any("nothing is conserved" in n for n in notes))
        self.assertFalse((self.dir / "log/directives.csv").exists())

    def test_sscna_dispatch(self):
        with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as stop:
                sscna_cli.main(["conserve", "--help"])
            self.assertEqual(stop.exception.code, 0)
            self.assertEqual(sscna_cli.main(["conserve", "--manifest", str(self.dir / "none.csv"), "--output-dir", str(self.dir / "o")]), 2)


if __name__ == "__main__":
    unittest.main()
