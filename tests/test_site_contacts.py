"""Tests for modules/site_contacts.py (atom-level site contacts with every alternate conformation)."""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from modules import site_contacts as sc


def atom_line(serial, name, resname, chain, resnum, xyz, *, altloc=" ", occupancy=1.0, element=None, hetero=False, icode=" "):
    element = element or name.strip()[0]
    field = f" {name:<3}" if len(element) == 1 and len(name) < 4 else f"{name:<4}"
    x, y, z = xyz
    record = "HETATM" if hetero else "ATOM  "
    return (f"{record}{serial:5d} {field}{altloc:1s}{resname:>3s} {chain:1s}{resnum:4d}{icode:1s}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{20.0:6.2f}          {element:>2s}")


class Builder:
    """Collects synthetic PDB records; every residue is one chain A residue unless stated."""

    def __init__(self):
        self.lines = []

    def add(self, *args, **kwargs):
        self.lines.append(atom_line(len(self.lines) + 1, *args, **kwargs))
        return self

    def write(self, folder, name="synthetic.pdb"):
        path = Path(folder) / name
        path.write_text("\n".join(self.lines + ["END", ""]))
        return path


def heme(builder, resnum=155):
    """A minimal HEM: the A and D propionate oxygens, an iron and one ring nitrogen."""
    builder.add("O1A", "HEM", "A", resnum, (0.0, 0.0, 0.0), hetero=True)
    builder.add("O2A", "HEM", "A", resnum, (1.2, 0.0, 0.0), hetero=True)
    builder.add("O1D", "HEM", "A", resnum, (0.0, 20.0, 0.0), hetero=True)
    builder.add("O2D", "HEM", "A", resnum, (1.2, 20.0, 0.0), hetero=True)
    builder.add("FE", "HEM", "A", resnum, (10.0, 10.0, 10.0), hetero=True, element="FE")
    return builder


PROPIONATE_A = {"label": "propionate_A", "resname": "HEM", "atoms": ["O1A", "O2A"]}
PROPIONATE_D = {"label": "propionate_D", "resname": "HEM", "atoms": ["O1D", "O2D"]}


class SiteContactsTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.folder = Path(self._tmp.name)

    def run_site(self, builder, sources, **kwargs):
        return sc.site_contacts(builder.write(self.folder), sources, **kwargs)

    # -- alternate conformations
    def two_conformer_his(self, occ_a=0.6, occ_b=0.4):
        b = heme(Builder())
        # NE2 of residue 64 has two conformers 2.9 A and 3.4 A from O1A; its CE1 is not disordered.
        b.add("NE2", "HIS", "A", 64, (0.0, 2.9, 0.0), altloc="A", occupancy=occ_a)
        b.add("NE2", "HIS", "A", 64, (0.0, 3.4, 0.0), altloc="B", occupancy=occ_b)
        b.add("CE1", "HIS", "A", 64, (0.0, 3.0, 1.0))
        return b

    def test_every_conformer_is_a_row_with_altloc_and_occupancy(self):
        result = self.run_site(self.two_conformer_his(), [PROPIONATE_A])
        rows = result["site_atom_contacts"]
        ne2 = rows[(rows.source_atom == "O1A") & (rows.partner_atom == "NE2")].sort_values("partner_altloc")
        self.assertEqual(list(ne2.partner_altloc), ["A", "B"])
        self.assertEqual(list(ne2.partner_occupancy), [0.6, 0.4])
        self.assertEqual(list(ne2.distance_A), [2.9, 3.4])
        self.assertTrue(ne2.polar_candidate.all())
        self.assertEqual(list(ne2.close_polar), [True, False])          # 2.9 <= 3.2 < 3.4
        self.assertEqual(result["manifest"]["counts"]["n_atoms_with_altlocs"], 1)

    def test_highest_occupancy_policy_picks_the_larger_occupancy_not_the_first_letter(self):
        result = self.run_site(self.two_conformer_his(occ_a=0.3, occ_b=0.7), [PROPIONATE_A], altloc_policy="highest_occupancy")
        rows = result["site_atom_contacts"]
        ne2 = rows[(rows.source_atom == "O1A") & (rows.partner_atom == "NE2")]
        self.assertEqual(list(ne2.partner_altloc), ["B"])
        self.assertEqual(set(result["site_atom_contacts"].altloc_policy), {"highest_occupancy"})

    def test_equal_occupancy_picks_the_first_listed_like_biopython(self):
        result = self.run_site(self.two_conformer_his(occ_a=0.5, occ_b=0.5), [PROPIONATE_A], altloc_policy="highest_occupancy")
        rows = result["site_atom_contacts"]
        self.assertEqual(list(rows[rows.partner_atom == "NE2"].partner_altloc.unique()), ["A"])

    def test_residue_summary_comparison_columns(self):
        b = heme(Builder())
        b.add("NE2", "HIS", "A", 64, (0.0, 2.9, 0.0), altloc="A", occupancy=0.3)
        b.add("NE2", "HIS", "A", 64, (0.0, 3.4, 0.0), altloc="B", occupancy=0.7)
        b.add("OG", "SER", "A", 92, (0.0, -4.4, 0.0))
        result = self.run_site(b, [PROPIONATE_A])
        summary = result["site_residue_summary"]
        his = summary[(summary.partner_resname == "HIS") & (summary.source_group == "propionate_A")].iloc[0]
        self.assertEqual(his.min_distance_A, 2.9)                       # the minor conformer is visible
        self.assertEqual(his.partner_altlocs_present, "A|B")
        self.assertTrue(his.in_highest_occupancy)
        self.assertEqual(his.min_distance_highest_occupancy_A, 3.4)     # what the highest-occupancy view sees
        ser = summary[summary.partner_resname == "SER"].iloc[0]
        self.assertEqual(ser.min_distance_A, ser.min_distance_highest_occupancy_A)

    def test_minor_only_contact_marks_not_in_highest_occupancy(self):
        # the major conformer (A, 0.7) lies beyond the neighbor cutoff; only the minor one (B, 0.3) is a neighbor
        b = heme(Builder())
        b.add("NE2", "HIS", "A", 64, (0.0, 6.0, 0.0), altloc="A", occupancy=0.7)
        b.add("NE2", "HIS", "A", 64, (0.0, 3.0, 0.0), altloc="B", occupancy=0.3)
        result = self.run_site(b, [PROPIONATE_A])
        summary = result["site_residue_summary"]
        his = summary[summary.partner_resname == "HIS"].iloc[0]
        self.assertFalse(his.in_highest_occupancy)
        self.assertTrue(pd.isna(his.min_distance_highest_occupancy_A))

    def test_altloc_compatibility_flag(self):
        b = Builder()
        b.add("O1A", "HEM", "A", 155, (0.0, 0.0, 0.0), altloc="A", occupancy=0.5, hetero=True)
        b.add("O1A", "HEM", "A", 155, (0.0, 0.2, 0.0), altloc="B", occupancy=0.5, hetero=True)
        b.add("NE2", "HIS", "A", 64, (0.0, 3.0, 0.0), altloc="A", occupancy=0.5)
        b.add("NE2", "HIS", "A", 64, (0.0, 3.1, 0.0), altloc="B", occupancy=0.5)
        b.add("OG", "SER", "A", 92, (0.0, -3.0, 0.0))
        rows = self.run_site(b, [PROPIONATE_A])["site_atom_contacts"]
        his = rows[rows.partner_atom == "NE2"]
        self.assertEqual(len(his), 4)                                    # 2 source x 2 partner conformers
        self.assertEqual(int(his.altloc_compatible.sum()), 2)            # A-A and B-B only
        self.assertTrue(rows[rows.partner_atom == "OG"].altloc_compatible.all())

    # -- chemistry flags, identity, selection
    def test_polar_close_and_ionic_flags(self):
        b = heme(Builder())
        b.add("NH1", "ARG", "A", 45, (0.0, 3.8, 0.0))                  # ionic (<= 4.0) but polar needs <= 3.5
        b.add("NZ", "LYS", "A", 10, (0.0, -3.0, 0.0))                  # polar, close polar and ionic
        b.add("CB", "ALA", "A", 20, (0.0, 0.0, 4.0))                   # neighbor only (carbon)
        b.add("NE", "ARG", "A", 46, (0.0, 0.0, -4.4))                  # ionic window exceeded (4.4 > 4.0), neighbor only
        rows = self.run_site(b, [PROPIONATE_A])["site_atom_contacts"]
        o1a = rows[rows.source_atom == "O1A"].set_index(["partner_resname", "partner_resnum", "partner_atom"])
        arg45 = o1a.loc[("ARG", 45, "NH1")]
        self.assertFalse(arg45.polar_candidate)
        self.assertTrue(arg45.ionic_candidate)
        lys = o1a.loc[("LYS", 10, "NZ")]
        self.assertTrue(lys.polar_candidate and lys.close_polar and lys.ionic_candidate)
        ala = o1a.loc[("ALA", 20, "CB")]
        self.assertFalse(ala.polar_candidate or ala.ionic_candidate)
        self.assertFalse(o1a.loc[("ARG", 46, "NE")].ionic_candidate)

    def test_same_residue_pairs_are_excluded_unless_requested(self):
        b = heme(Builder())
        default = self.run_site(b, [PROPIONATE_A])["site_atom_contacts"]
        self.assertFalse((default.partner_resname == "HEM").any())
        kept = self.run_site(b, [PROPIONATE_A], include_same_residue=True)["site_atom_contacts"]
        own = kept[kept.partner_resname == "HEM"]
        self.assertTrue(own.same_residue.all())
        self.assertFalse(((own.source_atom == own.partner_atom)).any())  # an atom is never its own neighbor

    def test_hydrogens_are_ignored_and_counted(self):
        b = heme(Builder())
        b.add("NE2", "HIS", "A", 64, (0.0, 3.0, 0.0))
        b.add("HE2", "HIS", "A", 64, (0.0, 2.0, 0.0), element="H")
        result = self.run_site(b, [PROPIONATE_A])
        self.assertNotIn("HE2", set(result["site_atom_contacts"].partner_atom))
        self.assertEqual(result["manifest"]["counts"]["n_hydrogens_ignored"], 1)

    def test_identity_is_reported_verbatim_with_insertion_code_and_hetero_flag(self):
        b = heme(Builder())
        b.add("NE2", "HIS", "A", 64, (0.0, 3.0, 0.0), icode="B")
        b.add("O", "HOH", "A", 301, (3.0, 0.0, 0.0), hetero=True)
        rows = self.run_site(b, [PROPIONATE_A])["site_atom_contacts"]
        his = rows[rows.partner_atom == "NE2"].iloc[0]
        self.assertEqual((his.partner_chain, his.partner_resnum, his.partner_icode, his.partner_hetero_flag), ("A", 64, "B", ""))
        water = rows[rows.partner_resname == "HOH"].iloc[0]
        self.assertEqual((water.partner_kind, water.partner_hetero_flag), ("water", "W"))
        self.assertEqual(rows[rows.source_atom == "O1A"].source_hetero_flag.iloc[0], "H_HEM")

    def test_selector_by_number_and_missing_source_is_a_status_row(self):
        b = heme(Builder())
        b.add("NE2", "HIS", "A", 64, (0.0, 3.0, 0.0))
        sources = [PROPIONATE_A,
                   {"label": "O2_ligand", "resname": "OXY"},
                   {"label": "residue_64", "resname": "HIS", "resnum": 64, "chain": "A", "atoms": ["NE2", "ND1"]},
                   {"label": "residue_64_as_ser", "resname": "SER", "resnum": 64}]
        result = self.run_site(b, sources)
        status = result["site_status"].set_index("source_group")
        self.assertEqual(status.loc["propionate_A", "status"], "matched")
        self.assertEqual(status.loc["O2_ligand", "status"], "not_found")
        self.assertEqual(status.loc["residue_64", "status"], "partial")
        self.assertEqual(status.loc["residue_64", "missing_atoms"], "ND1")
        self.assertEqual(status.loc["residue_64_as_ser", "status"], "not_found")
        self.assertIn("HIS", status.loc["residue_64_as_ser", "note"])
        self.assertFalse((result["site_atom_contacts"].source_group == "O2_ligand").any())
        self.assertIn("residue_64", set(result["site_atom_contacts"].source_group))

    def test_other_source_groups_stay_separate(self):
        b = heme(Builder())
        b.add("NE2", "HIS", "A", 64, (0.0, 3.0, 0.0))                  # near the A propionate only
        b.add("NH2", "ARG", "A", 45, (0.0, 17.0, 0.0))                 # near the D propionate only
        rows = self.run_site(b, [PROPIONATE_A, PROPIONATE_D])["site_atom_contacts"]
        self.assertEqual(set(rows[rows.partner_resname == "HIS"].source_group), {"propionate_A"})
        self.assertEqual(set(rows[rows.partner_resname == "ARG"].source_group), {"propionate_D"})

    def test_residue_with_two_names_at_one_number_does_not_crash(self):
        b = heme(Builder())
        b.add("OG", "SER", "A", 92, (0.0, 3.0, 0.0), altloc="A", occupancy=0.6)
        b.add("SG", "CYS", "A", 92, (0.0, 3.4, 0.0), altloc="B", occupancy=0.4, element="S")
        result = self.run_site(b, [PROPIONATE_A])
        names = set(result["site_atom_contacts"].partner_resname)
        self.assertIn("SER", names)
        self.assertIn("CYS", names)
        self.assertEqual(result["manifest"]["counts"]["n_disordered_residues"], 1)

    # -- water bridges
    def test_single_water_bridge_and_cutoff(self):
        b = heme(Builder())
        b.add("O", "HOH", "A", 301, (3.3, 0.0, 0.0), hetero=True)                     # 3.3 A from O1A
        b.add("NH1", "ARG", "A", 45, (3.3, 3.2, 0.0))                                  # 3.2 A from the water, > 3.5 from O1A
        b.add("OD1", "ASP", "A", 60, (3.3, -4.0, 0.0))                                 # 4.0 A from the water: too far
        bridges = self.run_site(b, [PROPIONATE_A])["site_water_bridges"]
        hit = bridges[(bridges.source_atom == "O1A") & (bridges.partner_resname == "ARG")]
        self.assertEqual(len(hit), 1)
        row = hit.iloc[0]
        self.assertEqual((row.water_resnum, row.source_water_distance_A, row.water_partner_distance_A), (301, 3.3, 3.2))
        self.assertFalse(row.direct_polar_candidate)
        self.assertFalse((bridges.partner_resname == "ASP").any())
        tight = self.run_site(b, [PROPIONATE_A], water_bridge_cutoff=3.0)["site_water_bridges"]
        self.assertEqual(len(tight), 0)

    def test_water_bridge_never_returns_the_sources_own_residue_or_other_waters(self):
        b = heme(Builder())
        b.add("O", "HOH", "A", 301, (3.3, 0.0, 0.0), hetero=True)
        b.add("O", "HOH", "A", 302, (3.3, 3.0, 0.0), hetero=True)
        bridges = self.run_site(b, [PROPIONATE_A])["site_water_bridges"]
        self.assertEqual(len(bridges), 0)

    # -- manifest and outputs
    def test_manifest_records_input_hash_version_parameters_and_caveats(self):
        path = self.two_conformer_his().write(self.folder)
        result = sc.site_contacts(path, [PROPIONATE_A])
        manifest = result["manifest"]
        self.assertEqual(manifest["input"]["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
        self.assertEqual(manifest["input"]["size_bytes"], path.stat().st_size)
        self.assertEqual(set(result["site_atom_contacts"].input_sha256), {manifest["input"]["sha256"]})
        for key in ("package_version", "git_commit", "git_working_tree_clean", "module_sha256", "python", "biopython", "numpy",
                    "pandas", "created_utc", "site_contacts_version", "symmetry", "caveats", "sources"):
            self.assertIn(key, manifest)
        self.assertEqual(manifest["package_version"], "0.2.1")
        self.assertEqual(manifest["parameters"]["neighbor_cutoff_A"], 4.5)
        self.assertEqual(manifest["parameters"]["altloc_policy"], "all")
        self.assertTrue(manifest["parameters"]["heavy_atoms_only"])
        self.assertIn("asymmetric unit only", manifest["symmetry"])
        json.dumps(manifest, default=str)

    def test_write_outputs_lists_file_hashes(self):
        result = self.run_site(self.two_conformer_his(), [PROPIONATE_A])
        out = self.folder / "out" / "entry"
        written = sc.write_site_contacts(result, out, suffix="__all")
        self.assertEqual({p.name for p in written},
                         {"site_atom_contacts__all.csv", "site_residue_summary__all.csv", "site_water_bridges__all.csv",
                          "site_status__all.csv", "site_manifest__all.json"})
        manifest = json.loads((out / "site_manifest__all.json").read_text())
        for name, info in manifest["output_files"].items():
            self.assertEqual(info["sha256"], hashlib.sha256((out / name).read_bytes()).hexdigest())
        again = pd.read_csv(out / "site_atom_contacts__all.csv", keep_default_na=False)
        self.assertEqual(list(again.columns), sc.CONTACT_COLUMNS)
        self.assertEqual({p.name for p in out.iterdir()}, {p.name for p in written})

    def test_inputs_are_not_modified(self):
        path = self.two_conformer_his().write(self.folder)
        before = path.read_bytes()
        sc.site_contacts(path, [PROPIONATE_A])
        self.assertEqual(path.read_bytes(), before)

    def test_bad_arguments(self):
        path = self.two_conformer_his().write(self.folder)
        with self.assertRaises(ValueError):
            sc.site_contacts(path, [PROPIONATE_A], altloc_policy="first")
        with self.assertRaises(ValueError):
            sc.site_contacts(path, [PROPIONATE_A], polar_cutoff=5.0)
        with self.assertRaises(ValueError):
            sc.site_contacts(path, [])
        with self.assertRaises(ValueError):
            sc.site_contacts(path, [{"label": "x"}])
        with self.assertRaises(FileNotFoundError):
            sc.site_contacts(self.folder / "missing.pdb", [PROPIONATE_A])


if __name__ == "__main__":
    unittest.main()
