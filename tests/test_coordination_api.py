"""Tests for the side-effect-free and batch coordination APIs."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from modules.coordination_api import _motif_contact_rows, analyze_structure, batch_analyze


ROOT = Path(__file__).resolve().parents[1]
PLASTOCYANIN = ROOT / "reference_structures/0_Plastocyanin/1ag6.cif"
OEX = ROOT / "reference_structures/1_OEX/0_Aligned_Reduced/4ub6.pdb"
ZINC = ROOT / "reference_structures/3_Metals/1ca2.pdb"
MYOGLOBIN = ROOT / "reference_structures/3_Myoglobin/1a6m.pdb"


class CoordinationApiTests(unittest.TestCase):
    def test_analyze_structure_returns_tables_without_files(self):
        with TemporaryDirectory(prefix="sscna-api-") as work_dir:
            with redirect_stdout(StringIO()):
                result = analyze_structure(
                    PLASTOCYANIN,
                    "CU",
                    exclude_moieties=["alanine_sidechain"],
                )

            self.assertEqual(set(result), {"residues", "atoms", "links", "contacts"})
            self.assertEqual(len(result["residues"]), 8)
            self.assertEqual(len(result["atoms"]), 8)
            self.assertEqual(len(result["links"]), 7)
            self.assertGreaterEqual(len(result["contacts"]), len(result["links"]))
            self.assertIn("motif", result["atoms"].columns)
            self.assertIn("coordination_role", result["atoms"].columns)
            self.assertIn("src_motif", result["contacts"].columns)
            self.assertIn("dst_motif", result["contacts"].columns)
            self.assertEqual(
                list(result["residues"].columns),
                [
                    "structure_id", "site_id", "cofactor_residue_name",
                    "cofactor_residue_number", "cofactor_chain",
                    "cofactor_insertion_code", "shell", "residue_name",
                    "residue_number", "chain", "insertion_code", "hetero_flag",
                    "atoms_involved", "minimum_distance_A",
                ],
            )
            self.assertEqual(
                set(result["atoms"]["coordination_role"]),
                {"primary_coordinator", "network_context", "active_site_component"},
            )
            self.assertEqual(result["residues"]["site_id"].unique().tolist(), ["all"])
            self.assertEqual(
                result["residues"]
                .loc[result["residues"]["shell"] == "PCS", "minimum_distance_A"]
                .isna()
                .sum(),
                0,
            )
            self.assertEqual(list(Path(work_dir).iterdir()), [])

    def test_batch_analyze_isolates_bad_inputs_and_combines_results(self):
        result = batch_analyze(
            [PLASTOCYANIN, OEX, ROOT / "does-not-exist.pdb"],
            "CU,OEX",
            workers=2,
            exclude_moieties=["alanine_sidechain"],
        )

        self.assertEqual(len(result["errors"]), 1)
        self.assertEqual(
            set(result["residues"]["structure_id"]),
            {"1ag6", "4ub6"},
        )
        self.assertEqual(set(result["links"]["structure_id"]), {"1ag6", "4ub6"})
        self.assertEqual(set(result["contacts"]["structure_id"]), {"1ag6", "4ub6"})

    def test_motif_contacts_preserve_atom_and_motif_identity(self):
        with redirect_stdout(StringIO()):
            result = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
            )

        copper_contacts = result["contacts"].query("src_resname == 'CU' and dst_resname == 'HIS'")
        self.assertGreaterEqual(len(copper_contacts), 2)
        self.assertTrue(set(copper_contacts["dst_atom"]).issubset({"ND1", "NE2"}))
        self.assertEqual(set(copper_contacts["dst_motif"]), {"imidazole"})
        self.assertTrue(copper_contacts["distance_A"].le(3.6).all())

        with redirect_stdout(StringIO()):
            zinc_contacts = analyze_structure(ZINC, "ZN")["contacts"]
        direct_histidines = zinc_contacts.query(
            "src_resname == 'ZN' and dst_resname == 'HIS' and direct_coordination == True"
        )
        self.assertEqual(len(direct_histidines), 3)
        self.assertTrue(set(direct_histidines["dst_atom"]).issubset({"ND1", "NE2"}))
        self.assertEqual(set(direct_histidines["dst_motif"]), {"imidazole"})
        secondary_carboxylates = zinc_contacts.query(
            "src_resname == 'HIS' and dst_resname == 'GLU'"
        )
        self.assertGreaterEqual(len(secondary_carboxylates), 1)
        self.assertTrue(set(secondary_carboxylates["dst_atom"]).issubset({"OE1", "OE2"}))
        self.assertEqual(set(secondary_carboxylates["dst_motif"]), {"COO"})

    def test_shells_three_adds_tcs_and_adjacent_links(self):
        with redirect_stdout(StringIO()):
            result = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
                shells=3,
            )

        self.assertIn("TCS", set(result["residues"]["shell"]))
        self.assertIn("scs->tcs", set(result["links"]["link_type"]))

    def test_heme_propionate_water_contact_is_inferred_primary_hbond(self):
        with redirect_stdout(StringIO()):
            result = analyze_structure(
                MYOGLOBIN,
                "HEM",
                exclude_moieties=["alanine_sidechain"],
                shells=3,
            )

        contacts = result["contacts"]
        inferred = contacts.query(
            "contact_role == 'primary_motif_contact' and "
            "src_resname == 'HEM' and dst_resname == 'HOH'"
        )
        self.assertGreaterEqual(len(inferred), 1)
        self.assertEqual(set(inferred["src_motif"]), {"heme_propionate"})
        self.assertEqual(set(inferred["dst_motif"]), {"water"})
        self.assertTrue(inferred["distance_A"].le(3.2).all())

    def test_heme_propionate_hydroxyl_is_primary_but_not_direct(self):
        cofactor = [{
            "residue": "HEM",
            "residue_number": 1500,
            "chain": "A",
            "name": "O1A",
            "element": "O",
            "model_id": 0,
            "coordinates": [0.0, 0.0, 0.0],
        }]
        shell_atoms = {1: [{
            "residue": "SER",
            "residue_number": 315,
            "chain": "A",
            "name": "OG",
            "element": "O",
            "model_id": 0,
            "coordinates": [2.55, 0.0, 0.0],
        }]}

        contacts = _motif_contact_rows(
            "katg",
            "site_1",
            cofactor,
            shell_atoms,
            distance_cutoff=3.6,
            direct_coordination_cutoff=2.6,
        )

        self.assertEqual(len(contacts), 1)
        self.assertEqual(contacts.iloc[0]["contact_role"], "primary_motif_contact")
        self.assertFalse(bool(contacts.iloc[0]["direct_coordination"]))
        self.assertEqual(contacts.iloc[0]["dst_motif"], "hydroxyl")

    def test_per_site_mode_assigns_stable_site_ids(self):
        with redirect_stdout(StringIO()):
            result = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
                site_mode="per-site",
            )

        self.assertEqual(set(result["residues"]["site_id"]), {"site_1"})
        self.assertEqual(set(result["atoms"]["site_id"]), {"site_1"})
        self.assertIn("site_id", result["links"].columns)

    def test_definition_options_are_opt_in(self):
        with redirect_stdout(StringIO()):
            default = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
            )
            carbon = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
                include_carbon_seeds=True,
            )
            direct = analyze_structure(
                PLASTOCYANIN,
                "CU",
                exclude_moieties=["alanine_sidechain"],
                direct_coordination=True,
            )

        default_pcs = len(default["atoms"].query("shell == 'PCS'"))
        carbon_pcs = len(carbon["atoms"].query("shell == 'PCS'"))
        self.assertEqual(default_pcs, 5)
        self.assertGreater(carbon_pcs, default_pcs)
        self.assertEqual(int(default["links"]["direct_coordination"].sum()), 0)
        self.assertEqual(int(direct["links"]["direct_coordination"].sum()), 3)


if __name__ == "__main__":
    unittest.main()
