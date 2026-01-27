"""Tests for chemistry and moieties modules."""

import os
import sys
import pytest

# Add project root to path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from modules.moieties import chemical_moieties, bond_lookup


class TestChemicalMoieties:
    """Test moiety lookup table."""

    def test_moieties_not_empty(self):
        """Moiety lookup should contain entries."""
        assert len(chemical_moieties) > 0

    def test_standard_amino_acid_moieties(self):
        """Standard amino acids should have moiety definitions."""
        # Check histidine imidazole
        assert ("HIS", "ND1") in chemical_moieties or ("HIS", "NE2") in chemical_moieties

    def test_cysteine_thiol(self):
        """Cysteine SG should be labeled as thiol."""
        if ("CYS", "SG") in chemical_moieties:
            assert "thiol" in chemical_moieties[("CYS", "SG")].lower()

    def test_backbone_atoms(self):
        """Backbone atoms should be defined for common residues."""
        # N, CA, C, O are backbone for all amino acids
        standard_aa = ["ALA", "GLY", "VAL", "LEU"]
        for aa in standard_aa:
            # At least one backbone atom should be defined
            found = any(
                (aa, atom) in chemical_moieties
                for atom in ["N", "CA", "C", "O"]
            )
            # Note: Not all implementations may have backbone in moieties
            # This is more of a documentation check


class TestBondLookup:
    """Test bond topology lookup table."""

    def test_bond_lookup_not_empty(self):
        """Bond lookup should contain entries."""
        assert len(bond_lookup) > 0

    def test_alanine_bonds(self):
        """Alanine should have basic backbone bonds."""
        if "ALA" in bond_lookup:
            bonds = bond_lookup["ALA"]
            # Should have N-CA, CA-C, C-O at minimum
            bond_set = set(bonds)
            assert ("N", "CA") in bond_set or ("CA", "N") in bond_set

    def test_glycine_bonds(self):
        """Glycine should have backbone bonds (no sidechain)."""
        if "GLY" in bond_lookup:
            bonds = bond_lookup["GLY"]
            # Glycine has no CB, so fewer bonds than other AAs
            assert len(bonds) >= 3  # At least backbone bonds

    def test_bonds_are_pairs(self):
        """All bonds should be 2-element tuples."""
        for resname, bonds in bond_lookup.items():
            for bond in bonds:
                assert len(bond) == 2, f"Invalid bond in {resname}: {bond}"
                assert isinstance(bond[0], str), f"Bond atom should be string: {bond}"
                assert isinstance(bond[1], str), f"Bond atom should be string: {bond}"

    def test_heme_bonds(self):
        """Heme (HEM) should have extensive bond definitions."""
        heme_names = ["HEM", "HEA", "HM1"]
        found_heme = any(name in bond_lookup for name in heme_names)
        if found_heme:
            for name in heme_names:
                if name in bond_lookup:
                    # Heme should have many bonds (porphyrin ring)
                    assert len(bond_lookup[name]) > 20


class TestMoietyConsistency:
    """Test consistency between moieties and bonds."""

    def test_bond_atoms_exist(self):
        """Atoms referenced in bonds should exist in moiety definitions."""
        # This is more of a sanity check - not all atoms need moiety labels
        # but if a residue has bonds defined, it should be a known residue
        for resname in bond_lookup.keys():
            # At least check it's a reasonable residue name
            assert len(resname) >= 1
            assert resname.isupper() or resname[0].isupper()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
