"""Tests for io_utils module."""

import os
import sys
import pytest

# Add project root to path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from modules.io_utils import unpack_pdb_file

# Test data paths (relative to project root)
TEST_DATA_DIR = os.path.join(PROJECT_ROOT, "reference_structures")
PLASTOCYANIN_CIF = os.path.join(TEST_DATA_DIR, "0_Plastocyanin", "1ag6.cif")
OEX_PDB = os.path.join(TEST_DATA_DIR, "1_OEX", "0_Aligned_Reduced", "4ub6.pdb")


class TestUnpackPdbFile:
    """Test structure file loading."""

    @pytest.mark.skipif(
        not os.path.exists(PLASTOCYANIN_CIF),
        reason="Test data file not found"
    )
    def test_load_cif_file(self):
        """Should load mmCIF file and return structure + atoms."""
        structure, atoms = unpack_pdb_file(PLASTOCYANIN_CIF)
        assert structure is not None
        assert len(atoms) > 0
        # Check atom dict structure
        first_atom = atoms[0]
        assert "coordinates" in first_atom
        assert "name" in first_atom
        assert "element" in first_atom
        assert "residue" in first_atom
        assert "residue_number" in first_atom
        assert "chain" in first_atom

    @pytest.mark.skipif(
        not os.path.exists(OEX_PDB),
        reason="Test data file not found"
    )
    def test_load_pdb_file(self):
        """Should load PDB file and return structure + atoms."""
        structure, atoms = unpack_pdb_file(OEX_PDB)
        assert structure is not None
        assert len(atoms) > 0

    @pytest.mark.skipif(
        not os.path.exists(PLASTOCYANIN_CIF),
        reason="Test data file not found"
    )
    def test_atoms_have_coordinates(self):
        """All atoms should have 3D coordinates."""
        _, atoms = unpack_pdb_file(PLASTOCYANIN_CIF)
        for atom in atoms[:10]:  # Check first 10
            coords = atom["coordinates"]
            assert len(coords) == 3
            assert all(isinstance(c, (int, float)) for c in coords)

    def test_nonexistent_file_raises(self):
        """Should raise error for nonexistent file."""
        with pytest.raises(Exception):
            unpack_pdb_file("/nonexistent/path/to/file.pdb")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
