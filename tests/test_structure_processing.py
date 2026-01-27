"""Tests for structure_processing module."""

import os
import sys
import numpy as np
import pytest

# Add project root to path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from modules.structure_processing import (
    make_residue_centroid_sphere,
    extract_query_box,
    generate_residue_bonds,
    _euclid2,
    _nearest_idx_and_dist,
)


class TestEuclideanDistance:
    """Test distance calculation utilities."""

    def test_euclid2_same_point(self):
        """Distance squared between same point should be 0."""
        pt = (1.0, 2.0, 3.0)
        assert _euclid2(pt, pt) == 0.0

    def test_euclid2_unit_distance(self):
        """Test unit distance along axis."""
        a = (0.0, 0.0, 0.0)
        b = (1.0, 0.0, 0.0)
        assert _euclid2(a, b) == 1.0

    def test_euclid2_diagonal(self):
        """Test 3D diagonal distance."""
        a = (0.0, 0.0, 0.0)
        b = (1.0, 1.0, 1.0)
        assert abs(_euclid2(a, b) - 3.0) < 1e-10


class TestNearestIdxAndDist:
    """Test nearest point finding."""

    def test_nearest_to_empty_cloud(self):
        """Should return -1 and inf for empty cloud."""
        pt = np.array([0.0, 0.0, 0.0])
        cloud = np.empty((0, 3))
        idx, dist = _nearest_idx_and_dist(pt, cloud)
        assert idx == -1
        assert dist == float("inf")

    def test_nearest_single_point(self):
        """Should find the only point in cloud."""
        pt = np.array([1.0, 0.0, 0.0])
        cloud = np.array([[0.0, 0.0, 0.0]])
        idx, dist = _nearest_idx_and_dist(pt, cloud)
        assert idx == 0
        assert abs(dist - 1.0) < 1e-10

    def test_nearest_multiple_points(self):
        """Should find closest point among multiple."""
        pt = np.array([0.0, 0.0, 0.0])
        cloud = np.array([
            [10.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
        ])
        idx, dist = _nearest_idx_and_dist(pt, cloud)
        assert idx == 1
        assert abs(dist - 2.0) < 1e-10


class TestMakeResidueCentroidSphere:
    """Test residue centroid calculation."""

    def test_empty_input(self):
        """Empty input should return empty list."""
        result = make_residue_centroid_sphere([])
        assert result == []

    def test_single_atom(self):
        """Single atom should return that atom's coordinates as centroid."""
        atoms = [
            {
                "residue_number": 1,
                "residue": "ALA",
                "coordinates": np.array([1.0, 2.0, 3.0]),
            }
        ]
        result = make_residue_centroid_sphere(atoms)
        assert len(result) == 1
        assert result[0]["residue_number"] == 1
        assert result[0]["residue_name"] == "ALA"
        np.testing.assert_array_almost_equal(
            result[0]["coordinates"], [1.0, 2.0, 3.0]
        )

    def test_multiple_atoms_same_residue(self):
        """Multiple atoms should be averaged to centroid."""
        atoms = [
            {"residue_number": 1, "residue": "ALA", "coordinates": np.array([0.0, 0.0, 0.0])},
            {"residue_number": 1, "residue": "ALA", "coordinates": np.array([2.0, 0.0, 0.0])},
            {"residue_number": 1, "residue": "ALA", "coordinates": np.array([1.0, 3.0, 0.0])},
        ]
        result = make_residue_centroid_sphere(atoms)
        assert len(result) == 1
        np.testing.assert_array_almost_equal(
            result[0]["coordinates"], [1.0, 1.0, 0.0]
        )

    def test_malformed_entries_skipped(self):
        """Entries missing required fields should be skipped."""
        atoms = [
            {"residue_number": 1, "residue": "ALA", "coordinates": np.array([1.0, 2.0, 3.0])},
            {"residue_number": None, "residue": "GLY", "coordinates": np.array([4.0, 5.0, 6.0])},
            {"residue_number": 2, "residue": None, "coordinates": np.array([7.0, 8.0, 9.0])},
            {"residue_number": 3, "residue": "VAL"},  # missing coordinates
        ]
        result = make_residue_centroid_sphere(atoms)
        assert len(result) == 1
        assert result[0]["residue_number"] == 1


class TestExtractQueryBox:
    """Test query box extraction."""

    def test_empty_template(self):
        """Empty template should use distance_cutoff as cube limit."""
        template = []
        query = [
            {"name": "CA", "coordinates": np.array([0.0, 0.0, 0.0])},
        ]
        result, limit = extract_query_box(template, query, distance_cutoff=5.0)
        assert limit == 5.0

    def test_filters_to_ca_cb(self):
        """Should only return CA and CB atoms."""
        template = [
            {"name": "CA", "coordinates": np.array([0.0, 0.0, 0.0])},
        ]
        query = [
            {"name": "CA", "coordinates": np.array([1.0, 0.0, 0.0])},
            {"name": "CB", "coordinates": np.array([1.0, 1.0, 0.0])},
            {"name": "N", "coordinates": np.array([0.5, 0.5, 0.0])},
            {"name": "O", "coordinates": np.array([0.5, 0.0, 0.5])},
        ]
        result, _ = extract_query_box(template, query, distance_cutoff=5.0)
        names = {a["name"] for a in result}
        assert names == {"CA", "CB"}

    def test_cube_filtering(self):
        """Atoms outside cube should be excluded."""
        template = [
            {"name": "CA", "coordinates": np.array([0.0, 0.0, 0.0])},
        ]
        query = [
            {"name": "CA", "coordinates": np.array([1.0, 0.0, 0.0])},  # inside
            {"name": "CA", "coordinates": np.array([100.0, 0.0, 0.0])},  # outside
        ]
        result, limit = extract_query_box(template, query, distance_cutoff=2.0)
        assert len(result) == 1
        np.testing.assert_array_almost_equal(
            result[0]["coordinates"], [1.0, 0.0, 0.0]
        )


class TestGenerateResidueBonds:
    """Test bond generation from lookup table."""

    def test_empty_residues(self):
        """Empty residues should return empty bonds."""
        result = generate_residue_bonds([], {"ALA": [("N", "CA")]})
        assert result == []

    def test_unknown_residue(self):
        """Unknown residue type should return empty bonds."""
        residues = [
            {"name": "N", "residue": "XXX", "residue_number": 1, "chain": "A", "coordinates": np.array([0.0, 0.0, 0.0])},
            {"name": "CA", "residue": "XXX", "residue_number": 1, "chain": "A", "coordinates": np.array([1.0, 0.0, 0.0])},
        ]
        result = generate_residue_bonds(residues, {"ALA": [("N", "CA")]})
        assert result == []

    def test_valid_bond(self):
        """Should generate bonds for known residue type."""
        residues = [
            {"name": "N", "residue": "ALA", "residue_number": 1, "chain": "A", "coordinates": np.array([0.0, 0.0, 0.0])},
            {"name": "CA", "residue": "ALA", "residue_number": 1, "chain": "A", "coordinates": np.array([1.0, 0.0, 0.0])},
        ]
        bond_lookup = {"ALA": [("N", "CA")]}
        result = generate_residue_bonds(residues, bond_lookup)
        assert len(result) == 1
        assert result[0] == ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))

    def test_missing_atom_for_bond(self):
        """Missing atom should skip that bond."""
        residues = [
            {"name": "N", "residue": "ALA", "residue_number": 1, "chain": "A", "coordinates": np.array([0.0, 0.0, 0.0])},
        ]
        bond_lookup = {"ALA": [("N", "CA"), ("CA", "C")]}
        result = generate_residue_bonds(residues, bond_lookup)
        assert result == []  # CA is missing, so no bonds


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
