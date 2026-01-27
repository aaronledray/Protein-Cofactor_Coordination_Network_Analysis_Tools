"""Pytest configuration and fixtures."""

import os
import sys
import pytest

# Add project root to path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)


@pytest.fixture
def project_root():
    """Return the project root directory."""
    return PROJECT_ROOT


@pytest.fixture
def test_data_dir():
    """Return the reference structures directory."""
    return os.path.join(PROJECT_ROOT, "reference_structures")


@pytest.fixture
def plastocyanin_cif(test_data_dir):
    """Return path to plastocyanin CIF file."""
    path = os.path.join(test_data_dir, "0_Plastocyanin", "1ag6.cif")
    if not os.path.exists(path):
        pytest.skip(f"Test file not found: {path}")
    return path


@pytest.fixture
def oex_pdb(test_data_dir):
    """Return path to OEX PDB file."""
    path = os.path.join(test_data_dir, "1_OEX", "0_Aligned_Reduced", "4ub6.pdb")
    if not os.path.exists(path):
        pytest.skip(f"Test file not found: {path}")
    return path
