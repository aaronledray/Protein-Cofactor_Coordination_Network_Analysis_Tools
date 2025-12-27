# # modules/moieties.py

# from __future__ import annotations



# # CHEMICAL MOIETIES LIVE IN CHEMISTRY MODULE
# from .chemistry import chemical_moieties


# try:
#     from .chemistry import bond_lookup  # noqa: F401
# # except ImportError:
# #     # Minimal, safe default you can expand later
# #     bond_lookup = {
# #         # Peptide backbone (very small starter set)
# #         "ALA": [("N", "CA"), ("CA", "C"), ("C", "O")],
# #         "GLY": [("N", "CA"), ("CA", "C"), ("C", "O")],
# #         # Add more residues / cofactors as needed…
# #     }






# # Canonical element colors used across plotting
# atom_type_colors = {
#     "C": "black",
#     "N": "blue",
#     "O": "red",
#     "S": "yellow",
#     "FE": "orange",
#     "MN": "purple",
#     "CA": "green",
#     "CU": "goldenrod",
#     "MO": "teal",
# }

# __all__ = ["chemical_moieties", "bond_lookup", "atom_type_colors"]





# modules/moieties.py
from __future__ import annotations

from typing import Dict, List, Tuple

# ----------------------------
# Optional imports from chemistry
# ----------------------------
chemical_moieties: Dict[tuple, str]
bond_lookup: Dict[str, List[Tuple[str, str]]]

try:
    # CHEMICAL MOIETIES LIVE IN CHEMISTRY MODULE
    from .chemistry import chemical_moieties as _chem_moieties  # type: ignore
    chemical_moieties = _chem_moieties
except Exception:
    # Safe default if chemistry module or variable missing
    chemical_moieties = {}

# bond_lookup is optional; provide defaults if absent
try:
    from .chemistry import bond_lookup as _bond_lookup  # type: ignore
    bond_lookup = _bond_lookup
except Exception:
    # Minimal, safe default you can expand later
    bond_lookup = {}

# ----------------------------
# Canonical element colors used across plotting
# ----------------------------
atom_type_colors: Dict[str, str] = {
    "C": "black",
    "N": "blue",
    "O": "red",
    "S": "yellow",
    "FE": "orange",
    "MN": "purple",
    "CA": "green",
    "CU": "goldenrod",
    "MO": "teal",
}

__all__ = ["chemical_moieties", "bond_lookup", "atom_type_colors"]
