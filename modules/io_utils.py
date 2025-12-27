# modules/io_utils.py
import os
from typing import Tuple, List, Dict, Any
import numpy as np

def _choose_parser(struct_path: str):
    """Return a Biopython parser appropriate for the structure file."""
    from Bio.PDB import PDBParser, MMCIFParser
    lp = struct_path.lower()
    if lp.endswith(".pdb") or lp.endswith(".pdb.gz"):
        return PDBParser(QUIET=True)
    if lp.endswith(".cif") or lp.endswith(".mmcif") or lp.endswith(".cif.gz") or lp.endswith(".mmcif.gz"):
        return MMCIFParser(QUIET=True)
    # Fallback: try CIF first, then PDB
    try:
        return MMCIFParser(QUIET=True)
    except Exception:
        return PDBParser(QUIET=True)

def _safe_resnum(residue):
    """Return an int residue number when possible, else the original value."""
    try:
        resnum = residue.get_id()[1]
    except Exception:
        return None
    # Already an int
    if isinstance(resnum, (int, np.integer)):
        return int(resnum)
    # Try to coerce strings like "123", leave insertion codes like "123A" as-is
    try:
        return int(str(resnum).strip())
    except Exception:
        return str(resnum)

def unpack_pdb_file(struct_path: str) -> Tuple[Any, List[Dict]]:
    """
    Load a PDB or mmCIF file and return:
      (structure, atoms_data)
    atoms_data: list of dicts with keys:
      name, element, coordinates (np.array shape (3,)),
      residue (resname), residue_number, chain
    """
    struct_path = os.path.abspath(os.path.expanduser(struct_path))
    parser = _choose_parser(struct_path)
    structure = parser.get_structure("structure", struct_path)

    atoms_data: List[Dict] = []
    for model in structure:
        for chain in model:
            chain_id = chain.id
            for residue in chain:
                resname = residue.get_resname()
                resnum = _safe_resnum(residue)
                for atom in residue:
                    atoms_data.append({
                        "name": atom.get_name(),
                        "element": getattr(atom, "element", "") or "",   # mmCIF/PDB both supported
                        "coordinates": np.array(atom.coord, dtype=float),
                        "residue": resname,
                        "residue_number": resnum,
                        "chain": chain_id,
                    })

    return structure, atoms_data
