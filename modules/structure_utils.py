# modules/structure_utils.py
"""
Utilities for working with Biopython structures and residues.
"""

from typing import List, Dict
from Bio.PDB.Residue import Residue

def get_residue_atoms(residue: Residue) -> List[Dict]:
    """
    Return all atoms in a residue as a list of dictionaries with atom details.
    Keys: 'coordinates', 'name', 'element', 'residue', 'residue_number', 'chain'
    """
    residue_number = residue.get_id()[1]          # integer residue number
    residue_name = residue.get_resname()          # e.g., "HEM", "HIS"
    chain_id = residue.get_full_id()[2]           # chain identifier

    atoms: List[Dict] = []
    for atom in residue:
        atoms.append({
            "coordinates": atom.coord,
            "name": atom.get_name(),
            "element": atom.element,
            "residue": residue_name,
            "residue_number": residue_number,
            "chain": chain_id,
        })
    return atoms



from typing import Iterable, List, Union
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom

def _normalize_resname_list(names: Union[str, Iterable[str]]) -> List[str]:
    """
    Normalize residue-name input to a list of uppercase strings with whitespace trimmed.
    """
    if isinstance(names, str):
        names = [names]
    return [str(n).strip().upper() for n in names]

def find_cofactor_atoms(
    structure: Structure,
    cofactor_resname: Union[str, Iterable[str]],
    alternative_resname: Union[None, str, Iterable[str]] = None,
) -> List[Atom]:
    """
    Return all Atom objects whose parent residue name matches any of the provided
    cofactor residue names (and optional alternative names).

    Parameters
    ----------
    structure : Biopython Structure
        Parsed structure from Bio.PDB (e.g., via PDBParser()).
    cofactor_resname : str | Iterable[str]
        Primary residue name(s) of the cofactor (e.g., "HEM" or ["HEA","HEM","HM1","HM2"]).
    alternative_resname : None | str | Iterable[str], optional
        Additional residue name(s) to accept.

    Returns
    -------
    List[Atom]
        List of Biopython Atom objects belonging to matching residues.
    """
    targets = set(_normalize_resname_list(cofactor_resname))
    if alternative_resname is not None:
        targets.update(_normalize_resname_list(alternative_resname))

    out: List[Atom] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip().upper()
                if resname in targets:
                    out.extend(list(residue.get_atoms()))
    return out
