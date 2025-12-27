# # modules/reporting.py

# from __future__ import annotations
# from typing import List, Dict
# import csv
# import numpy as np

# # Depends on your chemistry utilities for moiety expansion
# from .chemistry import get_moiety_atoms





# def generate_coordination_csv_with_moieties(
#     cofactor_sphere: List[Dict],
#     pcs_residues: List[Dict],
#     scs_residues: List[Dict],
#     pdb_file_name: str,
#     structure,  # Biopython structure object
#     bond_lookup,  # currently unused here but kept for signature stability
#     distance_cutoff: float = 3.0,
#     output_file_name: str | None = None,
# ) -> str:
#     """
#     Generate a CSV summarizing the Cofactor / PCS / SCS residues with:
#       - residue name, number, chain
#       - chemical moiety (expanded by get_moiety_atoms)
#       - intra-network interactors within `distance_cutoff`

#     Inputs are lists of atom dicts with keys:
#       'coordinates', 'name', 'element', 'residue', 'residue_number', 'chain'
#     """
#     if output_file_name is None:
#         output_file_name = f"{pdb_file_name}_Coord_Breakdown.csv"

#     rows: List[Dict] = []

#     def _get_residue_object(atom_dict):
#         """Return the Biopython Residue object matching an atom dict."""
#         target_res = atom_dict["residue"]
#         target_num = atom_dict["residue_number"]
#         target_chain = atom_dict["chain"]
#         for res in structure.get_residues():
#             if (
#                 res.get_resname() == target_res
#                 and res.get_id()[1] == target_num
#                 and res.get_full_id()[2] == target_chain
#             ):
#                 return res
#         return None

#     def _nearest_interactors(atom_dict, search_atoms):
#         """Names of atoms within cutoff in the provided list (excluding same coordinate object)."""
#         hits = []
#         a = np.array(atom_dict["coordinates"], dtype=float)
#         for other in search_atoms:
#             if other is atom_dict:
#                 continue
#             b = np.array(other["coordinates"], dtype=float)
#             if np.linalg.norm(a - b) <= distance_cutoff:
#                 hits.append(other["name"])
#         return hits

#     def _collect(residue_atoms: List[Dict], label: str):
#         """
#         Collapse per-atom dicts into per-residue rows with:
#           - residue info
#           - combined moiety (name list)
#           - combined interactors (names list)
#         """
#         by_key: Dict[tuple, Dict] = {}
#         for atom in residue_atoms:
#             key = (atom["residue"], atom["residue_number"], atom["chain"])
#             bucket = by_key.setdefault(
#                 key,
#                 {
#                     "residue_name": atom["residue"],
#                     "residue_number": atom["residue_number"],
#                     "chain": atom["chain"],
#                     "moiety": set(),
#                     "interactors": set(),
#                 },
#             )

#             # Find residue object and expand the atom to its moiety members
#             residue_obj = _get_residue_object(atom)
#             if residue_obj is not None:
#                 try:
#                     # find the matching Biopython Atom inside this residue
#                     biopy_atom = next(a for a in residue_obj if a.get_name() == atom["name"])
#                     moiety_atoms = get_moiety_atoms(biopy_atom, residue_obj)
#                     bucket["moiety"].update(a.get_name() for a in moiety_atoms)
#                 except StopIteration:
#                     # fallback: no matching atom by name in residue
#                     pass

#             # Interactors against PCS+SCS (you can include cofactor too if desired)
#             bucket["interactors"].update(
#                 _nearest_interactors(atom, pcs_residues + scs_residues)
#             )

#         # Emit rows
#         for (_, resnum, chain), info in by_key.items():
#             rows.append(
#                 {
#                     "Residue Category": label,
#                     "Residue Name": info["residue_name"],
#                     "Residue Number": resnum,
#                     "Chain": chain,
#                     "Chemical Moiety": ", ".join(sorted(info["moiety"])) if info["moiety"] else "",
#                     "Interactors": ", ".join(sorted(info["interactors"])) if info["interactors"] else "",
#                 }
#             )

#     # Build CSV rows: Cofactor, PCS, SCS
#     _collect(cofactor_sphere, "Cofactor")
#     _collect(pcs_residues, "PCS")
#     _collect(scs_residues, "SCS")

#     # Write CSV
#     with open(output_file_name, "w", newline="") as f:
#         writer = csv.DictWriter(
#             f,
#             fieldnames=[
#                 "Residue Category",
#                 "Residue Name",
#                 "Residue Number",
#                 "Chain",
#                 "Chemical Moiety",
#                 "Interactors",
#             ],
#         )
#         writer.writeheader()
#         writer.writerows(rows)

#     print(f"CSV file saved as {output_file_name}")
#     return output_file_name



# __all__ = ["generate_coordination_csv_with_moieties"]



from typing import List, Dict, Any, Tuple, Iterable, Optional
import csv
import os
import numpy as np

def _as_key(a: Dict[str, Any]) -> Tuple[str, str, str]:
    """Unique atom key as (chain, residue_number(str), atom_name). Robust for PDB/mmCIF."""
    chain = str(a.get("chain", "") or "")
    # residue_number can be int or '123A' etc; stringify consistently
    resnum = str(a.get("residue_number", "") or "")
    atom  = str(a.get("name", "") or "")
    return (chain, resnum, atom)

def _res_key(a: Dict[str, Any]) -> Tuple[str, str]:
    """Residue key as (chain, residue_number str)."""
    chain = str(a.get("chain", "") or "")
    resnum = str(a.get("residue_number", "") or "")
    return (chain, resnum)

def _collect_interactor_sets(
    cofactor_atoms: List[Dict],
    pcs_atoms: List[Dict],
    scs_atoms: List[Dict],
):
    """
    Build:
      - atom-level sets: keys are (chain,resnum,atom)
      - residue-level sets: keys are (chain,resnum)
      - a per-residue category map (Cofactor, PCS, SCS) with PCS priority over SCS if overlap.
    """
    cof_atom_set = {_as_key(a) for a in (cofactor_atoms or [])}
    pcs_atom_set = {_as_key(a) for a in (pcs_atoms or [])}
    scs_atom_set = {_as_key(a) for a in (scs_atoms or [])}

    cof_res_set = {_res_key(a) for a in (cofactor_atoms or [])}
    pcs_res_set = {_res_key(a) for a in (pcs_atoms or [])}
    scs_res_set = {_res_key(a) for a in (scs_atoms or [])}

    # Category per residue; if something is both PCS and SCS, call it PCS.
    res_category: Dict[Tuple[str, str], str] = {}
    for rk in scs_res_set:
        res_category[rk] = "SCS"
    for rk in pcs_res_set:
        res_category[rk] = "PCS"
    for rk in cof_res_set:
        res_category[rk] = "Cofactor"

    return {
        "cof_atom_set": cof_atom_set,
        "pcs_atom_set": pcs_atom_set,
        "scs_atom_set": scs_atom_set,
        "res_category": res_category,
    }



def _moiety_name(resname: str, atom_name: str, moiety_lookup: Dict[Tuple[str, str], str]) -> str:
    """Lookup moiety label with safe fallback."""
    key = (str(resname or "").upper(), str(atom_name or "").upper())
    return moiety_lookup.get(key, "unknown_moiety")

def write_coord_breakdown_v2(
    structure,                                  # Bio.PDB Structure
    cofactor_atoms: List[Dict],
    pcs_atoms: List[Dict],
    scs_atoms: List[Dict],
    moiety_lookup: Dict[Tuple[str, str], str],  # your moiety table (incl. HEM/ICS/HCA/etc.)
    output_csv_path: str = "Coord_Breakdown.csv",
):
    """
    Writes a per-atom CSV with columns:
      Category, Atom, Residue, Residue Number, Chain, Moiety, InteractorFlag, x, y, z

    Behavior:
      - Any residue that appears in cofactor_atoms / pcs_atoms / scs_atoms is included in full
        (all atoms of that residue are written).
      - Atoms that are explicitly present in those lists are flagged as 'interactor';
        others in the same residue are 'non-interactor'.
      - Category for each residue is Cofactor / PCS / SCS (PCS wins over SCS if overlapping).
      - Moiety is taken from `moiety_lookup[(RESNAME, ATOMNAME)]` with fallback 'unknown_moiety'.
    """
    coll = _collect_interactor_sets(cofactor_atoms, pcs_atoms, scs_atoms)
    cof_atom_set = coll["cof_atom_set"]
    pcs_atom_set = coll["pcs_atom_set"]
    scs_atom_set = coll["scs_atom_set"]
    res_category = coll["res_category"]

    # Convenience: build a quick index of residues we care about.
    target_residues = set(res_category.keys())  # (chain, resnum_str)

    # Prepare rows
    rows: List[List[Any]] = []
    header = [
        "Category",
        "Atom",
        "Residue",
        "Residue Number",
        "Chain",
        "Moiety",
        "InteractorFlag",
        "x", "y", "z",
    ]

    # Iterate structure and dump atoms for residues we care about
    for model in structure:
        for chain in model:
            chain_id = str(chain.id)
            for residue in chain:
                # residue number can be int or str; normalize to str
                try:
                    resnum = residue.get_id()[1]
                except Exception:
                    resnum = ""
                resnum_str = str(resnum)
                rk = (chain_id, resnum_str)
                if rk not in target_residues:
                    continue  # only residues that participate in interactions

                category = res_category.get(rk, "SCS")  # default SCS if somehow missing
                resname = residue.get_resname()

                for atom in residue:
                    atom_name = atom.get_name()
                    atom_key = (chain_id, resnum_str, atom_name)
                    # determine interactor flag by membership in the provided atom sets
                    if atom_key in cof_atom_set or atom_key in pcs_atom_set or atom_key in scs_atom_set:
                        flag = "interactor"
                    else:
                        flag = "non-interactor"

                    coords = np.asarray(atom.coord, dtype=float).tolist()
                    moiety = _moiety_name(resname, atom_name, moiety_lookup)

                    rows.append([
                        category,
                        atom_name,
                        resname,
                        resnum_str,
                        chain_id,
                        moiety,
                        flag,
                        coords[0], coords[1], coords[2],
                    ])

    # Write CSV
    os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)) or ".", exist_ok=True)
    with open(output_csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

    print(f"[INFO] Wrote {len(rows)} rows to '{output_csv_path}'")
