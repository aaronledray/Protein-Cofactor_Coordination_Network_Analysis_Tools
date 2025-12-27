# import numpy as np
# import pandas as pd
# import os


from typing import List, Dict, Tuple, Union



# import numpy as np
# from Bio.PDB import NeighborSearch



# # local helpers
# from .chemistry import get_moiety_atoms          # you placed this earlier
# from .structure_utils import get_residue_atoms   # returns list[dict] for a residue
# from .structure_utils import find_cofactor_atoms # finds atoms for given residue name(s)
from .structure_utils import find_cofactor_atoms

import csv
import os
from typing import Dict, List, Tuple, Iterable, Set, Optional
from collections import defaultdict, namedtuple
import numpy as np




import math

from typing import Dict, List, Tuple, Iterable, Set, Optional
from collections import defaultdict, namedtuple
import numpy as np

# Types
AtomDict = Dict[str, object]  # keys: name,residue,residue_number,chain,element,coordinates




# ---- Config: exclude certain elements from *seed coordinators* (not from expansion) ----
EXCLUDED_COORD_ELEMENTS: Set[str] = {"C"}   # <- your requested change



# ---- Multi-coordinator policy (per residue + moiety) ----
# For ASN/GLN amide: allow up to 2 seeds; for ARG guanidinium: up to 3.
# (If you later want to constrain by element, you can add a filter in _multi_pick_for_group().)
MULTI_COORD_POLICY: Dict[Tuple[str, str], int] = {
    ("ASN", "amide"): 2,
    ("GLN", "amide"): 2,
    ("ARG", "guanidinium"): 3,
    ("ASP", "COO"): 2,
    ("GLU", "COO"): 2,
    ("ANY", "backbone"): 2,     # e.g., N and O can both seed (C is excluded elsewhere)
    ("ANY", "c_terminus"): 2,   # e.g., O and OXT can both seed when present
}



def _lookup_k(policy: Dict[Tuple[str, str], int], res: str, moiety: str, default_k: int) -> int:
    return policy.get((res, moiety), policy.get(("ANY", moiety), default_k))



def _euclid2(a: Tuple[float,float,float], b: Tuple[float,float,float]) -> float:
    dx = a[0]-b[0]; dy = a[1]-b[1]; dz = a[2]-b[2]
    return dx*dx + dy*dy + dz*dz




def _nearest_idx_and_dist(pt: np.ndarray, cloud: np.ndarray) -> Tuple[int, float]:
    """Return (index, distance) of nearest point in cloud to pt."""
    if cloud.size == 0:
        return -1, float("inf")
    diffs = cloud - pt
    d2 = np.einsum("ij,ij->i", diffs, diffs)
    i = int(np.argmin(d2))
    return i, float(np.sqrt(d2[i]))





def _write_coord_links_csv(
    path: str,
    cofactor_atoms: List[AtomDict],
    pcs_seed_atoms: List[AtomDict],
    scs_seed_atoms: List[AtomDict],
) -> int:
    """
    Write link rows:
      - cofactor->pcs : for each PCS seed, link to its nearest cofactor atom
      - pcs->scs      : for each SCS seed,  link to its nearest *PCS seed* atom
    Returns number of rows written.
    """
    rows = []

    # 1) Cofactor -> PCS
    for pcs in pcs_seed_atoms:
        src, dist = _nearest_source_atom(cofactor_atoms, pcs)
        if src is None:
            continue
        rows.append({
            "link_type": "cofactor->pcs",
            "src_resname": src["residue"], "src_resnum": src["residue_number"], "src_chain": src.get("chain",""),
            "src_atom": src["name"], "src_moiety": src.get("moiety", ""),  # moiety optional if you’ve added it
            "dst_resname": pcs["residue"], "dst_resnum": pcs["residue_number"], "dst_chain": pcs.get("chain",""),
            "dst_atom": pcs["name"], "dst_moiety": pcs.get("moiety", ""),
            "distance_A": f"{dist:.3f}",
        })

    # 2) PCS -> SCS (nearest PCS seed)
    for scs in scs_seed_atoms:
        src, dist = _nearest_source_atom(pcs_seed_atoms, scs)
        if src is None:
            continue
        rows.append({
            "link_type": "pcs->scs",
            "src_resname": src["residue"], "src_resnum": src["residue_number"], "src_chain": src.get("chain",""),
            "src_atom": src["name"], "src_moiety": src.get("moiety", ""),
            "dst_resname": scs["residue"], "dst_resnum": scs["residue_number"], "dst_chain": scs.get("chain",""),
            "dst_atom": scs["name"], "dst_moiety": scs.get("moiety", ""),
            "distance_A": f"{dist:.3f}",
        })

    if not rows:
        # still create an empty file with header to make behavior predictable
        header = ["link_type",
                  "src_resname","src_resnum","src_chain","src_atom","src_moiety",
                  "dst_resname","dst_resnum","dst_chain","dst_atom","dst_moiety",
                  "distance_A"]
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=header)
            w.writeheader()
        print(f"[INFO] Coord links: wrote 0 links to '{path}'")
        return 0

    header = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)

    print(f"[INFO] Coord links: wrote {len(rows)} links to '{path}'")
    return len(rows)





# def _write_coord_links_csv(
#     outfile: str,
#     cofactor_atoms: List[AtomDict],
#     pcs_seeds: List[AtomDict],
#     scs_seeds: List[AtomDict],
# ):
#     """
#     Write a tidy edge list:
#       link_type ∈ {"cofactor->pcs","pcs->scs"}
#       src_* (resname,resnum,chain,atom,moiety), dst_* (...), distance_A
#     """
#     cof_coords = np.array([a["coordinates"] for a in cofactor_atoms]) if cofactor_atoms else np.empty((0,3))
#     pcs_coords = np.array([a["coordinates"] for a in pcs_seeds]) if pcs_seeds else np.empty((0,3))

#     rows = []
#     # Links: cofactor -> PCS (nearest cofactor to each PCS seed)
#     for p in pcs_seeds:
#         idx, d = _nearest_idx_and_dist(np.array(p["coordinates"]), cof_coords)
#         if idx >= 0:
#             c = cofactor_atoms[idx]
#             rows.append({
#                 "link_type": "cofactor->pcs",
#                 "src_resname": c["residue"], "src_resnum": c["residue_number"], "src_chain": c["chain"],
#                 "src_atom": c["name"], "src_moiety": _get_moiety_label(str(c["residue"]), str(c["name"])),
#                 "dst_resname": p["residue"], "dst_resnum": p["residue_number"], "dst_chain": p["chain"],
#                 "dst_atom": p["name"], "dst_moiety": _get_moiety_label(str(p["residue"]), str(p["name"])),
#                 "distance_A": f"{d:.3f}",
#             })

#     # Links: PCS -> SCS (nearest PCS to each SCS seed)
#     for s in scs_seeds:
#         idx, d = _nearest_idx_and_dist(np.array(s["coordinates"]), pcs_coords)
#         if idx >= 0:
#             p = pcs_seeds[idx]
#             rows.append({
#                 "link_type": "pcs->scs",
#                 "src_resname": p["residue"], "src_resnum": p["residue_number"], "src_chain": p["chain"],
#                 "src_atom": p["name"], "src_moiety": _get_moiety_label(str(p["residue"]), str(p["name"])),
#                 "dst_resname": s["residue"], "dst_resnum": s["residue_number"], "dst_chain": s["chain"],
#                 "dst_atom": s["name"], "dst_moiety": _get_moiety_label(str(s["residue"]), str(s["name"])),
#                 "distance_A": f"{d:.3f}",
#             })

#     fieldnames = [
#         "link_type",
#         "src_resname","src_resnum","src_chain","src_atom","src_moiety",
#         "dst_resname","dst_resnum","dst_chain","dst_atom","dst_moiety",
#         "distance_A",
#     ]
#     with open(outfile, "w", newline="") as f:
#         w = csv.DictWriter(f, fieldnames=fieldnames)
#         w.writeheader()
#         for r in rows:
#             w.writerow(r)
#     print(f"[INFO] Wrote {len(rows)} links → '{outfile}'")





# def calculate_centroid(residue_atoms):
#     coords = np.array([atom['coordinates'] for atom in residue_atoms])
#     return coords.mean(axis=0)

# def get_best_fit_plane(atom_coords):
#     centroid = np.mean(atom_coords, axis=0)
#     centered_coords = atom_coords - centroid
#     _, _, vh = np.linalg.svd(centered_coords)
#     normal = vh[-1]
#     return normal, centroid

# def extract_residue_centroids(atoms_data):
#     residues = {}
#     for atom in atoms_data:
#         residue_number = atom['residue_number']
#         if residue_number not in residues:
#             residues[residue_number] = []
#         residues[residue_number].append(atom)

#     centroids = {
#         residue_number: calculate_centroid(residue_atoms)
#         for residue_number, residue_atoms in residues.items()
#     }
#     return centroids

# def calculate_residue_distances(query_centroids, template_centroids):
#     distances = {}
#     for template_residue, template_coord in template_centroids.items():
#         query_distances = [
#             np.linalg.norm(template_coord - query_coord)
#             for query_coord in query_centroids.values()
#         ]
#         distances[template_residue] = min(query_distances)
#     return distances

# def process_pdb_file(file_path, template_centroids):
#     from structure_io import unpack_pdb_file

#     structure, atoms_data = unpack_pdb_file(file_path)
#     query_centroids = extract_residue_centroids(atoms_data)
#     distances = calculate_residue_distances(query_centroids, template_centroids)

#     residue_names = [
#         f"{residue_number}" for residue_number in query_centroids.keys()
#         if residue_number in distances
#     ]
#     return distances, residue_names

# def build_matrices(file_paths, template_centroids):
#     distance_matrix = []
#     residue_matrix = []

#     for file_path in file_paths:
#         file_name = os.path.basename(file_path)
#         distances, residue_names = process_pdb_file(file_path, template_centroids)

#         distance_row = [file_name] + [distances.get(res, np.nan) for res in template_centroids.keys()]
#         residue_row = [file_name] + residue_names

#         distance_matrix.append(distance_row)
#         residue_matrix.append(residue_row)

#     return distance_matrix, residue_matrix

# def save_matrices_to_csv(distance_matrix, residue_matrix, output_dir):
#     distance_df = pd.DataFrame(distance_matrix)
#     residue_df = pd.DataFrame(residue_matrix)

#     distance_df.to_csv(os.path.join(output_dir, "distance_matrix.csv"), index=False)
#     residue_df.to_csv(os.path.join(output_dir, "residue_matrix.csv"), index=False)



# # modules/structure_processing.py
# from collections import defaultdict
# import numpy as np




def _nearest_source_atom(source_atoms: List[AtomDict], target_atom: AtomDict) -> Tuple[Optional[AtomDict], float]:
    """Return (nearest_source_atom, distance_A) for the target_atom, searching in source_atoms."""
    if not source_atoms:
        return None, math.inf
    t = target_atom["coordinates"]
    best = None
    best_d2 = math.inf
    for s in source_atoms:
        d2 = _euclid2(s["coordinates"], t)
        if d2 < best_d2:
            best_d2 = d2
            best = s
    return best, math.sqrt(best_d2)







def make_residue_centroid_sphere(template_coord_residues):
    """
    Build a residue 'sphere' by grouping atom dicts by residue number and
    computing the centroid for each residue.

    Parameters
    ----------
    template_coord_residues : list of dict
        Each dict must at least contain:
            - "coordinates": np.ndarray-like (x, y, z)
            - "residue": residue name (e.g., "HEM", "ALA")
            - "residue_number": int

    Returns
    -------
    list of dict
        One entry per residue with:
            - "residue_number" : int
            - "residue_name"   : str
            - "coordinates"    : np.ndarray (centroid)
            - "atom_coords"    : list of np.ndarray (all atom coords in that residue)
    """
    residue_groups = defaultdict(list)

    # Group atoms by residue number (chain is ignored intentionally to
    # match previous behavior; add chain if your use-case requires it)
    for atom in template_coord_residues:
        resnum = atom.get("residue_number")
        resname = atom.get("residue")
        coords = atom.get("coordinates")
        if resnum is None or resname is None or coords is None:
            # Skip malformed entries silently
            continue
        residue_groups[resnum].append(atom)

    residue_sphere = []
    for resnum, atoms in residue_groups.items():
        if not atoms:
            continue
        # Collect coordinates
        coords = np.array([np.asarray(a["coordinates"], dtype=float) for a in atoms])
        centroid = coords.mean(axis=0)
        resname = atoms[0].get("residue", "UNK")
        residue_sphere.append({
            "residue_number": resnum,
            "residue_name": resname,
            "coordinates": centroid,
            "atom_coords": [np.asarray(a["coordinates"], dtype=float) for a in atoms],
        })

    return residue_sphere





# # modules/structure_processing.py
# import numpy as np





def extract_query_box(template_atoms, query_atoms, distance_cutoff):
    """
    From the template atoms (list of dicts), keep only CA/CB, compute the max |coord|
    among those atoms, inflate by distance_cutoff to form a cube centered at the origin,
    then return all query atoms inside that cube (also CA/CB only).

    Returns
    -------
    query_box_filtered : list[dict]
        Query atoms (only CA/CB) within the cube.
    cube_limit : float
        Half-length of the cube in each axis.
    """
    # Filter template atoms for CA and CB
    template_filtered = [a for a in template_atoms if a.get('name') in ('CA', 'CB')]
    if not template_filtered:
        # Fallback: if the template set has no CA/CB, just use all coords
        template_filtered = template_atoms

    # Max absolute coordinate across all dims for the selected template atoms
    all_abs = []
    for a in template_filtered:
        c = np.asarray(a['coordinates'], dtype=float)
        all_abs.extend(np.abs(c))
    cube_limit = (max(all_abs) if all_abs else 0.0) + float(distance_cutoff)

    # Inside cube filter
    def in_cube(a):
        c = np.asarray(a['coordinates'], dtype=float)
        return np.all(np.abs(c) <= cube_limit)

    query_box = [a for a in query_atoms if in_cube(a)]
    # Keep only CA/CB on the query side
    query_box_filtered = [a for a in query_box if a.get('name') in ('CA', 'CB')]

    return query_box_filtered, cube_limit






# from typing import Dict, List, Tuple, Iterable, Set, Optional
# from collections import defaultdict
# import numpy as np

# # If you centralize these elsewhere, you can import instead:
# try:
#     from modules.moieties import chemical_moieties  # your big mapping
# except Exception:
#     chemical_moieties: Dict[Tuple[str, str], str] = {}

# AtomDict = Dict[str, object]  # name,residue,residue_number,chain,element,coordinates

# # -----------------------------
# # Utilities / pretty logging
# # -----------------------------
# def _bioatom_to_dict(residue, atom, chain_id: str) -> AtomDict:
#     return {
#         "name": atom.get_name(),
#         "residue": residue.get_resname(),
#         "residue_number": residue.get_id()[1],
#         "chain": chain_id,
#         "element": getattr(atom, "element", ""),
#         "coordinates": np.array(atom.coord, dtype=float),
#     }

# def _residue_atoms_as_dicts(residue) -> List[AtomDict]:
#     chain_id = residue.get_full_id()[2]
#     return [_bioatom_to_dict(residue, at, chain_id) for at in residue]

# def _group_atoms_by_residue(atoms: Iterable[AtomDict]) -> Dict[Tuple[str, int, str], List[str]]:
#     grouped: Dict[Tuple[str, int, str], List[str]] = defaultdict(list)
#     for a in atoms:
#         key = (str(a["residue"]), int(a["residue_number"]), str(a["chain"]))
#         grouped[key].append(str(a["name"]))
#     for k in list(grouped.keys()):
#         grouped[k] = sorted(set(grouped[k]))
#     return dict(sorted(grouped.items(), key=lambda x: (x[0][2], x[0][0], x[0][1])))

# def _print_atom_set(label: str, atoms: Iterable[AtomDict]) -> None:
#     atoms = list(atoms)
#     print(f"[INFO] {label}: {len(atoms)} atoms")
#     grouped = _group_atoms_by_residue(atoms)
#     for (resname, resnum, chain), names in grouped.items():
#         print(f"[INFO]   {label} ▸ {resname} {resnum} {chain}: {', '.join(names)}")

# def _dedup_atoms(atoms: Iterable[AtomDict]) -> List[AtomDict]:
#     seen: Set[Tuple[str, int, str, str]] = set()
#     out: List[AtomDict] = []
#     for a in atoms:
#         key = (str(a["residue"]), int(a["residue_number"]), str(a["chain"]), str(a["name"]))
#         if key not in seen:
#             seen.add(key)
#             out.append(a)
#     return out

# def _dist(a: np.ndarray, b: np.ndarray) -> float:
#     return float(np.linalg.norm(a - b))

# _BACKBONE_NAMES = {"N", "H", "CA", "HA", "C", "O", "OXT"}

# def _expand_by_moiety_and_backbone_with_logging(
#     seed_atoms: List[AtomDict],
#     structure,
#     exclude_moieties: Optional[List[str]] = None,
# ) -> List[AtomDict]:
#     """
#     For each residue present in seed_atoms:
#       • Identify moiety labels on seed atoms (from chemical_moieties).
#       • Add all atoms in that residue sharing those moieties (except those in exclude_moieties).
#       • Add backbone atoms (N,H,CA,HA,C,O,OXT).
#       • Print [INFO] lines showing exactly what was added.
#     """
#     if exclude_moieties is None:
#         exclude_moieties = []
#     exclude_set = set(m.lower() for m in exclude_moieties)

#     seeds_by_res: Dict[Tuple[str, int, str], Dict[str, Set[str]]] = defaultdict(lambda: {"moieties": set(), "seed_names": set()})
#     for a in seed_atoms:
#         rk = (str(a["residue"]), int(a["residue_number"]), str(a["chain"]))
#         seeds_by_res[rk]["seed_names"].add(str(a["name"]))
#         moiety = chemical_moieties.get((str(a["residue"]), str(a["name"])))
#         if moiety and moiety.lower() not in exclude_set:
#             seeds_by_res[rk]["moieties"].add(moiety)

#     expanded: Dict[Tuple[str, int, str, str], AtomDict] = {}
#     for a in seed_atoms:
#         k = (str(a["residue"]), int(a["residue_number"]), str(a["chain"]), str(a["name"]))
#         expanded[k] = a

#     for model in structure:
#         for chain in model:
#             chain_id = chain.id
#             for residue in chain:
#                 resname = residue.get_resname()
#                 resnum = residue.get_id()[1]
#                 rk = (resname, resnum, chain_id)
#                 if rk not in seeds_by_res:
#                     continue

#                 all_atoms = _residue_atoms_as_dicts(residue)
#                 moieties = seeds_by_res[rk]["moieties"]
#                 seed_names = seeds_by_res[rk]["seed_names"]

#                 # (1) Moiety expansion
#                 added_moiety: List[str] = []
#                 if moieties:
#                     target = {m for m in moieties if m.lower() not in exclude_set}
#                     for atom in all_atoms:
#                         m = chemical_moieties.get((resname, str(atom["name"])))
#                         if m and (m in target) and (m.lower() not in exclude_set):
#                             key = (resname, resnum, chain_id, str(atom["name"]))
#                             if key not in expanded:
#                                 expanded[key] = atom
#                                 if atom["name"] not in seed_names:
#                                     added_moiety.append(str(atom["name"]))
#                     moiety_txt = ", ".join(sorted(target)) if target else "(all excluded)"
#                     if added_moiety:
#                         print(f"[INFO] Moiety expansion: {resname} {resnum} {chain_id} | moieties=[{moiety_txt}] → added: {', '.join(sorted(set(added_moiety)))}")
#                     else:
#                         print(f"[INFO] Moiety expansion: {resname} {resnum} {chain_id} | moieties=[{moiety_txt}] → added: (none)")
#                 else:
#                     print(f"[INFO] Moiety expansion: {resname} {resnum} {chain_id} → no moiety labels on seed atoms; no moiety additions")

#                 # (2) Backbone expansion
#                 added_backbone: List[str] = []
#                 for atom in all_atoms:
#                     if str(atom["name"]) in _BACKBONE_NAMES:
#                         key = (resname, resnum, chain_id, str(atom["name"]))
#                         if key not in expanded:
#                             expanded[key] = atom
#                             if atom["name"] not in seed_names:
#                                 added_backbone.append(str(atom["name"]))
#                 if added_backbone:
#                     print(f"[INFO] Backbone expansion: {resname} {resnum} {chain_id} → added backbone: {', '.join(sorted(set(added_backbone)))}")
#                 else:
#                     print(f"[INFO] Backbone expansion: {resname} {resnum} {chain_id} → added backbone: (none)")

#     return list(expanded.values())

# # ---------------------------------------------------------
# # Main API (keeps your exact call signature)
# # ---------------------------------------------------------
# def identify_coordination_network(
#     structure,
#     cofactor_resname: List[str],
#     distance_cutoff: float,
#     expand_residues: bool,
#     combinatorial_mode: bool,
#     combinatorial_cofactor_cutoff: Optional[float] = None,
#     cofactor_resname2: Optional[List[str]] = None,
#     exclude_moieties: Optional[List[str]] = None,
# ):
#     """
#     Returns: (cofactor_sphere, pcs_atoms, scs_atoms)

#     • Identifies cofactor atoms by residue name(s).
#     • If combinatorial_mode=True and cofactor_resname2 is provided:
#         - Adds atoms from residues in cofactor_resname2 that lie within
#           combinatorial_cofactor_cutoff of the primary cofactor atoms.
#     • Builds PCS seeds as atoms within distance_cutoff of any cofactor atom;
#       SCS seeds as atoms within 2*distance_cutoff of any cofactor atom
#       (excluding PCS & cofactor). Adjust as you like.
#     • If expand_residues=True, runs moiety + backbone expansion with [INFO] logs.
#     • exclude_moieties: list of moiety names to ignore during moiety expansion.
#     • Prints [INFO] summaries for identified and final sets.
#     """
#     if exclude_moieties is None:
#         exclude_moieties = []

#     primary_names = {r.upper() for r in (cofactor_resname or [])}
#     secondary_names = {r.upper() for r in (cofactor_resname2 or [])}

#     print(f"[INFO] Identifying coordination network for {sorted(primary_names)}")

#     # Collect all atoms
#     all_atoms: List[AtomDict] = []
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 for atom in residue:
#                     all_atoms.append(_bioatom_to_dict(residue, atom, chain.id))

#     # Primary cofactor atoms
#     cof_primary = [a for a in all_atoms if str(a["residue"]).upper() in primary_names]

#     # Optional combinatorial extension
#     cof_extended: List[AtomDict] = list(cof_primary)
#     if combinatorial_mode and cofactor_resname2 and secondary_names:
#         if not combinatorial_cofactor_cutoff:
#             combinatorial_cofactor_cutoff = max(3.0, distance_cutoff)  # sensible default
#         pcoords = np.array([a["coordinates"] for a in cof_primary]) if cof_primary else np.empty((0, 3))
#         added = 0
#         if pcoords.size:
#             for a in all_atoms:
#                 if str(a["residue"]).upper() not in secondary_names:
#                     continue
#                 ac = np.array(a["coordinates"])
#                 if np.any([_dist(ac, pc) <= combinatorial_cofactor_cutoff for pc in pcoords]):
#                     cof_extended.append(a)
#                     added += 1
#         cof_extended = _dedup_atoms(cof_extended)
#         print(f"[INFO] Combinatorial cofactor extension: added {added} atoms from {sorted(secondary_names)} within {combinatorial_cofactor_cutoff:.2f} Å")
#     else:
#         if combinatorial_mode and not cofactor_resname2:
#             print("[INFO] combinatorial_mode=True but no cofactor_resname2 provided; skipping combinatorial extension.")

#     cofactor_sphere = _dedup_atoms(cof_extended)
#     _print_atom_set("Cofactor (identified)", cofactor_sphere)

#     # Build PCS seeds: within distance_cutoff of any cofactor atom (non-cofactor)
#     ccoords = np.array([a["coordinates"] for a in cofactor_sphere]) if cofactor_sphere else np.empty((0, 3))
#     pcs_seed: List[AtomDict] = []
#     if ccoords.size:
#         for a in all_atoms:
#             if str(a["residue"]).upper() in primary_names or str(a["residue"]).upper() in secondary_names:
#                 continue
#             ac = np.array(a["coordinates"])
#             if np.any([_dist(ac, cc) <= distance_cutoff for cc in ccoords]):
#                 pcs_seed.append(a)
#     pcs_seed = _dedup_atoms(pcs_seed)
#     _print_atom_set("PCS (seed, pre-expansion)", pcs_seed)

#     # Build SCS seeds: within 2*distance_cutoff of any cofactor atom, excluding PCS & cofactor
#     scs_seed: List[AtomDict] = []
#     if ccoords.size:
#         pcs_keys = {(a["residue"], a["residue_number"], a["chain"], a["name"]) for a in pcs_seed}
#         cof_keys = {(a["residue"], a["residue_number"], a["chain"], a["name"]) for a in cofactor_sphere}
#         scs_cut = float(distance_cutoff) * 2.0
#         for a in all_atoms:
#             k = (a["residue"], a["residue_number"], a["chain"], a["name"])
#             if k in pcs_keys or k in cof_keys:
#                 continue
#             ac = np.array(a["coordinates"])
#             if np.any([_dist(ac, cc) <= scs_cut for cc in ccoords]):
#                 scs_seed.append(a)
#     scs_seed = _dedup_atoms(scs_seed)
#     _print_atom_set("SCS (seed, pre-expansion)", scs_seed)

#     # Optional moiety + backbone expansion
#     if expand_residues:
#         pcs_atoms = _expand_by_moiety_and_backbone_with_logging(pcs_seed, structure, exclude_moieties=exclude_moieties)
#         scs_atoms = _expand_by_moiety_and_backbone_with_logging(scs_seed, structure, exclude_moieties=exclude_moieties)
#     else:
#         pcs_atoms, scs_atoms = pcs_seed, scs_seed

#     # Final summaries
#     _print_atom_set("Cofactor (final)", cofactor_sphere)
#     _print_atom_set("PCS (final)", pcs_atoms)
#     _print_atom_set("SCS (final)", scs_atoms)

#     return cofactor_sphere, pcs_atoms, scs_atoms


















# # Add near the top if missing:
# from typing import List, Dict, Union
# import numpy as np



def identify_residues_of_interest(
    structure,
    cofactor_resname: Union[str, List[str]],
    alternative_resname: Union[str, List[str], None] = None,
    combinatorial_mode: bool = False,
    combinatorial_cofactor_cutoff: float = 10.0,
    cofactor_resname2: Union[str, List[str], None] = None,
    residues_of_interest_integers: Union[List[int], None] = None,
) -> (List[Dict], List[Dict]):
    """
    Identify cofactor atoms and atoms in template residues of interest.

    Returns:
      (cofactor_atoms, residues_of_interest_atoms)

      Each list contains dicts with:
      'coordinates', 'name', 'element', 'residue', 'residue_number', 'chain'
    """
    # Normalize input lists
    if isinstance(cofactor_resname, str):
        cofactor_res_list = [cofactor_resname]
    else:
        cofactor_res_list = list(cofactor_resname)

    all_cofactor_resnames = cofactor_res_list.copy()
    if alternative_resname:
        if isinstance(alternative_resname, str):
            all_cofactor_resnames.append(alternative_resname)
        else:
            all_cofactor_resnames.extend(alternative_resname)

    # Gather cofactor atoms (Bio.PDB Atom objects) using our helper
    # NOTE: this helper should be defined in modules/structure_utils and imported at top of this file


    bio_cofactor_atoms = []
    for rn in all_cofactor_resnames:
        bio_cofactor_atoms.extend(find_cofactor_atoms(structure, rn))

    # Combinatorial filtering (optional)
    if combinatorial_mode:
        if cofactor_resname2 is None:
            raise ValueError("In combinatorial_mode, cofactor_resname2 must be provided.")
        if isinstance(cofactor_resname2, str):
            c2_list = [cofactor_resname2]
        else:
            c2_list = list(cofactor_resname2)

        group1 = [a for a in bio_cofactor_atoms if a.get_parent().get_resname().upper() in {r.upper() for r in cofactor_res_list}]
        group2 = [a for a in bio_cofactor_atoms if a.get_parent().get_resname().upper() in {r.upper() for r in c2_list}]

        # Keep only atoms that have a neighbor in the other group within cutoff
        filtered_g1 = [
            a1 for a1 in group1
            if any(np.linalg.norm(a1.coord - a2.coord) <= combinatorial_cofactor_cutoff for a2 in group2)
        ]
        filtered_g2 = [
            a2 for a2 in group2
            if any(np.linalg.norm(a2.coord - a1.coord) <= combinatorial_cofactor_cutoff for a1 in group1)
        ]
        bio_cofactor_atoms = filtered_g1 + filtered_g2

    # Convert Biopython atoms → dicts (uniform return type)
    def atom_to_dict(a):
        res = a.get_parent()
        chain = res.get_parent()
        return {
            "coordinates": a.coord,
            "name": a.get_name(),
            "element": a.element,
            "residue": res.get_resname(),
            "residue_number": res.get_id()[1],
            "chain": chain.get_id(),
        }

    cofactor_atoms = [atom_to_dict(a) for a in bio_cofactor_atoms]

    # Build residues_of_interest list (template-relative ints) → +1 → strings
    residues_of_interest_atoms: List[Dict] = []
    roi_strs: List[str] = []
    if residues_of_interest_integers is not None:
        roi_strs = [str(x + 1) for x in residues_of_interest_integers]

    # Collect all atoms from those residues
    if roi_strs:
        for model in structure:
            for chain in model:
                for residue in chain:
                    resnum_str = str(residue.get_id()[1])
                    if resnum_str in roi_strs:
                        for atom in residue:
                            residues_of_interest_atoms.append({
                                "coordinates": atom.get_coord(),
                                "name": atom.get_name(),
                                "element": atom.element,
                                "residue": residue.get_resname(),
                                "residue_number": residue.get_id()[1],
                                "chain": chain.get_id(),
                            })

    return cofactor_atoms, residues_of_interest_atoms










def generate_coordination_csv_with_moieties(
    cofactor_sphere: List[Dict],
    pcs_residues: List[Dict],
    scs_residues: List[Dict],
    pdb_file_name: str,
    structure,
    bond_lookup: Dict,
    output_dir: Optional[str] = None,
) -> None:
    """
    Generate a CSV summarizing the PCS and SCS residue names, chemical moieties, and interactors.

    Arguments
    ---------
    cofactor_sphere : list[dict]
        Atoms belonging to the cofactor (as dicts with coordinates, name, element, residue, residue_number, chain).
    pcs_residues : list[dict]
        Primary coordination sphere atoms (same dict structure).
    scs_residues : list[dict]
        Secondary coordination sphere atoms (same dict structure).
    pdb_file_name : str
        Used as the CSV file prefix.
    structure : Bio.PDB structure
        Biopython structure object.
    bond_lookup : dict
        Present for signature parity (not used directly here).

    Output
    ------
    Writes `<pdb_file_name>_Coord_Breakdown.csv` to `output_dir` (or CWD if not provided).
    """
    out_dir = output_dir or "."
    os.makedirs(out_dir, exist_ok=True)
    output_file_name = os.path.join(out_dir, f"{pdb_file_name}_Coord_Breakdown.csv")
    rows = []

    # --- Helpers -------------------------------------------------------------

    def get_residue_object(atom_dict):
        """Return the Biopython Residue object matching an atom-dict origin."""
        target_resname = atom_dict["residue"]
        target_resnum  = atom_dict["residue_number"]
        target_chain   = atom_dict["chain"]
        for res in structure.get_residues():
            if (res.get_resname() == target_resname
                and res.get_id()[1] == target_resnum
                and res.get_full_id()[2] == target_chain):
                return res
        raise ValueError(
            f"Residue not found for {target_resname} {target_resnum} chain {target_chain}"
        )

    def get_interactors(atom_dict, atoms_list):
        """
        Find atom names in `atoms_list` within 3.0 Å of `atom_dict`.
        """
        interactors = []
        a_coord = np.asarray(atom_dict["coordinates"], dtype=float)
        for other in atoms_list:
            # don't compare to self-object by identity
            if other is atom_dict:
                continue
            b_coord = np.asarray(other["coordinates"], dtype=float)
            if np.linalg.norm(a_coord - b_coord) <= 3.0:
                interactors.append(other["name"])
        return interactors

    def extract_residue_data(atoms: List[Dict], sphere_label: str):
        """
        Aggregate per-residue moieties + interactors and append to rows.
        """
        per_res = {}
        for a in atoms:
            key = (a["residue"], a["residue_number"], a["chain"])
            if key not in per_res:
                per_res[key] = {
                    "residue_name": a["residue"],
                    "residue_number": a["residue_number"],
                    "chain": a["chain"],
                    "moiety": set(),
                    "interactors": set(),
                }

            # moiety expansion for this atom
            try:
                residue_obj = get_residue_object(a)
                bp_atom = next(x for x in residue_obj if x.get_name() == a["name"])
                moiety_atoms = get_moiety_atoms(bp_atom, residue_obj)
                per_res[key]["moiety"].update(x.get_name() for x in moiety_atoms)
            except Exception as e:
                # Keep robust; just skip moiety expansion if something is off
                # print(f"[WARN] Moiety expansion failed for {a}: {e}")
                pass

            # interactors among PCS+SCS sets (excluding cofactor set)
            interactors = get_interactors(a, pcs_residues + scs_residues)
            per_res[key]["interactors"].update(interactors)

        # materialize rows
        for (rname, rnum, ch), payload in per_res.items():
            rows.append({
                "Residue Category": sphere_label,
                "Residue Name": payload["residue_name"],
                "Residue Number": payload["residue_number"],
                "Chain": payload["chain"],
                "Chemical Moiety": ", ".join(sorted(payload["moiety"])) if payload["moiety"] else "",
                "Interactors": ", ".join(sorted(payload["interactors"])) if payload["interactors"] else "",
            })

    # --- Build rows for each sphere -----------------------------------------
    extract_residue_data(cofactor_sphere, "Cofactor")
    extract_residue_data(pcs_residues, "PCS")
    extract_residue_data(scs_residues, "SCS")

    # --- Write CSV -----------------------------------------------------------
    with open(output_file_name, mode="w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Residue Category", "Residue Name", "Residue Number",
            "Chain", "Chemical Moiety", "Interactors"
        ])
        writer.writeheader()
        writer.writerows(rows)

        # Append atom-level coordinates for all residues referenced above.
        plain_writer = csv.writer(f)
        plain_writer.writerow([])  # blank line between sections
        plain_writer.writerow([
            "Residue Category", "Residue Name", "Residue Number",
            "Chain", "Atom", "x", "y", "z"
        ])

        def _rkey(resname, resnum, chain):
            return (str(resname), str(resnum), str(chain))

        res_category = {}
        for a in scs_residues:
            res_category[_rkey(a["residue"], a["residue_number"], a["chain"])] = "SCS"
        for a in pcs_residues:
            res_category[_rkey(a["residue"], a["residue_number"], a["chain"])] = "PCS"
        for a in cofactor_sphere:
            res_category[_rkey(a["residue"], a["residue_number"], a["chain"])] = "Cofactor"

        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                for residue in chain:
                    resname = residue.get_resname()
                    resnum = residue.get_id()[1]
                    category = res_category.get(_rkey(resname, resnum, chain_id))
                    if category is None:
                        continue
                    for atom in residue:
                        x, y, z = map(float, atom.get_coord())
                        plain_writer.writerow([
                            category,
                            resname,
                            resnum,
                            chain_id,
                            atom.get_name(),
                            x, y, z
                        ])

    print(f"[INFO] Coordination breakdown written to: {output_file_name}")





def get_residue_bonds(residue_atoms: List[Dict], bond_lookup: Dict[str, List[Tuple[str, str]]]):
    """
    Given a list of atom dicts for *one residue* and a bond lookup table,
    return a list of bonds as tuples of 3D coordinates: [((x1,y1,z1), (x2,y2,z2)), ...].

    residue_atoms: [
        {'name': 'CA', 'element': 'C', 'coordinates': np.array([..,..,..]),
         'residue': 'GLY', 'residue_number': 42, 'chain': 'A'},
        ...
    ]
    bond_lookup: {'GLY': [('N','CA'), ('CA','C'), ('C','O')], ...}
    """
    # basic validation
    if not residue_atoms or not isinstance(residue_atoms, list) or not isinstance(residue_atoms[0], dict):
        return []

    residue_name = residue_atoms[0].get("residue")
    if not residue_name:
        return []

    # skip solvent/ions quickly
    if residue_name in ("HOH", "WAT", "NA"):
        return []

    rules = bond_lookup.get(residue_name)
    if not rules:
        # silently skip residues not in the lookup (common for ligands/hets)
        return []

    # map atom name -> coord
    atom_coords = {a.get("name"): a.get("coordinates") for a in residue_atoms if a.get("name") and a.get("coordinates") is not None}

    bonds = []
    for a, b in rules:
        ca = atom_coords.get(a)
        cb = atom_coords.get(b)
        if ca is not None and cb is not None:
            bonds.append((tuple(ca), tuple(cb)))
    return bonds





#     # modules/structure_processing.py  (add near your other bond helpers)

# from typing import Dict, List, Tuple



def generate_residue_bonds(
    residues: List[Dict],
    bond_lookup: Dict[str, List[Tuple[str, str]]]
) -> List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    """
    Build bonds for a set of atoms (from multiple residues) using a residue-specific lookup.

    residues: list of atom dicts with keys:
      - 'name', 'element', 'coordinates' (np.array-like), 'residue' (resname),
        'residue_number', 'chain'
    bond_lookup: {'RESNAME': [('N','CA'), ('CA','C'), ...], ...}

    Returns: list of bonds as pairs of 3D coords: [((x1,y1,z1),(x2,y2,z2)), ...]
    """
    bonds: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = []

    # group atoms by residue identity (resname, number, chain)
    groups: Dict[Tuple[str, int, str], List[Dict]] = {}
    for atom in residues:
        key = (str(atom.get('residue','')).upper(),
               int(atom.get('residue_number')),
               str(atom.get('chain')))
        groups.setdefault(key, []).append(atom)

    # within each residue, apply the bond lookup
    for (resname, _resnum, _chain), atoms in groups.items():
        if resname not in bond_lookup:
            continue
        atom_pos = {a['name'].strip().upper(): a['coordinates'] for a in atoms}
        for a1, a2 in bond_lookup[resname]:
            a1u, a2u = a1.strip().upper(), a2.strip().upper()
            if a1u in atom_pos and a2u in atom_pos:
                p1 = tuple(map(float, atom_pos[a1u]))
                p2 = tuple(map(float, atom_pos[a2u]))
                bonds.append((p1, p2))
    return bonds






def generate_all_bonds(
    structure,
    bond_lookup: Dict[str, List[Tuple[str, str]]]
) -> List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    """
    Build bonds for **every residue** in a Biopython structure using bond_lookup.

    structure: Biopython Structure
    bond_lookup: {'RESNAME': [('N','CA'), ('CA','C'), ...], ...}
    """
    all_bonds: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = []

    for residue in structure.get_residues():
        resname = residue.get_resname()
        if not resname:
            continue
        resname_u = resname.strip().upper()
        if resname_u not in bond_lookup:
            continue

        # build per-atom dicts in the format expected by get_residue_bonds
        atoms = [{
            'name': atom.get_name().strip().upper(),
            'element': atom.element,
            'coordinates': atom.coord,
            'residue': resname_u,
            'residue_number': residue.get_id()[1],
            'chain': residue.get_full_id()[2],
        } for atom in residue]

        if not atoms:
            continue

        all_bonds.extend(generate_residue_bonds(atoms, bond_lookup))

    return all_bonds






# # ensure these are exported
# try:
#     __all__
# except NameError:
#     __all__ = []

# __all__ += ["generate_residue_bonds", "generate_all_bonds"]







# If you centralize these, import instead:
try:
    from modules.moieties import chemical_moieties  # your big mapping
except Exception:
    chemical_moieties: Dict[Tuple[str, str], str] = {}

# -----------------------------
# Constants / small preferences
# -----------------------------
_BACKBONE_NAMES = {"N", "H", "CA", "HA", "C", "O", "OXT"}

# Optional tie-breaker per moiety label (first hits are preferred)
_PREFERRED_BY_MOIETY: Dict[str, List[str]] = {
    "thiol": ["SG", "CB"],
    "thioether": ["SD", "CG", "CE"],
    "imidazole": ["ND1", "NE2", "CG", "CE1", "CD2"],
    "COO": ["OE1", "OE2", "OD1", "OD2", "CD", "CG"],
    "amine": ["NZ", "CE", "CD"],
    "phenol": ["OH", "CZ", "CE1", "CE2"],
    "hydroxyl": ["OG", "OG1"],
    "backbone": ["O", "N", "C", "CA", "OXT", "H", "HA"],
    # add others as you like…
}

# -----------------------------
# Utilities / pretty logging
# -----------------------------
def _bioatom_to_dict(residue, atom, chain_id: str) -> AtomDict:
    return {
        "name": atom.get_name(),
        "residue": residue.get_resname(),
        "residue_number": residue.get_id()[1],
        "chain": chain_id,
        "element": getattr(atom, "element", ""),
        "coordinates": np.array(atom.coord, dtype=float),
    }

def _residue_atoms_as_dicts(residue) -> List[AtomDict]:
    chain_id = residue.get_full_id()[2]
    return [_bioatom_to_dict(residue, at, chain_id) for at in residue]

def _group_atoms_by_residue(atoms: Iterable[AtomDict]) -> Dict[Tuple[str, int, str], List[str]]:
    grouped: Dict[Tuple[str, int, str], List[str]] = defaultdict(list)
    for a in atoms:
        key = (str(a["residue"]), int(a["residue_number"]), str(a["chain"]))
        grouped[key].append(str(a["name"]))
    for k in list(grouped.keys()):
        grouped[k] = sorted(set(grouped[k]))
    return dict(sorted(grouped.items(), key=lambda x: (x[0][2], x[0][0], x[0][1])))

def _print_atom_set(label: str, atoms: Iterable[AtomDict]) -> None:
    atoms = list(atoms)
    print(f"[INFO] {label}: {len(atoms)} atoms")
    grouped = _group_atoms_by_residue(atoms)
    for (resname, resnum, chain), names in grouped.items():
        print(f"[INFO]   {label} ▸ {resname} {resnum} {chain}: {', '.join(names)}")

def _dedup_atoms(atoms: Iterable[AtomDict]) -> List[AtomDict]:
    seen: Set[Tuple[str, int, str, str]] = set()
    out: List[AtomDict] = []
    for a in atoms:
        key = (str(a["residue"]), int(a["residue_number"]), str(a["chain"]), str(a["name"]))
        if key not in seen:
            seen.add(key)
            out.append(a)
    return out

def _min_dist_to_cloud(pt: np.ndarray, cloud: np.ndarray) -> float:
    if cloud.size == 0:
        return float("inf")
    diffs = cloud - pt
    d2 = np.einsum('ij,ij->i', diffs, diffs)
    return float(np.sqrt(d2.min()))




# def _get_moiety_label(resname: str, atomname: str) -> str:
#     m = chemical_moieties.get((resname, atomname))
#     if m:
#         return m
#     if atomname in _BACKBONE_NAMES:
#         return "backbone"
#     return "unknown_moiety"




# -- Moiety lookup
def _get_moiety_label(resname: str, atom_name: str) -> str:
    try:
        from modules.moieties import chemical_moieties as _CHEM_MOIETIES  # type: ignore
    except Exception:
        _CHEM_MOIETIES = {}
    return _CHEM_MOIETIES.get((resname, atom_name), "unknown_moiety")



def _tie_break_prefer(atom_a: str, atom_b: str, moiety: str) -> bool:
    """Return True if atom_a is preferred over atom_b for this moiety."""
    pref = _PREFERRED_BY_MOIETY.get(moiety, [])
    try:
        ia = pref.index(atom_a)
    except ValueError:
        ia = 10_000
    try:
        ib = pref.index(atom_b)
    except ValueError:
        ib = 10_000
    if ia != ib:
        return ia < ib
    # fallback: lexicographic
    return atom_a < atom_b

Candidate = namedtuple("Candidate", ["atom", "distance"])

def _keep_closest(candidates: Dict[Tuple, Candidate], key: Tuple, atom: AtomDict, dist: float, moiety: str) -> None:
    cur = candidates.get(key)
    if cur is None:
        candidates[key] = Candidate(atom, dist)
        return
    eps = 1e-6
    if dist + eps < cur.distance:
        candidates[key] = Candidate(atom, dist)
    elif abs(dist - cur.distance) <= eps:
        # tie-breaker
        a_new = str(atom["name"])
        a_old = str(cur.atom["name"])
        if _tie_break_prefer(a_new, a_old, moiety):
            candidates[key] = Candidate(atom, dist)



# --------------------------------------------
# Expansion: moiety + backbone (with [INFO])
# --------------------------------------------
def _expand_by_moiety_and_backbone_with_logging(
    seed_atoms: List[AtomDict],
    structure,
    exclude_moieties: Optional[List[str]] = None,
) -> List[AtomDict]:
    if exclude_moieties is None:
        exclude_moieties = []
    exclude_set = {m.lower() for m in exclude_moieties}

    # index seeds per residue + gather their moiety set
    seeds_by_res: Dict[Tuple[str, int, str], Dict[str, Set[str]]] = defaultdict(lambda: {"moieties": set(), "seed_names": set()})
    for a in seed_atoms:
        rk = (str(a["residue"]), int(a["residue_number"]), str(a["chain"]))
        seeds_by_res[rk]["seed_names"].add(str(a["name"]))
        m = _get_moiety_label(str(a["residue"]), str(a["name"]))
        if m and m.lower() not in exclude_set:
            seeds_by_res[rk]["moieties"].add(m)

    expanded: Dict[Tuple[str, int, str, str], AtomDict] = {}
    for a in seed_atoms:
        k = (str(a["residue"]), int(a["residue_number"]), str(a["chain"]), str(a["name"]))
        expanded[k] = a

    for model in structure:
        for chain in model:
            chain_id = chain.id
            for residue in chain:
                resname = residue.get_resname()
                resnum = residue.get_id()[1]
                rk = (resname, resnum, chain_id)
                if rk not in seeds_by_res:
                    continue

                all_atoms = _residue_atoms_as_dicts(residue)
                moieties = {m for m in seeds_by_res[rk]["moieties"] if m.lower() not in exclude_set}
                seed_names = seeds_by_res[rk]["seed_names"]

                # (1) Moiety expansion
                added_moiety: List[str] = []
                if moieties:
                    for atom in all_atoms:
                        m = _get_moiety_label(resname, str(atom["name"]))
                        if m in moieties:
                            key = (resname, resnum, chain_id, str(atom["name"]))
                            if key not in expanded:
                                expanded[key] = atom
                                if atom["name"] not in seed_names:
                                    added_moiety.append(str(atom["name"]))
                    mtxt = ", ".join(sorted(moieties)) if moieties else "(all excluded)"
                    if added_moiety:
                        print(f"[INFO] Moiety expansion: {resname} {resnum} {chain_id} | moieties=[{mtxt}] → added: {', '.join(sorted(set(added_moiety)))}")
                    else:
                        print(f"[INFO] Moiety expansion: {resname} {resnum} {chain_id} | moieties=[{mtxt}] → added: (none)")
                else:
                    print(f"[INFO] Moiety expansion: {resname} {resnum} {chain_id} → no moiety labels (or all excluded); no moiety additions")

                # (2) Backbone expansion
                added_backbone: List[str] = []
                for atom in all_atoms:
                    if str(atom["name"]) in _BACKBONE_NAMES:
                        key = (resname, resnum, chain_id, str(atom["name"]))
                        if key not in expanded:
                            expanded[key] = atom
                            if atom["name"] not in seed_names:
                                added_backbone.append(str(atom["name"]))
                if added_backbone:
                    print(f"[INFO] Backbone expansion: {resname} {resnum} {chain_id} → added backbone: {', '.join(sorted(set(added_backbone)))}")
                else:
                    print(f"[INFO] Backbone expansion: {resname} {resnum} {chain_id} → added backbone: (none)")

    return list(expanded.values())



# def _akey(a: AtomDict) -> Tuple[str, int, str, str]:
#     return (str(a["residue"]), int(a["residue_number"]), str(a["chain"]), str(a["name"]))

# def _rkey(a: AtomDict) -> Tuple[str, int, str]:
#     return (str(a["residue"]), int(a["residue_number"]), str(a["chain"]))

# def _mk(a: AtomDict) -> Tuple[str, int, str, str]:
#     """Group by residue and moiety: (resname, resnum, chain, moiety_label)."""
#     return (str(a["residue"]), int(a["residue_number"]), str(a["chain"]),
#             _get_moiety_label(str(a["residue"]), str(a["name"])))

# def _np_coords(atoms: List[AtomDict]) -> np.ndarray:
#     return np.array([a["coordinates"] for a in atoms], dtype=float) if atoms else np.empty((0,3))




# def _pairwise_min_dist(A: np.ndarray, B: np.ndarray) -> np.ndarray:
#     """Return an array of shape (len(A),) with min distance from each row in A to ANY row in B."""
#     if A.size == 0 or B.size == 0:
#         return np.full((len(A),), np.inf)
#     # (i,j,3) diffs
#     diffs = A[:, None, :] - B[None, :, :]
#     d2 = np.einsum("ijk,ijk->ij", diffs, diffs)
#     return np.sqrt(np.min(d2, axis=1))

# def _closest_to_set(targets: np.ndarray, cloud: np.ndarray) -> Tuple[int, float]:
#     """Closest index in 'targets' to ANY point in 'cloud' (return index into targets, distance)."""
#     if targets.size == 0 or cloud.size == 0:
#         return -1, float("inf")
#     diffs = targets[:, None, :] - cloud[None, :, :]
#     d2 = np.einsum("ijk,ijk->ij", diffs, diffs)
#     j = np.argmin(d2, axis=1)
#     i = int(np.argmin(d2[np.arange(len(targets)), j]))
#     return i, float(np.sqrt(d2[i, j[i]]))



# -- Keys and utilities
def _akey(a: AtomDict) -> Tuple[str, int, str, str]:
    return (str(a["residue"]), int(a["residue_number"]), str(a["chain"]), str(a["name"]))

def _rkey(a: AtomDict) -> Tuple[str, int, str]:
    return (str(a["residue"]), int(a["residue_number"]), str(a["chain"]))

def _mk(a: AtomDict) -> Tuple[str, int, str, str]:
    """(resname, resnum, chain, moiety)"""
    res = str(a["residue"])
    return (res, int(a["residue_number"]), str(a["chain"]), _get_moiety_label(res, str(a["name"])))

def _np_coords(atoms: List[AtomDict]) -> np.ndarray:
    return np.array([a["coordinates"] for a in atoms], dtype=float) if atoms else np.empty((0,3))

def _pairwise_min_dist(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    if A.size == 0 or B.size == 0:
        return np.full((len(A),), np.inf)
    diffs = A[:, None, :] - B[None, :, :]
    d2 = np.einsum("ijk,ijk->ij", diffs, diffs)
    return np.sqrt(np.min(d2, axis=1))

def _closest_to_set(targets: np.ndarray, cloud: np.ndarray) -> Tuple[int, float]:
    if targets.size == 0 or cloud.size == 0:
        return -1, float("inf")
    diffs = targets[:, None, :] - cloud[None, :, :]
    d2 = np.einsum("ijk,ijk->ij", diffs, diffs)
    j = np.argmin(d2, axis=1)
    i = int(np.argmin(d2[np.arange(len(targets)), j]))
    return i, float(np.sqrt(d2[i, j[i]]))



# def _filter_seed_candidates(
#     atoms: List[AtomDict],
#     exclude_residue_keys: Set[Tuple[str,int,str]],
#     max_dist_from_set: float,
#     reference_set_coords: np.ndarray,
# ) -> List[AtomDict]:
#     """
#     Filter potential coordinator *seeds*:
#       - element NOT in EXCLUDED_COORD_ELEMENTS
#       - residue NOT in exclude_residue_keys
#       - within max_dist_from_set to the reference set (cofactor for PCS, PCS for SCS)
#     """
#     if reference_set_coords.size == 0:
#         return []

#     kept: List[AtomDict] = []
#     coords = _np_coords(atoms)
#     dmin = _pairwise_min_dist(coords, reference_set_coords)  # distance to nearest reference atom

#     for a, d in zip(atoms, dmin):
#         elem = str(a.get("element","")).upper()
#         if elem in EXCLUDED_COORD_ELEMENTS:
#             continue
#         if _rkey(a) in exclude_residue_keys:
#             continue
#         if d <= max_dist_from_set:
#             kept.append(a)
#     return kept


# -- Seed candidate filtering (excludes carbons, honors distance & residue exclusions)
def _filter_seed_candidates(
    atoms: List[AtomDict],
    exclude_residue_keys: Set[Tuple[str,int,str]],
    max_dist_from_set: float,
    reference_set_coords: np.ndarray,
) -> List[AtomDict]:
    if reference_set_coords.size == 0:
        return []
    kept: List[AtomDict] = []
    coords = _np_coords(atoms)
    dmin = _pairwise_min_dist(coords, reference_set_coords)
    for a, d in zip(atoms, dmin):
        elem = str(a.get("element","")).upper()
        if elem in EXCLUDED_COORD_ELEMENTS:
            continue
        if _rkey(a) in exclude_residue_keys:
            continue
        if d <= max_dist_from_set:
            kept.append(a)
    return kept




# -- Policy-aware picker: choose up to K seeds per (res,moiety) group, nearest first
def _multi_pick_for_group(
    candidates: List[AtomDict],
    reference_set_coords: np.ndarray,
    policy: Dict[Tuple[str, str], int],
    default_k: int = 1,
) -> List[AtomDict]:
    """
    candidates: already filtered list (by element, distance, residue exclusions)
    Group by (resname,resnum,chain,moiety). For each group, pick up to K atoms closest
    to the reference set (cofactor for PCS; PCS for SCS). K is looked up by (resname, moiety)
    from `policy`, defaulting to 1 if not found.
    """
    if not candidates or reference_set_coords.size == 0:
        return []

    grouped: Dict[Tuple[str,int,str,str], List[AtomDict]] = defaultdict(list)
    for a in candidates:
        grouped[_mk(a)].append(a)

    winners: List[AtomDict] = []
    # for (res, num, chain, moi), atoms in grouped.items():
    #     k = policy.get((res, moi), default_k)
    #     if k <= 0:
    #         continue

    for (res, num, chain, moi), atoms in grouped.items():
        k = _lookup_k(policy, res, moi, default_k)
        if k <= 0:
            continue
        

        coords = _np_coords(atoms)
        # sort indices by distance to reference set (ascending)
        diffs = coords[:, None, :] - reference_set_coords[None, :, :]
        d2 = np.einsum("ijk,ijk->ij", diffs, diffs).min(axis=1)
        order = np.argsort(d2)
        pick_idx = order[: min(k, len(order))]
        winners.extend(atoms[i] for i in pick_idx)
    return winners


def _choose_one_per_moiety_group(
    candidates: List[AtomDict],
    reference_set_coords: np.ndarray
) -> List[AtomDict]:
    """
    From candidate atoms, pick exactly ONE atom per (resname,resnum,chain,moiety) group:
    the atom closest to the reference_set (cofactor or PCS).
    """
    if not candidates or reference_set_coords.size == 0:
        return []

    by_group: Dict[Tuple[str,int,str,str], List[AtomDict]] = defaultdict(list)
    for a in candidates:
        by_group[_mk(a)].append(a)

    winners: List[AtomDict] = []
    for g, atoms in by_group.items():
        coords = _np_coords(atoms)
        idx, _ = _closest_to_set(coords, reference_set_coords)
        if idx >= 0:
            winners.append(atoms[idx])
    return winners






# ---------------------------------------------------------
# Main API (keeps your exact call signature)
# ---------------------------------------------------------
# def identify_coordination_network(
#     structure,
#     cofactor_resname: List[str],
#     distance_cutoff: float,
#     expand_residues: bool,
#     combinatorial_mode: bool,
#     combinatorial_cofactor_cutoff: Optional[float] = None,
#     cofactor_resname2: Optional[List[str]] = None,
#     exclude_moieties: Optional[List[str]] = None,
# ):


def identify_coordination_network(
    structure,
    cofactor_resname: List[str],
    distance_cutoff: float,
    expand_residues: bool,
    combinatorial_mode: bool,
    combinatorial_cofactor_cutoff: Optional[float] = None,
    cofactor_resname2: Optional[List[str]] = None,
    exclude_moieties: Optional[List[str]] = None,
    output_dir: Optional[str] = None,
    output_prefix: str = "",
):
    """
    SAME SIGNATURE.
    Changes:
      • PCS seeds: allow up to 2 for ASN/GLN amide, up to 3 for ARG guanidinium (others = 1);
        still excludes element 'C' from seeds.
      • SCS seeds: same multi-policy, chosen against PCS seeds (not 2× cutoff), and excludes cofactor+PCS residues.
    """
    # --- Collect all atoms from structure (same as before) ---
    all_atoms: List[AtomDict] = []
    for model in structure:
        for chain in model:
            chain_id = chain.id
            for residue in chain:
                resname = residue.get_resname()
                resnum  = residue.get_id()[1]
                for atom in residue:
                    all_atoms.append({
                        "name": atom.get_name(),
                        "residue": resname,
                        "residue_number": resnum,
                        "chain": chain_id,
                        "element": getattr(atom, "element", ""),
                        "coordinates": np.array(atom.coord, dtype=float),
                    })

    # --- Cofactor atoms ---
    cof_resnames = set(x.upper() for x in (cofactor_resname or []))
    if cofactor_resname2:
        cof_resnames |= set(x.upper() for x in cofactor_resname2)
    cofactor_sphere = [a for a in all_atoms if str(a["residue"]).upper() in cof_resnames]
    cof_coords = _np_coords(cofactor_sphere)

    print(f"[INFO] Cofactor (identified): {len(cofactor_sphere)} atoms")
    for a in cofactor_sphere:
        print(f"[INFO]   Cofactor ▸ {a['residue']} {a['residue_number']} {a['chain']}: {a['name']}")

    # --- PCS seeds: filter & pick (policy-aware) ---
    pcs_candidates = _filter_seed_candidates(
        atoms=[a for a in all_atoms if str(a["residue"]).upper() not in cof_resnames],
        exclude_residue_keys={_rkey(a) for a in cofactor_sphere},
        max_dist_from_set=distance_cutoff,
        reference_set_coords=cof_coords,
    )
    pcs_seed = _multi_pick_for_group(
        pcs_candidates,
        reference_set_coords=cof_coords,
        policy=MULTI_COORD_POLICY,
        default_k=1,
    )

    print(f"[INFO] PCS seeds (policy-applied): {len(pcs_seed)} atoms")
    pcs_by_res = defaultdict(list)
    for a in pcs_seed: pcs_by_res[_rkey(a)].append(a["name"])
    for (res, num, ch), names in sorted(pcs_by_res.items()):
        print(f"[INFO]   PCS seed ▸ {res} {num} {ch}: {', '.join(sorted(names))}")

    # --- SCS seeds: filter & pick (policy-aware) against PCS; exclude cofactor+PCS residues ---
    pcs_coords = _np_coords(pcs_seed)
    scs_candidates = _filter_seed_candidates(
        atoms=[a for a in all_atoms if str(a["residue"]).upper() not in cof_resnames],
        exclude_residue_keys=({_rkey(a) for a in cofactor_sphere} | {_rkey(a) for a in pcs_seed}),
        max_dist_from_set=distance_cutoff,             # SCS threshold is to PCS seed set
        reference_set_coords=pcs_coords,
    )
    scs_seed = _multi_pick_for_group(
        scs_candidates,
        reference_set_coords=pcs_coords,
        policy=MULTI_COORD_POLICY,
        default_k=1,
    )




    print(f"[INFO] SCS seeds (policy-applied): {len(scs_seed)} atoms")
    scs_by_res = defaultdict(list)
    for a in scs_seed: scs_by_res[_rkey(a)].append(a["name"])
    for (res, num, ch), names in sorted(scs_by_res.items()):
        print(f"[INFO]   SCS seed ▸ {res} {num} {ch}: {', '.join(sorted(names))}")




    # After pcs_seed and scs_seed are known
    out_dir = output_dir or "."
    os.makedirs(out_dir, exist_ok=True)
    link_path = os.path.join(out_dir, f"{output_prefix}Coord_Links.csv")
    try:
        _write_coord_links_csv(
            link_path,          # <-- positional instead of path="..."
            cofactor_atoms=cofactor_sphere,
            pcs_seed_atoms=pcs_seed,
            scs_seed_atoms=scs_seed,
        )
    except Exception as e:
        print(f"[WARNING] Failed to write Coord_Links.csv ({link_path}): {e}")



    # --- Expansion as before (can include carbons etc.) ---
    if expand_residues:
        pcs_atoms = _expand_by_moiety_and_backbone_with_logging(pcs_seed, structure, exclude_moieties=exclude_moieties)
        scs_atoms = _expand_by_moiety_and_backbone_with_logging(scs_seed, structure, exclude_moieties=exclude_moieties)
    else:
        pcs_atoms, scs_atoms = pcs_seed, scs_seed

    print(f"[INFO] PCS (final): {len(pcs_atoms)} atoms")
    print(f"[INFO] SCS (final): {len(scs_atoms)} atoms")

    return cofactor_sphere, pcs_atoms, scs_atoms




# def identify_coordination_network(
#     structure,
#     cofactor_resname: List[str],
#     distance_cutoff: float,
#     expand_residues: bool,
#     combinatorial_mode: bool,
#     combinatorial_cofactor_cutoff: Optional[float] = None,
#     cofactor_resname2: Optional[List[str]] = None,
#     exclude_moieties: Optional[List[str]] = None,
# ):
#     """
#     Unchanged signature. Now excludes atoms with element 'C' from being selected
#     as PCS/SCS *seed* coordinators. Expansion (if enabled) can still include carbons.
#     """
#     # ---- Gather atoms from structure (you likely already have a helper for this) ----
#     all_atoms: List[AtomDict] = []
#     for model in structure:
#         for chain in model:
#             chain_id = chain.id
#             for residue in chain:
#                 resname = residue.get_resname()
#                 resnum  = residue.get_id()[1]
#                 for atom in residue:
#                     all_atoms.append({
#                         "name": atom.get_name(),
#                         "residue": resname,
#                         "residue_number": resnum,
#                         "chain": chain_id,
#                         "element": getattr(atom, "element", ""),
#                         "coordinates": np.array(atom.coord, dtype=float)
#                     })

#     # ---- Build cofactor atom set ----
#     cof_resnames = set(x.upper() for x in (cofactor_resname or []))
#     if cofactor_resname2:
#         cof_resnames |= set(x.upper() for x in cofactor_resname2)

#     cofactor_sphere = [a for a in all_atoms if str(a["residue"]).upper() in cof_resnames]
#     cof_coords = _np_coords(cofactor_sphere)

#     print(f"[INFO] Cofactor (identified): {len(cofactor_sphere)} atoms")
#     for a in cofactor_sphere:
#         print(f"[INFO]   Cofactor ▸ {a['residue']} {a['residue_number']} {a['chain']}: {a['name']}")

#     # ---- PCS seeds (exclude C; moiety-wise pick one; distance to cofactor) ----
#     # Exclude residues that are themselves cofactor
#     exclude_res_keys = {_rkey(a) for a in cofactor_sphere}

#     pcs_candidates = _filter_seed_candidates(
#         atoms=[a for a in all_atoms if str(a["residue"]).upper() not in cof_resnames],
#         exclude_residue_keys=exclude_res_keys,
#         max_dist_from_set=distance_cutoff,
#         reference_set_coords=cof_coords,
#     )
#     pcs_seed = _choose_one_per_moiety_group(pcs_candidates, reference_set_coords=cof_coords)

#     print(f"[INFO] PCS seeds (post C-exclusion, moiety-wise): {len(pcs_seed)} atoms")
#     # Log by residue
#     pcs_by_res = defaultdict(list)
#     for a in pcs_seed: pcs_by_res[_rkey(a)].append(a["name"])
#     for (res, num, ch), names in sorted(pcs_by_res.items()):
#         print(f"[INFO]   PCS seed ▸ {res} {num} {ch}: {', '.join(sorted(names))}")

#     # ---- SCS seeds (exclude C; moiety-wise pick one; distance to PCS seeds; exclude cofactor & PCS residues) ----
#     pcs_coords = _np_coords(pcs_seed)
#     exclude_res_keys_for_scs = exclude_res_keys | {_rkey(a) for a in pcs_seed}

#     scs_candidates = _filter_seed_candidates(
#         atoms=[a for a in all_atoms if str(a["residue"]).upper() not in cof_resnames],
#         exclude_residue_keys=exclude_res_keys_for_scs,
#         max_dist_from_set=distance_cutoff,       # SCS distance is to PCS (not 2× cutoff)
#         reference_set_coords=pcs_coords,
#     )
#     scs_seed = _choose_one_per_moiety_group(scs_candidates, reference_set_coords=pcs_coords)

#     print(f"[INFO] SCS seeds (post C-exclusion, moiety-wise): {len(scs_seed)} atoms")
#     scs_by_res = defaultdict(list)
#     for a in scs_seed: scs_by_res[_rkey(a)].append(a["name"])
#     for (res, num, ch), names in sorted(scs_by_res.items()):
#         print(f"[INFO]   SCS seed ▸ {res} {num} {ch}: {', '.join(sorted(names))}")

#     # ---- Expansion (unchanged; can include carbons) ----
#     if expand_residues:
#         # You already have your expansion helpers; call them as before
#         pcs_atoms = _expand_by_moiety_and_backbone_with_logging(pcs_seed, structure, exclude_moieties=exclude_moieties)
#         scs_atoms = _expand_by_moiety_and_backbone_with_logging(scs_seed, structure, exclude_moieties=exclude_moieties)
#     else:
#         pcs_atoms, scs_atoms = pcs_seed, scs_seed

#     # Final summaries (unchanged)
#     print(f"[INFO] PCS (final): {len(pcs_atoms)} atoms")
#     print(f"[INFO] SCS (final): {len(scs_atoms)} atoms")

#     return cofactor_sphere, pcs_atoms, scs_atoms





# def identify_coordination_network(
#     structure,
#     cofactor_resname: List[str],
#     distance_cutoff: float,
#     expand_residues: bool,
#     combinatorial_mode: bool,
#     combinatorial_cofactor_cutoff: Optional[float] = None,
#     cofactor_resname2: Optional[List[str]] = None,
#     exclude_moieties: Optional[List[str]] = None,
# ):


#     """
#     Returns: (cofactor_sphere, pcs_atoms, scs_atoms)

#     • PCS = one atom per (residue,moiety), chosen by closest distance to final cofactor set.
#     • SCS = one atom per (residue,moiety), chosen by closest distance to the selected PCS
#       coordinators only; excludes anything in cofactor or PCS.
#     • Optional expansion: add all atoms of those moieties + backbone for their residues.
#     """
#     if exclude_moieties is None:
#         exclude_moieties = []

#     primary_names = {r.upper() for r in (cofactor_resname or [])}
#     secondary_names = {r.upper() for r in (cofactor_resname2 or [])}

#     print(f"[INFO] Identifying coordination network for {sorted(primary_names)}")

#     # --- Gather all atoms as dicts ---
#     all_atoms: List[AtomDict] = []
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 for atom in residue:
#                     # optional: skip hydrogens here entirely to speed up
#                     elem = getattr(atom, "element", "")
#                     if str(elem).upper() in {"H", "D"}:
#                         continue
#                     all_atoms.append(_bioatom_to_dict(residue, atom, chain.id))

#     # --- Cofactor atoms (primary) ---
#     cof_primary = [a for a in all_atoms if str(a["residue"]).upper() in primary_names]

#     # --- Combinatorial extension (optional) ---
#     cofactor_sphere: List[AtomDict] = list(cof_primary)
#     if combinatorial_mode and secondary_names:
#         if not combinatorial_cofactor_cutoff:
#             combinatorial_cofactor_cutoff = max(3.0, distance_cutoff)
#         pcoords = np.array([a["coordinates"] for a in cof_primary]) if cof_primary else np.empty((0, 3))
#         added = 0
#         if pcoords.size:
#             for a in all_atoms:
#                 if str(a["residue"]).upper() not in secondary_names:
#                     continue
#                 d = _min_dist_to_cloud(a["coordinates"], pcoords)
#                 if d <= combinatorial_cofactor_cutoff:
#                     cofactor_sphere.append(a)
#                     added += 1
#         cofactor_sphere = _dedup_atoms(cofactor_sphere)
#         print(f"[INFO] Combinatorial cofactor extension: added {added} atoms from {sorted(secondary_names)} within {combinatorial_cofactor_cutoff:.2f} Å")
#     elif combinatorial_mode and not secondary_names:
#         print("[INFO] combinatorial_mode=True but no cofactor_resname2 provided; skipping combinatorial extension.")

#     _print_atom_set("Cofactor (identified)", cofactor_sphere)

#     cof_coords = np.array([a["coordinates"] for a in cofactor_sphere]) if cofactor_sphere else np.empty((0, 3))

#     # ------------------------------------------------
#     # PCS seeds: per (res,moiety) closest to COFACTOR
#     # ------------------------------------------------
#     pcs_candidates: Dict[Tuple[str, int, str, str], Candidate] = {}
#     for a in all_atoms:
#         # exclude cofactor atoms themselves
#         if str(a["residue"]).upper() in primary_names or str(a["residue"]).upper() in secondary_names:
#             continue
#         d = _min_dist_to_cloud(a["coordinates"], cof_coords)
#         if d <= float(distance_cutoff):
#             moiety = _get_moiety_label(str(a["residue"]), str(a["name"]))
#             key = (str(a["residue"]), int(a["residue_number"]), str(a["chain"]), moiety)
#             _keep_closest(pcs_candidates, key, a, d, moiety)

#     # Freeze PCS seeds (one atom per (res,moiety))
#     pcs_seed_items = sorted(pcs_candidates.items(), key=lambda kv: kv[1].distance)
#     pcs_seed: List[AtomDict] = [cand.atom for (_, cand) in pcs_seed_items]

#     # Log PCS seeds with distances
#     for (resname, resnum, chain, moi), cand in pcs_seed_items:
#         print(f"[INFO] PCS (seed) ▸ {resname} {resnum} {chain}, moiety={moi}: {cand.atom['name']} ({cand.distance:.2f} Å to cofactor)")

#     # ------------------------------------------------
#     # SCS seeds: per (res,moiety) closest to **PCS**
#     #           exclude cofactor and any (res,moiety) in PCS
#     # ------------------------------------------------
#     pcs_coords = np.array([a["coordinates"] for a in pcs_seed]) if pcs_seed else np.empty((0, 3))
#     pcs_groups = {(a["residue"], a["residue_number"], a["chain"], _get_moiety_label(str(a["residue"]), str(a["name"]))) for a in pcs_seed}
#     cof_keys = {(a["residue"], a["residue_number"], a["chain"], a["name"]) for a in cofactor_sphere}

#     # choose a reasonable SCS cutoff; you can pass another param if you prefer
#     scs_cutoff = float(distance_cutoff) * 1.5  # e.g., a bit larger than PCS→cofactor
#     scs_candidates: Dict[Tuple[str, int, str, str], Candidate] = {}

#     if pcs_coords.size:
#         for a in all_atoms:
#             # exclude cofactor atoms (by individual atoms) and atoms whose (res,moiety) already appears in PCS
#             if (a["residue"], a["residue_number"], a["chain"], a["name"]) in cof_keys:
#                 continue
#             moiety = _get_moiety_label(str(a["residue"]), str(a["name"]))
#             g = (str(a["residue"]), int(a["residue_number"]), str(a["chain"]), moiety)
#             if g in pcs_groups:
#                 continue  # do not allow same (res,moiety) to appear in SCS

#             d = _min_dist_to_cloud(a["coordinates"], pcs_coords)  # distance to PCS coordinators
#             if d <= scs_cutoff:
#                 _keep_closest(scs_candidates, g, a, d, moiety)

#     scs_seed_items = sorted(scs_candidates.items(), key=lambda kv: kv[1].distance)
#     scs_seed: List[AtomDict] = [cand.atom for (_, cand) in scs_seed_items]

#     # Log SCS seeds with distances
#     for (resname, resnum, chain, moi), cand in scs_seed_items:
#         print(f"[INFO] SCS (seed) ▸ {resname} {resnum} {chain}, moiety={moi}: {cand.atom['name']} ({cand.distance:.2f} Å to PCS)")

#     # -------------------------------------------
#     # Optional expansion (moiety + backbone)
#     # -------------------------------------------
#     if expand_residues:
#         pcs_atoms = _expand_by_moiety_and_backbone_with_logging(pcs_seed, structure, exclude_moieties=exclude_moieties)
#         scs_atoms = _expand_by_moiety_and_backbone_with_logging(scs_seed, structure, exclude_moieties=exclude_moieties)
#     else:
#         pcs_atoms, scs_atoms = pcs_seed, scs_seed

#     # Final summaries
#     _print_atom_set("Cofactor (final)", cofactor_sphere)
#     _print_atom_set("PCS (final)", pcs_atoms)
#     _print_atom_set("SCS (final)", scs_atoms)


#     # -------------------------------------------
#     # NEW: write linked edge CSV from seeds
#     # -------------------------------------------
#     _write_coord_links_csv(
#         outfile="Coord_Links.csv",
#         cofactor_atoms=cofactor_sphere,
#         pcs_seeds=pcs_seed,
#         scs_seeds=scs_seed,
#     )

#     return cofactor_sphere, pcs_atoms, scs_atoms





# ensure it’s exported
__all__ = [name for name in globals().keys() if name in ("find_cofactor_atoms",)]
