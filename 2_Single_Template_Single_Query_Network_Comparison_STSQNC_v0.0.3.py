#!/usr/bin/env python3
"""
v0.0.3 (STSQNC): validates accessory names, supports two-section or single-section
template CSVs, aligns template atoms to query using Kabsch, writes summary +
atom-level alignment files, and adds a network match table that maps every
template network atom to its nearest query atom (distance in Å).
"""
import argparse
import csv
import os
import sys
from typing import Dict, Any, Tuple, List, Optional
import numpy as np




# Ensure project modules import
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from modules.io_utils import unpack_pdb_file


def load_template_csv(path: str) -> Dict[str, Any]:
    """
    Load a coordination CSV.
    Supports:
      - two-section format (residue summary, blank, atom section with Atom column)
      - single-section per-atom format with header containing Atom and x,y,z
    """
    with open(path, newline="") as f:
        rows = list(csv.reader(f))

    if not rows:
        raise ValueError("Template CSV is empty")

    header_idx = None
    atom_section_idx = None
    # find first header row
    for idx, row in enumerate(rows):
        if row and row[0] == "Residue Category":
            header_idx = idx
            break
    if header_idx is None:
        raise ValueError("Could not find header in template CSV")

    header_row = rows[header_idx]
    has_atom_in_header = "Atom" in header_row
    # if first header has no Atom col, search for second header with Atom
    if not has_atom_in_header:
        for idx in range(header_idx + 1, len(rows)):
            row = rows[idx]
            if row and row[0] == "Residue Category" and "Atom" in row:
                atom_section_idx = idx
                break
        if atom_section_idx is None:
            raise ValueError("Template CSV lacks atom section; regenerate with atom coordinates before alignment.")
    else:
        atom_section_idx = header_idx

    # Count rows
    if atom_section_idx == header_idx:
        # single-section atom table
        residue_rows = 0
        atom_rows = max(0, len(rows) - atom_section_idx - 1)
    else:
        residue_rows = max(0, atom_section_idx - header_idx - 1)
        atom_rows = max(0, len(rows) - atom_section_idx - 1)

    return {
        "residue_rows": residue_rows,
        "atom_rows": atom_rows,
        "has_atom_section": True,
        "rows": rows,
        "header_idx": header_idx,
        "atom_section_idx": atom_section_idx,
        "single_atom_section": atom_section_idx == header_idx,
    }


def parse_accessory(path: str) -> Dict[str, List[List[str]]]:
    """Read accessory_names.txt into residue and atom synonym groups."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Accessory file not found: {path}")

    residues: List[List[str]] = []
    atoms: List[List[str]] = []
    section: Optional[str] = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower() == "[residuenames]":
                section = "residues"
                continue
            if line.lower() == "[atomnames]":
                section = "atoms"
                continue
            parts = line.split()
            if not parts:
                continue
            if section == "residues":
                residues.append(parts)
            elif section == "atoms":
                atoms.append(parts)
    if not residues:
        raise ValueError("No residue names found in accessory file")
    if not atoms:
        raise ValueError("No atom names found in accessory file")
    return {"residues": residues, "atoms": atoms}


def _build_synonym_map(groups: List[List[str]]) -> Dict[str, str]:
    """Map every synonym (uppercased) to the group's canonical (first) name uppercased."""
    mapping: Dict[str, str] = {}
    for group in groups:
        if not group:
            continue
        canon = group[0].upper()
        for name in group:
            mapping[name.upper()] = canon
    return mapping


def _collect_template_atoms(rows: List[List[str]], atom_section_idx: int, single_atom_section: bool, res_map: Dict[str, str], atom_map: Dict[str, str]):
    """
    Collect template atoms keyed by residue.
    Returns { (canon_res, resnum, chain): {canon_atom: coord(np.array)} }
    """
    by_res: Dict[Tuple[str, str, str], Dict[str, np.ndarray]] = {}
    start_idx = atom_section_idx + 1
    for row in rows[start_idx:]:
        if not row:
            continue
        # Determine indices depending on format
        if single_atom_section:
            # Expected header: Category,Atom,Residue,Residue Number,Chain,Moiety,InteractorFlag,x,y,z
            if len(row) < 10:
                continue
            resname_raw = row[2].upper()
            atom_raw = row[1].upper()
            resnum = str(row[3])
            chain = str(row[4])
            coord_slice = row[7:10]
        else:
            # Two-section atom header: Residue Category,Residue Name,Residue Number,Chain,Atom,x,y,z
            if len(row) < 7:
                continue
            resname_raw = row[1].upper()
            atom_raw = row[4].upper()
            resnum = str(row[2])
            chain = str(row[3])
            coord_slice = row[5:8]

        if resname_raw not in res_map or atom_raw not in atom_map:
            continue
        canon_res = res_map[resname_raw]
        canon_atom = atom_map[atom_raw]
        try:
            x, y, z = map(float, coord_slice)
        except Exception:
            continue
        key = (canon_res, resnum, chain)
        by_res.setdefault(key, {})
        by_res[key][canon_atom] = np.array([x, y, z], dtype=float)
    return by_res


def _collect_query_atoms(structure, res_map: Dict[str, str], atom_map: Dict[str, str]):
    """Collect atoms from query structure keyed by residue."""
    by_res: Dict[Tuple[str, str, str], Dict[str, np.ndarray]] = {}
    for model in structure:
        for chain in model:
            chain_id = str(chain.id)
            for residue in chain:
                resname_raw = residue.get_resname().upper()
                if resname_raw not in res_map:
                    continue
                canon_res = res_map[resname_raw]
                resnum = str(residue.get_id()[1])
                key = (canon_res, resnum, chain_id)
                for atom in residue:
                    atom_raw = atom.get_name().upper()
                    if atom_raw not in atom_map:
                        continue
                    canon_atom = atom_map[atom_raw]
                    by_res.setdefault(key, {})
                    by_res[key][canon_atom] = np.array(atom.coord, dtype=float)
    return by_res


def _collect_all_query_atoms(structure) -> List[Dict[str, Any]]:
    """Collect all non-water atoms from the query structure."""
    water_names = {"HOH", "WAT", "H2O"}
    atoms: List[Dict[str, Any]] = []
    for model in structure:
        for chain in model:
            chain_id = str(chain.id)
            for residue in chain:
                resname = residue.get_resname()
                if str(resname).upper() in water_names:
                    continue
                resnum = residue.get_id()[1]
                for atom in residue:
                    atoms.append({
                        "residue": str(resname).upper(),
                        "residue_number": str(resnum),
                        "chain": chain_id,
                        "atom": str(atom.get_name()).upper(),
                        "coords": np.array(atom.coord, dtype=float),
                        "element": getattr(atom, "element", ""),
                    })
    return atoms


def _infer_atoms_csv_path(template_csv: str) -> str:
    """
    Infer path to the per-atom template CSV (Coord_Breakdown_atoms).
    Example: foo/4ub6.pdb_Coord_Breakdown.csv -> foo/4ub6.pdb_Coord_Breakdown_atoms.csv
    """
    base, ext = os.path.splitext(template_csv)
    if base.endswith("_Coord_Breakdown"):
        return base + "_atoms" + ext
    return base + "_atoms" + ext


def _load_template_atoms(atoms_csv_path: str) -> List[Dict[str, Any]]:
    """Load template network atoms from *_Coord_Breakdown_atoms.csv."""
    if not os.path.isfile(atoms_csv_path):
        raise FileNotFoundError(f"Template atoms file not found: {atoms_csv_path}")
    atoms: List[Dict[str, Any]] = []
    with open(atoms_csv_path, newline="") as f:
        reader = csv.DictReader(f)
        required = {"Category", "Atom", "Residue", "Residue Number", "Chain", "x", "y", "z"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Atoms CSV missing columns: {missing}")
        for row in reader:
            try:
                x = float(row["x"]); y = float(row["y"]); z = float(row["z"])
            except Exception:
                continue
            atoms.append({
                "category": row.get("Category", ""),
                "atom": row.get("Atom", "").upper(),
                "residue": row.get("Residue", "").upper(),
                "residue_number": str(row.get("Residue Number", "")),
                "chain": str(row.get("Chain", "")),
                "moiety": row.get("Moiety", ""),
                "interactor": row.get("InteractorFlag", ""),
                "coord": np.array([x, y, z], dtype=float),
            })
    return atoms


def _pick_best_residue(res_atoms: Dict[Tuple[str, str, str], Dict[str, np.ndarray]]) -> Optional[Tuple[Tuple[str, str, str], Dict[str, np.ndarray]]]:
    """Pick residue with the most mapped atoms."""
    best = None
    best_len = -1
    for key, atoms in res_atoms.items():
        l = len(atoms)
        if l > best_len:
            best = (key, atoms)
            best_len = l
    return best


def _kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return rotation matrix R and translation t such that R@P.T + t aligns P to Q."""
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc
    H = P0.T @ Q0
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = Qc - R @ Pc
    return R, t


def _rmsd(P_aligned: np.ndarray, Q: np.ndarray) -> float:
    diff = P_aligned - Q
    return float(np.sqrt((diff * diff).sum() / len(P_aligned)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Align template coordination CSV to a query structure and report per-atom matches."
    )
    parser.add_argument("--template-csv", required=True, help="Path to *_Coord_Breakdown.csv file.")
    parser.add_argument("--template-atoms-csv", required=False, help="Path to *_Coord_Breakdown_atoms.csv (defaults to inferred sibling).")
    parser.add_argument("--query-structure", required=True, help="Path to query PDB or CIF file.")
    parser.add_argument(
        "--job-name",
        required=True,
        help="Job name (creates/uses matches/<job-name> as working directory)."
    )
    parser.add_argument("--match-cutoff", type=float, default=2.0, help="Distance cutoff (Å) to consider a template atom matched to the query.")
    return parser.parse_args()


def ensure_job_dir(job_name: str) -> str:
    """Create matches/<job_name> and return its path."""
    base = os.path.join(os.getcwd(), "matches")
    job_dir = os.path.join(base, job_name)
    os.makedirs(job_dir, exist_ok=True)
    return job_dir


def write_accessory_file(job_dir: str) -> str:
    """
    Write accessory file with interchangeable residue names and atom names.
    Residue section: same-line synonyms (OEX and MY-OEX).
    Atom section: atom names, with synonyms on the same line (CA1 and CAL).
    """
    residue_groups: List[List[str]] = [["OEX", "MY-OEX"]]
    atom_groups: List[List[str]] = [
        ["O1"], ["O2"], ["O3"], ["O4"], ["O5"],
        ["MN1"], ["MN2"], ["MN3"], ["MN4"],
        ["CA1", "CAL"],
    ]

    path = os.path.join(job_dir, "accessory_names.txt")
    with open(path, "w", newline="") as f:
        f.write("# Allowed residue and atom names (synonyms on same line)\n")
        f.write("[ResidueNames]\n")
        for group in residue_groups:
            f.write(" ".join(group) + "\n")
        f.write("\n[AtomNames]\n")
        for group in atom_groups:
            f.write(" ".join(group) + "\n")
    return path


def main() -> None:
    args = parse_args()

    if not os.path.isfile(args.template_csv):
        raise FileNotFoundError(f"Template CSV not found: {args.template_csv}")
    if not os.path.isfile(args.query_structure):
        raise FileNotFoundError(f"Query structure not found: {args.query_structure}")

    job_dir = ensure_job_dir(args.job_name)
    accessory_path = write_accessory_file(job_dir)

    tpl_info = load_template_csv(args.template_csv)
    accessory = parse_accessory(accessory_path)
    res_map = _build_synonym_map(accessory["residues"])
    atom_map = _build_synonym_map(accessory["atoms"])

    template_atoms_all = _collect_template_atoms(
        tpl_info["rows"],
        tpl_info["atom_section_idx"],
        tpl_info["single_atom_section"],
        res_map,
        atom_map
    )
    tpl_pick = _pick_best_residue(template_atoms_all)
    if tpl_pick is None:
        raise ValueError("No template atoms matched allowed residue/atom names.")
    tpl_res_key, tpl_atoms = tpl_pick

    structure, atoms = unpack_pdb_file(args.query_structure)
    atom_count = len(atoms) if atoms is not None else 0

    query_atoms_all = _collect_query_atoms(structure, res_map, atom_map)
    qry_pick = _pick_best_residue(query_atoms_all)
    if qry_pick is None:
        raise ValueError("No query atoms matched allowed residue/atom names.")
    qry_res_key, qry_atoms = qry_pick

    # Pair atoms by canonical name
    common_atoms = sorted(set(tpl_atoms.keys()) & set(qry_atoms.keys()))
    if len(common_atoms) < 3:
        raise ValueError(f"Need at least 3 shared atoms for alignment; found {len(common_atoms)} ({common_atoms}).")

    P = np.vstack([tpl_atoms[a] for a in common_atoms])  # template coords
    Q = np.vstack([qry_atoms[a] for a in common_atoms])  # query coords
    R, t = _kabsch(P, Q)
    P_aligned = (R @ P.T).T + t
    rmsd_val = _rmsd(P_aligned, Q)

    # Load full template network atoms (Coord_Breakdown_atoms) and map to query
    atoms_csv_path = args.template_atoms_csv or _infer_atoms_csv_path(args.template_csv)
    template_network_atoms = _load_template_atoms(atoms_csv_path)
    query_atoms_full = _collect_all_query_atoms(structure)
    query_coords = np.vstack([a["coords"] for a in query_atoms_full]) if query_atoms_full else np.empty((0, 3))
    match_cutoff = float(args.match_cutoff)

    total_atoms = len(template_network_atoms)
    matched_total = 0
    category_stats: Dict[str, Dict[str, int]] = {}
    unmatched_atoms: List[str] = []

    network_path = os.path.join(job_dir, f"{args.job_name}_network_match.csv")
    with open(network_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "job_name",
            "category",
            "template_residue", "template_resnum", "template_chain", "template_atom", "template_moiety", "template_interactor",
            "query_residue", "query_resnum", "query_chain", "query_atom", "query_element",
            "distance_A",
            "same_residue_name", "same_atom_name",
            "matched_within_cutoff",
        ])
        for ta in template_network_atoms:
            tpl_aligned_coord = (R @ ta["coord"]) + t
            if query_coords.size:
                diffs = query_coords - tpl_aligned_coord
                d2 = np.einsum("ij,ij->i", diffs, diffs)
                idx = int(np.argmin(d2))
                min_dist = float(np.sqrt(d2[idx]))
                qa = query_atoms_full[idx]
            else:
                min_dist = float("nan")
                qa = None
            cat = ta["category"] or "UNK"
            category_stats.setdefault(cat, {"total": 0, "matched": 0})
            category_stats[cat]["total"] += 1
            is_matched = (min_dist == min_dist) and (min_dist <= match_cutoff)
            if is_matched:
                matched_total += 1
                category_stats[cat]["matched"] += 1
            else:
                unmatched_atoms.append(f"{ta['chain']}:{ta['residue']}{ta['residue_number']}:{ta['atom']} (dist={min_dist:.2f} Å)")
            writer.writerow([
                args.job_name,
                ta["category"],
                ta["residue"], ta["residue_number"], ta["chain"], ta["atom"], ta["moiety"], ta["interactor"],
                qa["residue"] if qa else "",
                qa["residue_number"] if qa else "",
                qa["chain"] if qa else "",
                qa["atom"] if qa else "",
                qa["element"] if qa else "",
                f"{min_dist:.4f}" if min_dist == min_dist else "",
                "Y" if qa and qa["residue"] == ta["residue"] else "N",
                "Y" if qa and qa["atom"] == ta["atom"] else "N",
                "Y" if is_matched else "N",
            ])

    # Summary markdown
    percent = (matched_total / total_atoms * 100.0) if total_atoms else 0.0
    summary_lines: List[str] = []
    summary_lines.append(f"# Alignment Summary — {args.job_name}")
    summary_lines.append("")
    summary_lines.append(f"- Template CSV: {args.template_csv}")
    summary_lines.append(f"- Template atoms: {atoms_csv_path}")
    summary_lines.append(f"- Query structure: {args.query_structure}")
    summary_lines.append(f"- Match cutoff: {match_cutoff:.2f} Å")
    summary_lines.append("")
    summary_lines.append("## Overall")
    summary_lines.append(f"- Matched {matched_total}/{total_atoms} network atoms ({percent:.1f}%).")
    summary_lines.append("")
    summary_lines.append("## By Category")
    summary_lines.append("| Category | Matched / Total | Percent |")
    summary_lines.append("| --- | --- | --- |")
    for cat, stats in sorted(category_stats.items()):
        tot = stats["total"]
        mat = stats["matched"]
        pct = (mat / tot * 100.0) if tot else 0.0
        summary_lines.append(f"| {cat} | {mat} / {tot} | {pct:.1f}% |")
    summary_lines.append("")
    summary_lines.append("## Unmatched template atoms (within cutoff = N)")
    if unmatched_atoms:
        for u in unmatched_atoms:
            summary_lines.append(f"- {u}")
    else:
        summary_lines.append("- None")

    summary_path = os.path.join(job_dir, f"{args.job_name}_results.md")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary_lines))

    # Write atom-by-atom alignment matrix (one file per run)
    matrix_path = os.path.join(job_dir, f"{args.job_name}_alignment_matrix.csv")
    with open(matrix_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "job_name",
            "atom_name",
            "template_residue", "template_resnum", "template_chain",
            "query_residue", "query_resnum", "query_chain",
            "tpl_aligned_x", "tpl_aligned_y", "tpl_aligned_z",
            "qry_x", "qry_y", "qry_z",
            "per_atom_rmsd"
        ])
        for idx, atom_name in enumerate(common_atoms):
            tpl_coords_aligned = P_aligned[idx]
            qry_coords = Q[idx]
            dist = float(np.linalg.norm(tpl_coords_aligned - qry_coords))
            writer.writerow([
                args.job_name,
                atom_name,
                tpl_res_key[0], tpl_res_key[1], tpl_res_key[2],
                qry_res_key[0], qry_res_key[1], qry_res_key[2],
                f"{tpl_coords_aligned[0]:.4f}", f"{tpl_coords_aligned[1]:.4f}", f"{tpl_coords_aligned[2]:.4f}",
                f"{qry_coords[0]:.4f}", f"{qry_coords[1]:.4f}", f"{qry_coords[2]:.4f}",
                f"{dist:.4f}",
            ])

    report_path = os.path.join(job_dir, "alignment_report.csv")
    write_header = not os.path.isfile(report_path)
    with open(report_path, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow([
                "job_name",
                "template_residue", "template_resnum", "template_chain",
                "query_residue", "query_resnum", "query_chain",
                "n_atoms", "atoms_used", "rmsd"
            ])
        writer.writerow([
            args.job_name,
            tpl_res_key[0], tpl_res_key[1], tpl_res_key[2],
            qry_res_key[0], qry_res_key[1], qry_res_key[2],
            len(common_atoms),
            " ".join(common_atoms),
            rmsd_val,
        ])

    print("[OK] Template CSV loaded.")
    print(f"      Residue rows: {tpl_info['residue_rows']}")
    print(f"      Atom rows: {tpl_info['atom_rows']}")
    print(f"      Format: {'single atom section' if tpl_info['single_atom_section'] else 'two-section'}")
    print("[OK] Accessory names parsed.")
    print(f"      Allowed residues: {accessory['residues']}")
    print(f"      Allowed atoms: {accessory['atoms']}")
    print("[OK] Query structure loaded.")
    print(f"      Model count: {len(list(structure)) if structure is not None else 0}")
    print(f"      Atom count: {atom_count}")
    print("[OK] Job directory prepared.")
    print(f"      Job dir: {job_dir}")
    print(f"      Accessory names: {accessory_path}")
    print("[OK] Alignment complete.")
    print(f"      Atoms used: {common_atoms}")
    print(f"      RMSD: {rmsd_val:.4f}")
    print(f"      Match cutoff: {match_cutoff:.2f} Å")
    print(f"      Network matched: {matched_total}/{total_atoms} ({percent:.1f}%)")
    print(f"[OK] Alignment matrix: {matrix_path}")
    print(f"[OK] Network match table: {network_path}")
    print(f"[OK] Results summary: {summary_path}")
    print(f"[OK] Report: {report_path}")
    print("[DONE] Load + validation + alignment successful.")


if __name__ == "__main__":
    main()
