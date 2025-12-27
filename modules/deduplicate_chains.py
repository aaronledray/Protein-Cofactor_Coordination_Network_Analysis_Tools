

#!/usr/bin/env python3
"""
Deduplicate chains in PDB / ENT / mmCIF structures, with optional **mask-based** canonical selection.

- Detects duplicate chains by **sequence** (default) or by **exact atom coordinate signature** (`--by coords`).
- Optional **mask** lets you keep *one* replicate of the chain that contains specific residues/ligands (e.g., HEM and CU1).
- If **multiple chains match the mask** after de-duplication, the tool can **interactively reduce to a single chain** (defaults to the lexicographically smallest chain id, e.g., 'A').
- Works on a **single file**, or **batch** for a directory (recursive).
- Interactive CLI by default; use `--yes` for unattended batch runs.
- Always writes a new **PDB** file, even if input was `.cif` or `.ent`.

Examples
--------
Single file (interactive):
    python dedupe_chains.py 5jqr.cif

Single file, pick the chain that contains HEM *and* CU1, auto-yes:
    python dedupe_chains.py 5jqr.cif --mask-res HEM,CU1 --mask-mode all --yes

Batch a directory (recursive), coordinate-based signature, prefer chain with HEM, auto-yes:
    python dedupe_chains.py ./structures --by coords --mask-res HEM --yes

Notes
-----
- Only the **first model** (model id 0) is processed/written.
- "Duplicate" means chains share the same signature according to the chosen method.
  * sequence: tuple of residue 3-letter names for **standard** residues (waters/hetero excluded)
  * coords:   ordered list of (residue id, atom name, rounded coords) across all **non-water** residues
- For masks, residue names are matched **case-insensitively** against the residue name reported by Biopython (includes hetero residues like HEM, CU1). Waters (HOH/WAT/H2O) are ignored.
- For mmCIF, the chain identifier is the **asym_id** as parsed by Biopython.
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Optional, Set

try:
    from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
except Exception as e:
    print("[ERROR] This tool requires Biopython. Install with: pip install biopython", file=sys.stderr)
    raise


WATER_NAMES = {"HOH", "WAT", "H2O"}


def is_water(res: Residue) -> bool:
    if res.id[0] == 'W':
        return True
    return res.get_resname().strip().upper() in WATER_NAMES


def is_standard_amino_or_nuc(res: Residue) -> bool:
    # Biopython marks standard residues with hetflag ' '
    return res.id[0] == ' '


def chain_sequence_signature(chain: Chain) -> Tuple[str, ...]:
    """Return a tuple of residue 3-letter names for **standard** residues (excludes waters/hetero)."""
    sig: List[str] = []
    for res in chain.get_residues():
        if not is_standard_amino_or_nuc(res):
            continue
        if is_water(res):
            continue
        sig.append(res.get_resname().strip().upper())
    return tuple(sig)


def chain_coordinate_signature(chain: Chain, ndp: int = 3) -> Tuple[Tuple, ...]:
    """Coordinate-based signature for a chain.

    The signature is an ordered tuple of entries:
        ((resseq, icode, resname), atom_name, (x,y,z))
    with coordinates rounded to `ndp` decimals, excluding **waters only** (hetero kept).
    Atoms are sorted stably by (resseq, icode, atom_name) to ensure determinism.
    """
    entries: List[Tuple] = []
    for res in chain.get_residues():
        if is_water(res):
            continue
        resseq = res.id[1]
        icode = res.id[2].strip() if res.id[2] else ''
        resname = res.get_resname().strip().upper()
        for atom in res.get_atoms():
            x, y, z = atom.coord
            entries.append(((resseq, icode, resname), atom.get_name().strip(), (round(float(x), ndp), round(float(y), ndp), round(float(z), ndp))))
    entries.sort(key=lambda t: (t[0][0], t[0][1], t[1]))
    return tuple(entries)


def parse_structure(path: Path) -> Structure:
    ext = path.suffix.lower()
    if ext in {".pdb", ".ent"}:
        parser = PDBParser(QUIET=True)
    elif ext == ".cif" or ext == ".mmcif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported file extension: {ext}. Use .pdb, .ent, or .cif")
    return parser.get_structure(path.stem, str(path))


class ChainKeepSelect(Select):
    def __init__(self, model_id_to_keep: int, chain_ids_to_keep: Set[str]):
        super().__init__()
        self.model_id_to_keep = model_id_to_keep
        self.keep = chain_ids_to_keep

    def accept_model(self, model: Model) -> bool:
        return model.id == self.model_id_to_keep

    def accept_chain(self, chain: Chain) -> bool:
        return chain.id in self.keep


def chain_matches_mask(chain: Chain, mask_resnames: Set[str], mode: str = "all") -> bool:
    """Return True if the chain contains the masked residues by **residue name**.

    - mask_resnames: set of residue names to match (e.g., {"HEM","CU1"}). Case-insensitive.
    - mode: "all" requires all names present at least once; "any" requires at least one.
    - Waters are ignored; hetero residues (ligands, metals) are eligible.
    """
    if not mask_resnames:
        return False
    present: Set[str] = set()
    for res in chain.get_residues():
        if is_water(res):
            continue
        present.add(res.get_resname().strip().upper())
    if mode == "all":
        return mask_resnames.issubset(present)
    else:
        return len(mask_resnames.intersection(present)) > 0


def find_duplicate_groups(model: Model, method: str) -> List[List[str]]:
    """Return groups of chain IDs that are duplicates under the given method.
    The first chain in each group is considered the canonical one to keep (may be remapped later by mask preference).
    """
    if method not in {"sequence", "coords"}:
        raise ValueError("method must be 'sequence' or 'coords'")

    sig_map: Dict[Tuple, List[str]] = {}

    # Preserve encounter order of chains as they appear in the file
    for chain in model.get_chains():
        sig = chain_sequence_signature(chain) if method == "sequence" else chain_coordinate_signature(chain)
        sig_map.setdefault(sig, []).append(chain.id)

    groups = [ids for ids in sig_map.values() if len(ids) > 1]
    # Sort groups by the first chain id for determinism
    groups.sort(key=lambda g: g[0])
    return groups


def apply_mask_preference(model: Model, groups: List[List[str]], mask_resnames: Set[str], mask_mode: str) -> List[List[str]]:
    """Within each duplicate group, reorder so that the **canonical** (index 0) is the first chain that matches the mask.
    If none match, leave the group as-is. If multiple match, keep the earliest encountered.
    """
    if not mask_resnames:
        return groups

    # Map chain id -> Chain object for quick lookup
    chain_by_id: Dict[str, Chain] = {c.id: c for c in model.get_chains()}

    new_groups: List[List[str]] = []
    for grp in groups:
        chosen = None
        for cid in grp:
            ch = chain_by_id.get(cid)
            if ch and chain_matches_mask(ch, mask_resnames, mask_mode):
                chosen = cid
                break
        if chosen and chosen != grp[0]:
            ordered = [chosen] + [cid for cid in grp if cid != chosen]
            new_groups.append(ordered)
        else:
            new_groups.append(grp)
    return new_groups


def prompt_yes_no(msg: str, default: bool = True) -> bool:
    if not sys.stdin.isatty():
        return default
    suffix = "[Y/n]" if default else "[y/N]"
    while True:
        resp = input(f"{msg} {suffix} ").strip().lower()
        if resp == "" and default is not None:
            return default
        if resp in {"y", "yes"}:
            return True
        if resp in {"n", "no"}:
            return False
        print("Please answer y or n.")


def dedupe_file(path: Path, by: str, outdir: Optional[Path], auto_yes: bool, dry_run: bool,
                mask_res: Set[str], mask_mode: str) -> Tuple[Path, int]:
    structure = parse_structure(path)
    model0: Model = structure[0]

    groups = find_duplicate_groups(model0, by)
    groups = apply_mask_preference(model0, groups, mask_res, mask_mode)

    if not groups:
        print(f"[INFO] No duplicate chains found in {path.name} (method={by}).")
        out_path = build_outpath(path, outdir)
        if not dry_run:
            write_pdb(structure, keep_chains={c.id for c in model0.get_chains()}, out_path=out_path)
        return out_path, 0

    print(f"[INFO] {path.name}: found {sum(len(g)-1 for g in groups)} duplicate chain(s) across {len(groups)} group(s) (method={by}).")
    for i, grp in enumerate(groups, 1):
        print(f"  Group {i}: keep '{grp[0]}' (canonical), remove {grp[1:]}")

    proceed = auto_yes or prompt_yes_no("Remove duplicates and write output?", default=True)
    if not proceed:
        print("[SKIP] User chose not to modify this file.")
        return path, 0

    # Compute chains to keep: drop non-canonical duplicates
    duplicates = {cid for grp in groups for cid in grp[1:]}
    keep = [c.id for c in model0.get_chains() if c.id not in duplicates]

    # If a mask is provided and multiple remaining chains match it, offer interactive reduction to a single chain.
    if mask_res:
        # Identify keepers that match mask
        chain_by_id: Dict[str, Chain] = {c.id: c for c in model0.get_chains()}
        matching = [cid for cid in keep if chain_matches_mask(chain_by_id[cid], mask_res, mask_mode)]
        if len(matching) > 1:
            if auto_yes:
                # In non-interactive mode, do not collapse to one; keep all matches (safer for batch).
                print(f"[INFO] {path.name}: {len(matching)} chains match mask {sorted(mask_res)}; non-interactive mode keeps them all.")
            else:
                chosen = min(matching)  # lexicographically smallest, e.g., 'A'
                msg = (
                    f"Multiple chains match the mask ({', '.join(matching)}). "
                    f"Reduce to a single chain '{chosen}'?"
                )
                if prompt_yes_no(msg, default=True):
                    keep = [chosen]
                    print(f"[INFO] Reduced to single chain '{chosen}' per user choice.")

    out_path = build_outpath(path, outdir)
    if not dry_run:
        write_pdb(structure, keep_chains=set(keep), out_path=out_path)
    return out_path, len(duplicates)


def write_pdb(structure: Structure, keep_chains: Set[str], out_path: Path) -> None:
    io = PDBIO()
    io.set_structure(structure)
    sel = ChainKeepSelect(model_id_to_keep=0, chain_ids_to_keep=keep_chains)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(out_path), select=sel)
    print(f"[WRITE] {out_path}")


def build_outpath(src: Path, outdir: Optional[Path]) -> Path:
    base = src.stem
    out_dir = outdir if outdir else src.parent
    out_name = f"{base}_dedup.pdb"
    out_path = out_dir / out_name
    i = 1
    while out_path.exists() and out_path.resolve() != src.resolve():
        out_name = f"{base}_dedup({i}).pdb"
        out_path = out_dir / out_name
        i += 1
    return out_path


def gather_files(path: Path) -> List[Path]:
    if path.is_file():
        return [path]
    if path.is_dir():
        files: List[Path] = []
        for ext in ("*.pdb", "*.ent", "*.cif", "*.mmcif"):
            files.extend(path.rglob(ext))
        return sorted(files)
    raise FileNotFoundError(path)


def parse_mask_resnames(arg: Optional[str]) -> Set[str]:
    if not arg:
        return set()
    parts = [p.strip().upper() for p in arg.split(',') if p.strip()]
    return set(parts)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Remove duplicate chains from PDB/ENT/mmCIF files and write a new PDB.")
    p.add_argument("path", type=Path, help="Input file or directory")
    p.add_argument("--by", choices=["sequence", "coords"], default="sequence", help="Duplicate detection method")
    p.add_argument("--outdir", type=Path, default=None, help="Output directory (default: alongside input)")
    p.add_argument("--yes", action="store_true", help="Do not prompt; proceed with removals (good for batch)")
    p.add_argument("--dry-run", action="store_true", help="Scan and report, but do not write output files")

    # Mask options
    p.add_argument("--mask-res", type=str, default=None,
                   help="Comma-separated residue names to prefer when choosing which duplicate chain to keep (e.g., 'HEM,CU1'). Case-insensitive. Waters ignored.")
    p.add_argument("--mask-mode", choices=["all", "any"], default="all",
                   help="'all' = chain must contain all mask residue names; 'any' = at least one.")

    args = p.parse_args(argv)

    mask_res = parse_mask_resnames(args.mask_res)

    files = gather_files(args.path)
    if not files:
        print("[WARN] No matching files found.")
        return 1

    total_removed = 0
    wrote = 0

    for f in files:
        try:
            out_path, removed = dedupe_file(
                f, by=args.by, outdir=args.outdir, auto_yes=args.yes, dry_run=args.dry_run,
                mask_res=mask_res, mask_mode=args.mask_mode
            )
            total_removed += removed
            if not args.dry_run:
                wrote += 1
        except Exception as e:
            print(f"[ERROR] {f}: {e}")

    print(f"[DONE] Processed {len(files)} file(s); chains removed: {total_removed}; outputs written: {wrote if not args.dry_run else 0}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
