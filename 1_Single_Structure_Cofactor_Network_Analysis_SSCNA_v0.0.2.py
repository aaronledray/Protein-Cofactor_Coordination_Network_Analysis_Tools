#!/usr/bin/env python3
# -*- coding: utf-8 -*-






# - [ ] Well, this next run will be for just OEC..... or OEX....







"""
SSCNA — Single Structure Cofactor Network Analysis

Usage (non-interactive):
    python scripts/sscna.py \
        --template 3u7q_monomer.pdb \
        --cofactor ICS,CLF \
        --cofactor2 HCA \
        --distance 3.6 \
        --combinatorial \
        --combinatorial-cutoff 20.0 \
        --exclude-moieties alanine_sidechain \
        --mode Coord_Network

Or:
    python scripts/sscna.py --template 3u7q_monomer.pdb
and it will discover config / prompt for the rest.



Config discovery (optional; overrides defaults, can be overridden by CLI):
- sscna.yaml or sscna.json located in the *same directory* as the template PDB.

Outputs:
- "<template>_Coord_Breakdown.csv"
- "1_static_*" PNG images
- "1_template_coordination_network.html"



"""






import os
import sys
import json
import argparse
from typing import List, Dict, Any, Optional

# Optional YAML support
try:
    import yaml  # type: ignore
    HAS_YAML = True
except Exception:
    HAS_YAML = False

# ---- Project modules (relative import path assumes project root execution) ----
# If you run from repository root: `python scripts/sscna.py ...`
# adjust sys.path so `modules` is importable.
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from modules.io_utils import unpack_pdb_file


from modules.analysis import (
    identify_coordination_network,
    identify_residues_of_interest,
    generate_coordination_csv_with_moieties,
)


from modules.plotting import (
    static_plots_2d,
    plot_interactive_modes_with_network,
    plot_interactive_modes_with_roi,
    plot_template_heatmap_interactive,  # available for later
)



from modules.reporting import write_coord_breakdown_v2


from modules.moieties import bond_lookup, atom_type_colors, chemical_moieties





# structure_utils provides helpers used by analysis/plotting flows
# from modules.structure_utils import (
#     # exposed helpers if you need them directly later
#)

# --------------------------
# Defaults (sane fallbacks)
# --------------------------
DEFAULTS = {
    "mode": "Coord_Network",                  # or "Residues_of_Interest"
    "distance_cutoff": 3.6,                  # moiety interaction cutoff
    "expand_residues": False,
    "combinatorial": False,
    "combinatorial_cofactor_cutoff": 20.0,
    "exclude_moieties": ["alanine_sidechain"],

    # Residues-of-Interest mode (template-relative indices before +1 step)
    "residues_of_interest_integers": [29, 32, 39, 42, 43, 45, 46,
                                      64, 67, 68, 71, 72, 89, 93,
                                      97, 99, 103, 104, 107, 138, 142],
}








import os
from typing import Optional, Tuple

def resolve_structure_path(user_path: str) -> Optional[str]:
    """
    Return a real path to the structure file if it exists.
    Tries:
      - Exact path (after expanding ~ and making absolute)
      - Same basename with common structure extensions swapped in
      - .gz variants
    """
    if not user_path:
        return None

    p = os.path.abspath(os.path.expanduser(user_path))
    if os.path.isfile(p):
        return p

    base, ext = os.path.splitext(p)

    # Handle double extensions like .cif.gz
    if ext == ".gz":
        base2, ext2 = os.path.splitext(base)
        primary_ext = ext2 + ext   # e.g., ".cif.gz"
        base = base2
    else:
        primary_ext = ext

    # Try alternative extensions
    candidates = []
    # Preserve original first (already checked)
    # Common structure extensions (and compressed)
    exts = [".pdb", ".cif", ".mmcif"]
    for e in exts:
        candidates.append(base + e)
        candidates.append(base + e + ".gz")

    # If the user provided a bare PDB ID like "1ag6", try adding extensions
    if os.path.sep not in user_path and len(os.path.basename(base)) in (4, 5):
        for e in exts:
            candidates.append(os.path.join(os.getcwd(), os.path.basename(base) + e))
            candidates.append(os.path.join(os.getcwd(), os.path.basename(base) + e + ".gz"))

    for c in candidates:
        if os.path.isfile(c):
            return c

    return None


def load_structure_generic(struct_path: str):
    """
    Load either PDB or mmCIF with Biopython and return a Structure object.
    """
    from Bio.PDB import PDBParser, MMCIFParser
    sp = struct_path.lower()
    if sp.endswith(".pdb") or sp.endswith(".pdb.gz"):
        parser = PDBParser(QUIET=True)
    elif sp.endswith(".cif") or sp.endswith(".mmcif") or sp.endswith(".cif.gz") or sp.endswith(".mmcif.gz"):
        parser = MMCIFParser(QUIET=True)
    else:
        # Fallback: try CIF first, then PDB
        try:
            parser = MMCIFParser(QUIET=True)
            return parser.get_structure("structure", struct_path)
        except Exception:
            parser = PDBParser(QUIET=True)

    return parser.get_structure("structure", struct_path)













# --------------------------
# Small utilities
# --------------------------
def str_to_list(s: Optional[str]) -> List[str]:
    if not s:
        return []
    return [x.strip() for x in s.split(",") if x.strip()]

def load_sidecar_config(template_path: str) -> Dict[str, Any]:
    """
    Look for sscna.yaml or sscna.json in the directory of the template PDB.
    Return {} if none found.
    """
    base_dir = os.path.dirname(os.path.abspath(template_path))
    yaml_path = os.path.join(base_dir, "sscna.yaml")
    json_path = os.path.join(base_dir, "sscna.json")

    if HAS_YAML and os.path.isfile(yaml_path):
        try:
            with open(yaml_path, "r") as f:
                return yaml.safe_load(f) or {}
        except Exception as e:
            print(f"[WARN] Could not parse {yaml_path}: {e}")
    if os.path.isfile(json_path):
        try:
            with open(json_path, "r") as f:
                return json.load(f) or {}
        except Exception as e:
            print(f"[WARN] Could not parse {json_path}: {e}")
    return {}

def prompt_interactive(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Prompt the user for any missing critical fields.
    """
    print("\n[SSCNA] Interactive setup (press Enter to accept defaults)")

    # Template PDB (required)
    while not params.get("template_pdb_file"):
        tpl = input("Template PDB path: ").strip()
        if tpl:
            params["template_pdb_file"] = tpl

    # Cofactor 1
    if not params.get("cofactor_resname"):
        cf1 = input("Cofactor residue name(s) (comma-separated, e.g., ICS,CLF): ").strip()
        params["cofactor_resname"] = str_to_list(cf1) if cf1 else []

    # Optional cofactor2 list
    if params.get("cofactor_resname2") is None:
        cf2 = input("Second cofactor residue(s) (comma-separated) [optional]: ").strip()
        params["cofactor_resname2"] = str_to_list(cf2) if cf2 else None

    # Mode
    mode = params.get("mode") or DEFAULTS["mode"]
    mode_in = input(f"Mode [Coord_Network / Residues_of_Interest] (default {mode}): ").strip()
    if mode_in:
        params["mode"] = mode_in

    # Numbers
    def ask_float(key, label, default):
        cur = params.get(key, default)
        raw = input(f"{label} (default {cur}): ").strip()
        if raw:
            try:
                params[key] = float(raw)
            except ValueError:
                print(f"[WARN] Invalid number for {label}; keeping {cur}")

    ask_float("distance_cutoff", "Distance cutoff (Å)", DEFAULTS["distance_cutoff"])
    ask_float("combinatorial_cofactor_cutoff", "Combinatorial cofactor cutoff (Å)", DEFAULTS["combinatorial_cofactor_cutoff"])

    # Booleans
    def ask_bool(key, label, default):
        cur = params.get(key, default)
        raw = input(f"{label} [y/n] (default {'y' if cur else 'n'}): ").strip().lower()
        if raw in ("y", "yes"):
            params[key] = True
        elif raw in ("n", "no"):
            params[key] = False

    ask_bool("expand_residues", "Expand PCS/SCS to include entire residues", DEFAULTS["expand_residues"])
    ask_bool("combinatorial", "Use combinatorial mode", DEFAULTS["combinatorial"])

    # Exclude moieties
    if params.get("exclude_moieties") is None:
        em = input(f"Exclude moieties (comma-separated) [default {DEFAULTS['exclude_moieties']}]: ").strip()
        params["exclude_moieties"] = str_to_list(em) if em else DEFAULTS["exclude_moieties"]

    # ROI list if selecting Residues_of_Interest
    if params.get("mode") == "Residues_of_Interest" and not params.get("residues_of_interest_integers"):
        raw = input("Residues-of-interest (template-relative indices, comma-separated) "
                    f"[default: {DEFAULTS['residues_of_interest_integers']}]: ").strip()
        if raw:
            try:
                params["residues_of_interest_integers"] = [int(x) for x in raw.split(",")]
            except Exception:
                print("[WARN] Invalid ROI list; using defaults")
                params["residues_of_interest_integers"] = DEFAULTS["residues_of_interest_integers"]
        else:
            params["residues_of_interest_integers"] = DEFAULTS["residues_of_interest_integers"]

    print("")
    return params

def merge(a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
    """b overrides a; shallow merge is sufficient for our flat params."""
    out = a.copy()
    out.update({k: v for k, v in b.items() if v is not None})
    return out

# --------------------------
# Main run
# --------------------------
def run_coord_network(params: Dict[str, Any]) -> None:
    tpl = params["template_pdb_file"]
    tpl_basename = os.path.basename(tpl)
    cofactor1 = params["cofactor_resname"]
    cofactor2 = params.get("cofactor_resname2")
    distance_cutoff = params["distance_cutoff"]
    expand_residues = params["expand_residues"]
    combinatorial = params["combinatorial"]
    comb_cutoff = params["combinatorial_cofactor_cutoff"]
    exclude_moieties = params["exclude_moieties"] or []
    output_dir = os.path.join(os.getcwd(), "SSCNA_output")
    os.makedirs(output_dir, exist_ok=True)
    file_prefix = f"{tpl_basename}_"

    if isinstance(cofactor1, str):
        cofactor1 = [cofactor1]
    if cofactor2 and isinstance(cofactor2, str):
        cofactor2 = [cofactor2]

    print(f"[INFO] Loading structure: {tpl}")
    structure, atoms_data = unpack_pdb_file(tpl)

    print(f"[INFO] Identifying coordination network for {cofactor1}"
          + (f" with cofactor2 {cofactor2}" if cofactor2 else ""))







    cofactor_sphere, pcs_atoms, scs_atoms = identify_coordination_network(
        structure=structure,
        cofactor_resname=cofactor1,
        distance_cutoff=distance_cutoff,
        expand_residues=expand_residues,
        combinatorial_mode=combinatorial,
        combinatorial_cofactor_cutoff=comb_cutoff,
        cofactor_resname2=cofactor2,
        exclude_moieties=exclude_moieties,
        output_dir=output_dir,
        output_prefix=file_prefix,
    )









    # Collect coordinate lists for static plotting
    cofactor_template_coords = [a['coordinates'] for a in cofactor_sphere]
    pcs_template_coords      = [a['coordinates'] for a in pcs_atoms]
    scs_template_coords      = [a['coordinates'] for a in scs_atoms]

    # CSV report
    print("[INFO] Writing CSV breakdown...")
    generate_coordination_csv_with_moieties(
        cofactor_sphere=cofactor_sphere,
        pcs_residues=pcs_atoms,
        scs_residues=scs_atoms,
        pdb_file_name=tpl_basename,
        structure=structure,
        bond_lookup=bond_lookup,
        output_dir=output_dir,
    )

    write_coord_breakdown_v2(
        structure=structure,
        cofactor_atoms=cofactor_sphere,
        pcs_atoms=pcs_atoms,
        scs_atoms=scs_atoms,
        moiety_lookup=chemical_moieties,
        output_csv_path=os.path.join(output_dir, f"{file_prefix}Coord_Breakdown_atoms.csv"),
    )

    # Static plots
    print("[INFO] Rendering static plots...")
    static_plots_2d(
        cofactor_coords=cofactor_template_coords,
        pcs_coords=pcs_template_coords,
        scs_coords=scs_template_coords,
        cofactor_atoms=cofactor_sphere,
        pcs_atoms=pcs_atoms,
        scs_atoms=scs_atoms,
        structure=structure,
        bond_lookup=bond_lookup,
        pdb_name=tpl_basename,
        cofactor_resname=cofactor1,
        include_bonds=True,
        focused_bonds=True,
        output_prefix=os.path.join(output_dir, f"{file_prefix}1_")
    )

    # Interactive plot
    print("[INFO] Rendering interactive Plotly graph...")
    plot_interactive_modes_with_network(
        structure=structure,
        cofactor_atoms=cofactor_sphere,
        pcs_atoms=pcs_atoms,
        scs_atoms=scs_atoms,
        bond_lookup_table=bond_lookup,
        pdb_name=tpl_basename,
        cofactor_resname=cofactor1,
        atom_type_colors=atom_type_colors,
        output_filename=os.path.join(output_dir, f"{file_prefix}1_template_coordination_network.html"),
        links_csv_path=os.path.join(output_dir, f"{file_prefix}Coord_Links.csv"),
    )

    print("[DONE] Coord_Network complete.")

def run_roi(params: Dict[str, Any]) -> None:
    tpl = params["template_pdb_file"]
    tpl_basename = os.path.basename(tpl)
    distance_cutoff = 10.0  # for ROI pipeline (as in your notes)
    cofactor_for_roi = params.get("roi_cofactor_resname") or ["HM1"]
    alternative = None
    combinatorial = False
    comb_cutoff = 10.0
    cofactor2 = None
    output_dir = os.path.join(os.getcwd(), "SSCNA_output")
    os.makedirs(output_dir, exist_ok=True)

    roi_integers = params.get("residues_of_interest_integers") or DEFAULTS["residues_of_interest_integers"]

    print(f"[INFO] Loading structure: {tpl}")
    structure, atoms_data = unpack_pdb_file(tpl)

    print(f"[INFO] Identifying ROI with cofactor {cofactor_for_roi}...")
    cofactor_atoms, roi_atoms = identify_residues_of_interest(
        structure=structure,
        cofactor_resname=cofactor_for_roi,
        alternative_resname=alternative,
        combinatorial_mode=combinatorial,
        combinatorial_cofactor_cutoff=comb_cutoff,
        cofactor_resname2=cofactor2,
        residues_of_interest_integers=roi_integers
    )

    print(f"[INFO] Cofactor atoms: {len(cofactor_atoms)} | ROI atoms: {len(roi_atoms)}")

    print("[INFO] Rendering interactive Plotly (ROI)...")
    plot_interactive_modes_with_roi(
        structure=structure,
        cofactor_atoms=cofactor_atoms,
        roi_atoms=roi_atoms,
        bond_lookup_table=bond_lookup,
        pdb_name=tpl_basename,
        cofactor_resname=cofactor_for_roi,
        atom_type_colors=atom_type_colors,
        output_filename=os.path.join(output_dir, f"{tpl_basename}_roi_coordination_network.html")
    )

    print("[DONE] Residues_of_Interest complete.")

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="SSCNA — Single Structure Cofactor Network Analysis")
    p.add_argument("--template", dest="template_pdb_file", help="Path to template PDB file")
    p.add_argument("--cofactor", dest="cofactor_resname", help="Cofactor residue name(s), comma-separated")
    p.add_argument("--cofactor2", dest="cofactor_resname2", help="Second cofactor residue name(s), comma-separated")
    p.add_argument("--distance", dest="distance_cutoff", type=float, help="Distance cutoff for moiety interactions (Å)")
    p.add_argument("--expand-residues", action="store_true", help="Expand PCS/SCS to entire residues")
    p.add_argument("--no-expand-residues", action="store_true", help="Do not expand residues")
    p.add_argument("--combinatorial", action="store_true", help="Enable combinatorial mode")
    p.add_argument("--no-combinatorial", action="store_true", help="Disable combinatorial mode")
    p.add_argument("--combinatorial-cutoff", dest="combinatorial_cofactor_cutoff", type=float, help="Cutoff for cofactor combinatorial filtering (Å)")
    p.add_argument("--exclude-moieties", help="Comma-separated moieties to exclude")
    p.add_argument("--mode", choices=["Coord_Network", "Residues_of_Interest"], help="Run mode")
    p.add_argument("--interactive", action="store_true", help="Force interactive prompting")
    return p.parse_args()

def main():
    args = parse_args()

    # Build params from CLI
    cli_params: Dict[str, Any] = {}
    if args.template_pdb_file:
        cli_params["template_pdb_file"] = args.template_pdb_file
    if args.cofactor_resname:
        cli_params["cofactor_resname"] = str_to_list(args.cofactor_resname)
    if args.cofactor_resname2:
        cli_params["cofactor_resname2"] = str_to_list(args.cofactor_resname2)
    if args.distance_cutoff is not None:
        cli_params["distance_cutoff"] = args.distance_cutoff
    if args.combinatorial_cofactor_cutoff is not None:
        cli_params["combinatorial_cofactor_cutoff"] = args.combinatorial_cofactor_cutoff
    if args.exclude_moieties:
        cli_params["exclude_moieties"] = str_to_list(args.exclude_moieties)
    if args.mode:
        cli_params["mode"] = args.mode

    if args.expand_residues and not args.no_expand_residues:
        cli_params["expand_residues"] = True
    elif args.no_expand_residues:
        cli_params["expand_residues"] = False

    if args.combinatorial and not args.no_combinatorial:
        cli_params["combinatorial"] = True
    elif args.no_combinatorial:
        cli_params["combinatorial"] = False



    # Merge: defaults <- sidecar config <- CLI
    # If we have template already, we can try to load sidecar config from its folder.
    sidecar = {}
    if cli_params.get("template_pdb_file"):
        sidecar = load_sidecar_config(cli_params["template_pdb_file"])

    params = merge(DEFAULTS, sidecar)
    params = merge(params, cli_params)

    # If still missing critical info or --interactive, prompt
    need_prompt = args.interactive or not params.get("template_pdb_file") or not params.get("cofactor_resname")
    if need_prompt:
        params = prompt_interactive(params)

    # # Final sanity
    # if not os.path.isfile(params["template_pdb_file"]):
    #     sys.exit(f"[ERROR] Template PDB not found: {params['template_pdb_file']}")


    # Final sanity + resolve structure path (supports .pdb / .cif / .mmcif / .gz)
    resolved_path = resolve_structure_path(params.get("template_pdb_file"))
    if not resolved_path:
        sys.exit(f"[ERROR] Template structure file not found (tried PDB/CIF variants): {params.get('template_pdb_file')}")
    params["template_pdb_file"] = resolved_path  # normalize

    # Optional: Load once here and pass the Structure down
    try:
        params["structure_obj"] = load_structure_generic(resolved_path)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to parse structure '{resolved_path}': {e}")


    if not params.get("cofactor_resname"):
        sys.exit("[ERROR] At least one cofactor residue name is required.")

    # Run selected mode
    if params.get("mode") == "Residues_of_Interest":
        run_roi(params)
    else:
        run_coord_network(params)

if __name__ == "__main__":
    main()
