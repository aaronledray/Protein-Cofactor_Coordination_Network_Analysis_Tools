"""Conserved-interaction log from a set of structures (``sscna conserve``; experimental).

A thin layer over :mod:`modules.conserved_interactions` and :mod:`modules.design_log`: it reads a manifest and an optional
position map, runs the sequence vote and writes the log folder. It adds no analysis of its own.

Manifest CSV (one row per structure): ``structure_id``, ``path`` (relative paths are taken from the manifest's folder),
``sequence_id`` (structures of one sequence get one vote); optional ``group`` (default ``all``) and ``chain`` (default ``A``).
Position-map CSV (long format): ``structure_id``, ``residue_number``, optional ``insertion_code``, ``position`` (the shared
label, for example an alignment column). Residues of a structure without a row are left out, as in :func:`load_atoms`. With no
map the residue number is the position, which is right only when the numbering is already common to the set.

With ``--single-protein`` the manifest must name one ``sequence_id`` (an optional ``state`` column labels the structures) and the log reports,
per contact, in how many structures of that protein it is seen and how close the residues come, instead of the vote (see :mod:`modules.single_protein`).
With ``--stratify-state`` the default sequence vote is retained and an additional ``state_contacts`` table applies it independently within each
state in the manifest; this is descriptive and does not test for a state effect.

Defaults are the manuscript's rules (90% support, five voting sequences per group, 4.0 A contacts). The sulfur acceptor, neighbor
contacts, coordination tiers, metal ligands, directives and audits are opt-in. The vote holds every structure's atoms in memory
(the 54-sequence, 249-structure oxidase set peaked near 2.5 GB); it runs in one process.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pandas as pd

from . import conserved_interactions as ci
from . import design_log
from . import single_protein
from . import state_stratification

REQUIRED_MANIFEST = ("structure_id", "path", "sequence_id")
COFACTOR_FLAGS = ("coord", "polar", "ionic", "contact")
TIER_FLAGS = ("bond", "weak", "ligated")


def read_manifest(path) -> pd.DataFrame:
    """The manifest with ``structure_id, path, sequence_id, group, chain``; raises ``ValueError`` on a problem in it."""
    path = Path(path)
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_MANIFEST if c not in table.columns]
    if missing:
        raise ValueError(f"manifest lacks column(s) {missing}")
    for column, default in (("group", "all"), ("chain", "A")):
        if column not in table.columns:
            table[column] = default
        table[column] = table[column].replace("", default)
    blank = table[list(REQUIRED_MANIFEST)].eq("").any(axis=1)
    if blank.any():
        raise ValueError(f"manifest row(s) {(table.index[blank] + 2).tolist()} have an empty structure_id, path or sequence_id")
    duplicated = table.loc[table["structure_id"].duplicated(), "structure_id"].unique().tolist()
    if duplicated:
        raise ValueError(f"duplicate structure_id in the manifest: {duplicated}")
    table["path"] = [str(p if (p := Path(v).expanduser()).is_absolute() else path.parent / p) for v in table["path"]]
    absent = [p for p in table["path"] if not Path(p).is_file()]
    if absent:
        raise ValueError(f"{len(absent)} structure file(s) not found, first: {absent[0]}")
    return table


def read_position_map(path) -> Dict[str, Dict[tuple, object]]:
    """``structure_id`` -> ``{(residue_number, insertion_code): position}`` from the long-format CSV."""
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [c for c in ("structure_id", "residue_number", "position") if c not in table.columns]
    if missing:
        raise ValueError(f"position map lacks column(s) {missing}")
    if "insertion_code" not in table.columns:
        table["insertion_code"] = ""
    keyed = table.assign(residue_number=pd.to_numeric(table["residue_number"], errors="raise").astype(int))
    if keyed.duplicated(["structure_id", "residue_number", "insertion_code"]).any():
        raise ValueError("position map has a residue listed twice for one structure")
    clash = keyed.groupby("structure_id")["position"].apply(lambda s: s.duplicated().any())
    if clash.any():
        raise ValueError(f"position map assigns one position to two residues in {clash[clash].index.tolist()[:5]}")
    numeric = pd.to_numeric(keyed["position"], errors="coerce")
    keyed["position"] = numeric.astype(int) if numeric.notna().all() else keyed["position"]
    return {sid: {(int(r.residue_number), r.insertion_code.strip()): r.position for r in frame.itertuples()} for sid, frame in keyed.groupby("structure_id")}


def split_names(value: Optional[str]) -> List[str]:
    return [item.strip() for item in (value or "").split(",") if item.strip()]


def load_set(manifest: pd.DataFrame, position_map: Optional[Dict[str, Dict[tuple, object]]], *, metals: bool):
    """(atoms, metals table) of every structure in the manifest."""
    atoms, found = [], []
    for row in manifest.itertuples():
        mapping = None
        if position_map is not None:
            if row.structure_id not in position_map:
                raise ValueError(f"the position map has no rows for structure {row.structure_id}")
            mapping = position_map[row.structure_id]
        table = ci.load_atoms(row.path, row.chain, mapping, structure_id=row.structure_id)
        if table.empty:
            raise ValueError(f"no residues of chain {row.chain} in {row.path} are mapped to a position")
        atoms.append(table)
        if metals:
            found.append(ci.load_metals(row.path, structure_id=row.structure_id))
    return pd.concat(atoms, ignore_index=True), (pd.concat(found, ignore_index=True) if found else None)


def metal_features(metals: pd.DataFrame, elements: Sequence[str]) -> pd.DataFrame:
    """One ``feature`` per metal element (all atoms of the element in a structure), ``kind`` ``metal``."""
    wanted = [e.upper() for e in elements]
    chosen = metals[metals["element"].isin(wanted)]
    return pd.DataFrame({"structure_id": chosen["structure_id"], "feature": chosen["element"], "atom": chosen["element"], "element": chosen["element"],
                         "kind": "metal", "x": chosen["x"], "y": chosen["y"], "z": chosen["z"]}).reset_index(drop=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="sscna conserve", description="Write a conserved-interaction log for a set of structures (experimental).")
    parser.add_argument("--manifest", required=True, help="CSV with structure_id, path, sequence_id and optionally group, chain")
    parser.add_argument("--position-map", help="CSV with structure_id, residue_number, [insertion_code], position; default: the residue number")
    parser.add_argument("--output-dir", required=True, help="log folder to write")
    parser.add_argument("--overwrite", action="store_true", help="replace an existing log in --output-dir")
    parser.add_argument("--title", default=None, help="title of the log (default: a title for the mode)")
    parser.add_argument("--support", type=float, default=None, help="share of voting sequences that must show a contact (default 0.90)")
    parser.add_argument("--min-sequences", type=int, default=None, help="voting sequences a group needs to count (default 5)")
    parser.add_argument("--cutoff", type=float, default=4.0, help="residue contact distance in Angstrom (default 4.0)")
    parser.add_argument("--metals", help="comma-separated metal elements to treat as cofactor features, e.g. CU,FE (adds cofactor_contacts)")
    parser.add_argument("--coordination-tiers", action="store_true", help="opt-in: bond/weak tiers per metal and donor atom for the cofactor contacts (needs --metals)")
    parser.add_argument("--metal-ligands", action="store_true", help="opt-in: add the descriptive table of non-protein atoms bound to metals")
    parser.add_argument("--sulfur-acceptor", action="store_true", help="opt-in: sulfur hydrogen-bond acceptor flags")
    parser.add_argument("--neighbor-contacts", action="store_true", help="opt-in: polar and ionic contacts between sequence neighbors")
    parser.add_argument("--anchor", help="position of an anchor residue to audit against the first --metals element (adds audit_anchor)")
    parser.add_argument("--directives", action="store_true", help="opt-in: write the directives (directives.csv; descriptive, or observed with --single-protein)")
    parser.add_argument("--single-protein", action="store_true", help="the structures are of one protein (one sequence_id): report how many show each contact and the distances, not the vote")
    parser.add_argument("--stratify-state", action="store_true", help="opt-in: add residue-contact votes calculated independently within each manifest state")
    parser.add_argument("--near-margin", type=float, default=single_protein.NEAR_MARGIN, help="--single-protein: flag a contact whose median closest distance is within this many A of the cutoff (default 0.3)")
    return parser


def _stack(pieces: List[pd.DataFrame]) -> Optional[pd.DataFrame]:
    """Directives of several kinds as one table; integral positions stay integers where another kind leaves them empty."""
    return pd.concat(pieces, ignore_index=True).convert_dtypes() if pieces else None


def run_single(args: argparse.Namespace, manifest: pd.DataFrame, atoms: pd.DataFrame, metals: Optional[pd.DataFrame], elements: List[str]) -> Path:
    """The single-protein log: persistence across the structures of one sequence, no vote."""
    structures = manifest[["structure_id", "group", "sequence_id"] + (["state"] if "state" in manifest.columns else [])]
    notes = []
    if structures["structure_id"].nunique() == 1:
        notes.append("one structure: every contact is observed in one structure and nothing can be said about persistence")
    tables = {"contact_persistence": single_protein.contact_persistence(atoms, structures, cutoff=args.cutoff, near_margin=args.near_margin),
              "audit_backbone_only": ci.backbone_only_audit(atoms)}
    directives = [design_log.observed_contact_directives(tables["contact_persistence"], structures)] if args.directives else []
    if elements:
        features = metal_features(metals, elements)
        if features.empty:
            raise ValueError(f"no {', '.join(elements)} atoms found in the structures")
        model = "by_pair" if args.coordination_tiers else None
        flags = list(COFACTOR_FLAGS) + (list(TIER_FLAGS) if model else [])
        cofactor = ci.cofactor_contacts(atoms, features, coordination_model=model)
        if model:
            cofactor["bond"] = cofactor["coord_tier"] == "bond"
            cofactor["weak"] = cofactor["coord_tier"] == "weak"
            cofactor["ligated"] = cofactor["bond"] | cofactor["weak"]
        tables["cofactor_persistence"] = single_protein.cofactor_persistence(cofactor, structures, flags)
        if args.directives:
            directives.append(design_log.observed_cofactor_directives(tables["cofactor_persistence"], structures, flags))
    if args.metal_ligands:
        tables["metal_ligands"] = pd.concat([ci.load_metal_ligands(row.path, structure_id=row.structure_id) for row in manifest.itertuples()], ignore_index=True)
    parameters = {"mode": "single_protein", "cutoff": args.cutoff, "near_margin": args.near_margin, "metals": elements, "coordination_tiers": args.coordination_tiers,
                  "metal_ligands": args.metal_ligands, "directives": args.directives, "manifest": str(Path(args.manifest).name),
                  "position_map": "given" if args.position_map else "residue numbers"}
    pieces = [d for d in directives if len(d)]
    if args.directives and not pieces:
        notes.append("directives were requested and no contact was observed, so directives.csv is not written")
    return design_log.write_log(args.output_dir, tables, structures=structures, title=args.title or "Single-protein contact log", parameters=parameters, notes=notes, overwrite=args.overwrite,
                                directives=_stack(pieces), mode="single_protein")


def run(args: argparse.Namespace) -> Path:
    manifest = read_manifest(args.manifest)
    if args.stratify_state:
        if args.single_protein:
            raise ValueError("--stratify-state applies to the sequence vote, not --single-protein")
        if args.directives:
            raise ValueError("--directives currently describe the pooled vote and cannot be combined with --stratify-state")
        if "state" not in manifest.columns:
            raise ValueError("--stratify-state needs a non-empty state column in the manifest")
        if manifest["state"].eq("").any():
            raise ValueError("--stratify-state needs a non-empty state value on every manifest row")
        if len(single_protein.states_of(manifest)) < 2:
            raise ValueError("--stratify-state needs at least two distinct state labels")
    if args.single_protein:
        given = [flag for flag, value in (("--support", args.support), ("--min-sequences", args.min_sequences)) if value is not None]
        given += [flag for flag, value in (("--sulfur-acceptor", args.sulfur_acceptor), ("--neighbor-contacts", args.neighbor_contacts), ("--anchor", args.anchor)) if value]
        if given:
            raise ValueError(f"{', '.join(given)} do not apply with --single-protein (they belong to the sequence vote or the alignment)")
        manifest["group"] = "all"
        single_protein.require_one_sequence(manifest)
        single_protein.states_of(manifest)
    args.support = 0.90 if args.support is None else args.support
    args.min_sequences = 5 if args.min_sequences is None else args.min_sequences
    position_map = read_position_map(args.position_map) if args.position_map else None
    elements = split_names(args.metals)
    if args.coordination_tiers and not elements:
        raise ValueError("--coordination-tiers needs --metals")
    if args.anchor and not elements:
        raise ValueError("--anchor needs --metals")
    need_metals = bool(elements or args.sulfur_acceptor or args.neighbor_contacts)
    atoms, metals = load_set(manifest, position_map, metals=need_metals)
    if args.single_protein:
        return run_single(args, manifest, atoms, metals, elements)
    structures = manifest[["structure_id", "group", "sequence_id"] + (["state"] if args.stratify_state else [])]
    depth = design_log.group_depth(structures, support=args.support, min_sequences=args.min_sequences)
    small = depth.loc[~depth["meets_minimum"], "group"].tolist()
    notes = [f"group {g} has fewer than {args.min_sequences} sequences and cannot make a contact conserved" for g in small]
    if args.stratify_state:
        notes.append("state_contacts applies the sequence vote independently within each caller-supplied state; state results are descriptive and do not test for a state effect")

    contacts = ci.conserved_contacts(atoms, structures, cutoff=args.cutoff, support=args.support, min_sequences=args.min_sequences,
                                     sulfur_acceptor=args.sulfur_acceptor, neighbor_contacts=args.neighbor_contacts,
                                     metals=metals if args.sulfur_acceptor else None, coverage=True)
    tables = {"conserved_contacts": contacts, "audit_backbone_only": ci.backbone_only_audit(atoms)}
    if args.stratify_state:
        state_structures = structures
        tables["state_structure_votes"] = state_stratification.state_structure_votes(
            atoms, state_structures, cutoff=args.cutoff, sulfur_acceptor=args.sulfur_acceptor,
            metals=metals if args.sulfur_acceptor else None, neighbor_contacts=args.neighbor_contacts)
        tables["state_contacts"] = state_stratification.state_stratified_contacts(
            atoms, state_structures, cutoff=args.cutoff, support=args.support, min_sequences=args.min_sequences,
            sulfur_acceptor=args.sulfur_acceptor, metals=metals if args.sulfur_acceptor else None,
            neighbor_contacts=args.neighbor_contacts)
        tables["state_contrasts"] = state_stratification.state_contrasts(
            tables["state_contacts"], min_sequences=args.min_sequences)
    directives = [design_log.contact_directives(contacts, structures, support=args.support, min_sequences=args.min_sequences)] if args.directives else []

    if elements:
        features = metal_features(metals, elements)
        if features.empty:
            raise ValueError(f"no {', '.join(elements)} atoms found in the structures")
        model = "by_pair" if args.coordination_tiers else None
        flags = list(COFACTOR_FLAGS) + (list(TIER_FLAGS) if model else [])
        cofactor = ci.cofactor_contacts(atoms, features, coordination_model=model)
        if model:                                       # the tiers are voted as flags; `ligated` (bond or weak) as in the plastocyanin check
            cofactor["bond"] = cofactor["coord_tier"] == "bond"
            cofactor["weak"] = cofactor["coord_tier"] == "weak"
            cofactor["ligated"] = cofactor["bond"] | cofactor["weak"]
        vote = ci.sequence_vote(cofactor, structures, ["feature", "position"], flags, support=args.support, min_sequences=args.min_sequences)
        tables["cofactor_contacts"] = vote
        if args.directives:
            directives.append(design_log.cofactor_directives(vote, structures, flags, support=args.support, min_sequences=args.min_sequences))
        if args.anchor:
            anchor = int(args.anchor) if args.anchor.lstrip("-").isdigit() else args.anchor
            sites = features[features["feature"] == elements[0].upper()]
            tables["audit_anchor"] = ci.anchor_audit(atoms, sites[["structure_id", "x", "y", "z"]], anchor, structures)
    if args.metal_ligands:
        ligands = [ci.load_metal_ligands(row.path, structure_id=row.structure_id) for row in manifest.itertuples()]
        tables["metal_ligands"] = pd.concat(ligands, ignore_index=True)

    parameters = {"support": args.support, "min_sequences": args.min_sequences, "cutoff": args.cutoff, "metals": elements,
                  "coordination_tiers": args.coordination_tiers, "sulfur_acceptor": args.sulfur_acceptor,
                  "neighbor_contacts": args.neighbor_contacts, "metal_ligands": args.metal_ligands, "directives": args.directives,
                  "position_map": "given" if position_map is not None else "residue numbers", "manifest": str(Path(args.manifest).name),
                  "stratify_state": args.stratify_state, **({"states": single_protein.states_of(manifest)} if args.stratify_state else {})}
    pieces = [d for d in directives if len(d)]
    if args.directives and not pieces:
        notes.append("directives were requested and nothing is conserved in any group, so directives.csv is not written")
    return design_log.write_log(args.output_dir, tables, structures=structures, title=args.title or "Conserved-interaction log", parameters=parameters, notes=notes,
                                overwrite=args.overwrite, directives=_stack(pieces))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(list(sys.argv[1:] if argv is None else argv))
    try:
        out = run(args)
    except (ValueError, OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as error:      # OSError: a missing file, FileExistsError
        print(f"sscna conserve: {error}", file=sys.stderr)
        return 2
    print(f"Wrote the log to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
