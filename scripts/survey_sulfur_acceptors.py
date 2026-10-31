#!/usr/bin/env python3
"""Survey N-H...S geometry around Cys SG and Met SD in atomic-resolution PDB structures.

Evidence for the sulfur-acceptor windows of ``residue_contacts(sulfur_acceptor=True)`` (TODO item C,
docs/scope_checks/03_sulfur_acceptor_survey.md). Same two steps as ``survey_cu_type1.py``:

    python scripts/survey_sulfur_acceptors.py fetch     --cache survey_cache/sulfur
    python scripts/survey_sulfur_acceptors.py summarize --cache survey_cache/sulfur --out docs/scope_checks

``fetch`` takes one representative per 30% sequence-identity cluster (RCSB grouping) of X-ray entries at or
below ``--max-resolution`` (1.2 A by default, about 2,100 structures) and downloads the mmCIF files (git-ignored).
For the first chain of each representative entity, ``summarize`` classifies every Cys SG (disulfide, metal-bound,
free) and Met SD (metal-bound, free) and lists the backbone-N and side-chain donor atoms within 4.3 A that are
at least two residues away in sequence. For a backbone N the amide hydrogen is placed from the backbone
(on the bisector of the C(i-1)-N and CA-N bonds, 1.01 A from N), and a contact is H-bond-like when
H...S <= 2.9 A and the N-H...S angle >= 120 degrees. Side-chain donors have no placed hydrogen and are reported by
distance only.
"""
import argparse
import json
import sys
import urllib.request
import warnings
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from survey_cu_type1 import FILES_URL, GRAPHQL_URL, SEARCH_URL, _post_json  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from modules import coordination_table  # noqa: E402
from modules.protonation import (AMIDE_DONORS, NH_LENGTH, RING_DONORS, amide_nh2_hydrogens, backbone_amide_hydrogen,  # noqa: E402
                                 hydrogen_bond_geometry, ring_nh_hydrogen, SULFUR_ANGLE_MIN, SULFUR_H_MAX)

amide_hydrogen = backbone_amide_hydrogen
nh_geometry = hydrogen_bond_geometry

METALS = {"LI", "NA", "MG", "K", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "MO", "W", "V", "CD", "HG", "AG", "AU", "PT", "PD",
          "RU", "IR", "RH", "OS", "RE", "TI", "CR", "FE2"}
SC_DONORS = {("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("ASN", "ND2"), ("GLN", "NE2"), ("ARG", "NE"), ("ARG", "NH1"),
             ("ARG", "NH2"), ("LYS", "NZ"), ("TRP", "NE1"), ("HIS", "ND1"), ("HIS", "NE2")}
GROUPS = {("TRP", "NE1"): "Trp NE1", ("ASN", "ND2"): "Asn ND2", ("GLN", "NE2"): "Gln NE2", ("ARG", "NE"): "Arg N", ("ARG", "NH1"): "Arg N",
          ("ARG", "NH2"): "Arg N", ("HIS", "ND1"): "His N (as if protonated)", ("HIS", "NE2"): "His N (as if protonated)",
          ("SER", "OG"): "hydroxyl (H not placed)", ("THR", "OG1"): "hydroxyl (H not placed)", ("TYR", "OH"): "hydroxyl (H not placed)",
          ("LYS", "NZ"): "Lys NZ (H not placed)"}
DISULFIDE = 2.3
CYS_METAL = 2.7
MET_METAL = 3.3
NEAR = 4.3
HS_MAX = SULFUR_H_MAX
ANGLE_MIN = SULFUR_ANGLE_MIN
BINS = np.round(np.arange(3.0, 4.3, 0.1), 1)
KINDS = ("Cys metal-bound", "Cys disulfide", "Cys free", "Met metal-bound", "Met free")


def term(attribute, operator, value):
    return {"type": "terminal", "service": "text", "parameters": {"attribute": attribute, "operator": operator, "value": value}}


def search_entities(max_resolution=1.2, identity=30):
    """``ENTRY_entity`` ids, one representative per ``identity``% sequence cluster of X-ray entities."""
    query = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        term("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        term("rcsb_entry_info.resolution_combined", "less_or_equal", max_resolution)]},
        "return_type": "polymer_entity",
        "request_options": {"group_by": {"aggregation_method": "sequence_identity", "similarity_cutoff": identity},
                            "group_by_return_type": "representatives", "paginate": {"start": 0, "rows": 10000}}}
    return sorted(item["identifier"] for item in _post_json(SEARCH_URL, query, timeout=120)["result_set"])


def first_chains(entity_ids):
    """{entry: first author chain of the representative entity}."""
    wanted = {}
    for entity in entity_ids:
        entry, number = entity.split("_")
        wanted[entry] = number
    query = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id polymer_entities{
        rcsb_polymer_entity_container_identifiers{entity_id auth_asym_ids}}}}"""
    out = {}
    entries = sorted(wanted)
    for start in range(0, len(entries), 100):
        reply = _post_json(GRAPHQL_URL, {"query": query, "variables": {"ids": entries[start:start + 100]}})
        for entry in reply["data"]["entries"]:
            for entity in entry["polymer_entities"]:
                ids = entity["rcsb_polymer_entity_container_identifiers"]
                if ids["entity_id"] == wanted[entry["rcsb_id"]] and ids["auth_asym_ids"]:
                    out[entry["rcsb_id"]] = ids["auth_asym_ids"][0]
    return out


def fetch(cache, max_resolution=1.2, workers=8):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    chains = first_chains(search_entities(max_resolution))
    todo = [e for e in chains if not (cache / f"{e}.cif").exists()]

    def get(entry):
        with urllib.request.urlopen(FILES_URL.format(entry), timeout=120) as response:
            (cache / f"{entry}.cif").write_bytes(response.read())

    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(get, todo))
    (cache / "chains.json").write_text(json.dumps({"max_resolution": max_resolution, "chains": chains}))
    print(f"{len(chains)} representative structures at <= {max_resolution} A; downloaded {len(todo)} new into {cache}")


def planar_hydrogens(residue, name):
    """Hydrogen positions of a planar donor placed from its heavy atoms, or None when the group is not planar or atoms are missing."""
    key = (residue.resname, name)
    if key in RING_DONORS:
        a, b = RING_DONORS[key]
        if name not in residue or a not in residue or b not in residue:
            return None
        return [ring_nh_hydrogen(residue[name].coord, residue[a].coord, residue[b].coord)]
    if key in AMIDE_DONORS:
        carbon, reference = AMIDE_DONORS[key]
        if name not in residue or carbon not in residue or reference not in residue:
            return None
        return amide_nh2_hydrogens(residue[name].coord, residue[carbon].coord, residue[reference].coord)
    return None


def sulfur_contacts(cif_path, chain_id, min_separation=1):
    """Donor contacts of every Cys SG and Met SD in one chain (first model): one dict per donor atom within ``NEAR``.

    ``offset`` is the donor residue's sequence index minus the sulfur residue's; donors fewer than ``min_separation`` residues away
    are skipped (the report uses 2 for its main tables, as ``residue_contacts`` does, and reports the +-1 neighbors apart).
    """
    from Bio.PDB import MMCIFParser, NeighborSearch

    warnings.filterwarnings("ignore")
    model = MMCIFParser(QUIET=True).get_structure("x", str(cif_path))[0]
    if chain_id not in model:
        return []
    residues, order, previous, index = [], 0, None, {}
    for residue in model[chain_id]:
        if residue.id[0] != " ":
            continue
        number = int(residue.id[1])
        order += 1 if previous is None or number <= previous else number - previous
        previous = number
        residues.append(residue)
        index[residue.id] = order
    everything = [a for a in model.get_atoms() if a.get_parent().id[0] != "W"]
    metals = [a for a in everything if (a.element or "").upper() in METALS]
    sulfurs = [a for a in everything if a.get_id() == "SG" and a.get_parent().resname == "CYS"]
    search = NeighborSearch([a for r in residues for a in r])
    position = {r.id: n for n, r in enumerate(residues)}
    rows = []
    for s in (a for r in residues for a in r):
        residue = s.get_parent()
        if residue.resname == "CYS" and s.get_id() == "SG":
            if any(o is not s and np.linalg.norm(o.coord - s.coord) < DISULFIDE for o in sulfurs):
                kind = "Cys disulfide"
            elif any(np.linalg.norm(m.coord - s.coord) <= CYS_METAL for m in metals):
                kind = "Cys metal-bound"
            else:
                kind = "Cys free"
        elif residue.resname == "MET" and s.get_id() == "SD":
            kind = "Met metal-bound" if any(np.linalg.norm(m.coord - s.coord) <= MET_METAL for m in metals) else "Met free"
        else:
            continue
        for atom in search.search(s.coord, NEAR):
            other = atom.get_parent()
            if other is residue or abs(index[other.id] - index[residue.id]) < min_separation:
                continue
            name = atom.get_id()
            backbone = name == "N" and other.resname != "PRO"
            if not (backbone or (other.resname, name) in SC_DONORS):
                continue
            if not backbone and any(coordination_table.coordination_tier(m.element.upper(), other.resname, name, float(np.linalg.norm(m.coord - atom.coord))).tier
                                    in ("bond", "default") for m in metals if np.linalg.norm(m.coord - atom.coord) <= 2.8):
                continue                                                  # a metal coordinates this atom: deprotonated, not a donor
            row = {"pdb": Path(cif_path).stem[:4].upper(), "kind": kind, "donor": "backbone N" if backbone else "side chain",
                   "donor_atom": f"{other.resname} {name}", "sulfur": f"{residue.resname}{residue.id[1]}",
                   "donor_res": f"{other.resname}{other.id[1]}", "offset": index[other.id] - index[residue.id], "distance": float(np.linalg.norm(atom.coord - s.coord)),
                   "h_s": None, "angle": None}
            i = position[other.id]
            if backbone and i > 0 and "C" in residues[i - 1] and "CA" in other \
                    and np.linalg.norm(residues[i - 1]["C"].coord - atom.coord) < 1.6:
                row["h_s"], row["angle"] = nh_geometry(amide_hydrogen(atom.coord, residues[i - 1]["C"].coord, other["CA"].coord),
                                                       atom.coord, s.coord)
            if not backbone:
                row["group"] = GROUPS[(other.resname, name)]
                placed = planar_hydrogens(other, name)
                if placed:
                    best = min((nh_geometry(h, atom.coord, s.coord) for h in placed), key=lambda hg: hg[0])
                    row["h_s"], row["angle"] = best
            rows.append(row)
        rows.append({"pdb": Path(cif_path).stem[:4].upper(), "kind": kind, "donor": "sulfur", "donor_atom": "", "distance": 0.0,
                     "h_s": None, "angle": None})
    return rows


def _job(args):
    path, chain = args
    return sulfur_contacts(path, chain)


def hbond_like(frame):
    return (frame["h_s"] <= HS_MAX) & (frame["angle"] >= ANGLE_MIN)


def report(frame, n_structures, max_resolution):
    sulfur = frame[frame["donor"] == "sulfur"]
    everything = frame[frame["donor"] != "sulfur"]
    donors = everything[everything["offset"].abs() >= 2]
    backbone = donors[(donors["donor"] == "backbone N") & donors["h_s"].notna()].copy()
    backbone["like"] = hbond_like(backbone)
    lines = ["# N-H...S geometry around Cys SG and Met SD (PDB survey)", "",
             f"Generated by `scripts/survey_sulfur_acceptors.py`. One representative per 30% sequence-identity cluster of X-ray "
             f"entries at <= {max_resolution} A: {n_structures} structures (first chain of the representative entity). Donors are "
             f"backbone N (not Pro) and side-chain N/OH atoms at least two residues away in sequence, within {NEAR} A of the sulfur. "
             f"The amide H is placed from the backbone; a contact is H-bond-like when H...S <= {HS_MAX} A and the N-H...S angle "
             f">= {ANGLE_MIN:g} degrees.", "",
             "## Backbone N donors, by sulfur class", "",
             "| Sulfur | S atoms | backbone N within 4.0 A | H-bond-like | N...S of H-bond-like: min / median / max |", "|---|---|---|---|---|"]
    for kind in KINDS:
        sub = backbone[(backbone["kind"] == kind) & (backbone["distance"] <= 4.0)]
        like = sub[sub["like"]]
        spread = f"{like['distance'].min():.2f} / {like['distance'].median():.2f} / {like['distance'].max():.2f}" if len(like) else "none"
        lines.append(f"| {kind} | {int((sulfur['kind'] == kind).sum())} | {len(sub)} | {len(like)} | {spread} |")
    lines += ["", "## Fraction H-bond-like by N...S distance (backbone N donors)", "",
              "| Sulfur | " + " | ".join(f"{lo:.1f}-{lo + 0.1:.1f}" for lo in BINS if lo >= 3.2) + " |",
              "|---|" + "---|" * len([lo for lo in BINS if lo >= 3.2])]
    for kind in KINDS:
        sub = backbone[backbone["kind"] == kind]
        cells = []
        for lo in BINS:
            if lo < 3.2:
                continue
            bin_ = sub[(sub["distance"] >= lo) & (sub["distance"] < lo + 0.1)]
            cells.append(f"{int(bin_['like'].sum())}/{len(bin_)}" if len(bin_) else "-")
        lines.append(f"| {kind} | " + " | ".join(cells) + " |")
    lines += ["", "## Side-chain donors", "",
              "Planar groups get a placed hydrogen (Trp NE1, Asn ND2, Gln NE2, Arg N, and each His ring N as if it carried the H); "
              "hydroxyl and Lys NZ hydrogens rotate and are not placed, so those rows are counts by distance only. A donor within "
              "2.8 A of a metal (bond tier, else 2.8 A) is dropped. Cells are H-bond-like / contacts.", ""]
    side = donors[donors["donor"] == "side chain"].copy()
    side["like"] = hbond_like(side)
    classes = (("Cys metal-bound or free", ["Cys metal-bound", "Cys free"], 3.7), ("Cys disulfide", ["Cys disulfide"], 3.7), ("Met free", ["Met free"], 3.6))
    for group in ("Trp NE1", "Asn ND2", "Gln NE2", "Arg N", "His N (as if protonated)"):
        lines += [f"### {group}", "", "| Sulfur | within its window | <= 3.4 | 3.4-3.6 | 3.6-3.8 | 3.8-4.0 |", "|---|---|---|---|---|---|"]
        for label, kinds, window in classes:
            sub = side[(side["group"] == group) & side["kind"].isin(kinds) & side["h_s"].notna()]
            cells = [f"{int(x['like'].sum())}/{len(x)}" for x in (sub[sub["distance"] <= window],
                     *[sub[(sub["distance"] > lo) & (sub["distance"] <= hi)] for lo, hi in ((0, 3.4), (3.4, 3.6), (3.6, 3.8), (3.8, 4.0))])]
            lines.append(f"| {label} (window {window}) | " + " | ".join(cells) + " |")
        lines.append("")
    lines += ["### Not placed (distance only)", "", "| Donor | Sulfur | contacts <= 4.0 A | <= 3.4 | 3.4-3.6 | 3.6-3.8 | 3.8-4.0 |", "|---|---|---|---|---|---|---|"]
    for group in ("hydroxyl (H not placed)", "Lys NZ (H not placed)"):
        for label, kinds, _ in classes:
            sub = side[(side["group"] == group) & side["kind"].isin(kinds) & (side["distance"] <= 4.0)]
            counts = [int(((sub["distance"] > lo) & (sub["distance"] <= hi)).sum()) for lo, hi in ((0, 3.4), (3.4, 3.6), (3.6, 3.8), (3.8, 4.0))]
            lines.append(f"| {group} | {label} | {len(sub)} | " + " | ".join(map(str, counts)) + " |")
    histidines = side[side["group"] == "His N (as if protonated)"]
    lines += ["", "### Histidine tautomer sensitivity", "",
              "For each His within the window of a sulfur (non-disulfide Cys 3.7 A, Met 3.6 A), which ring N would make an H-bond-like "
              "contact if it carried the proton. Only one tautomer working means the protonation state decides whether the contact is real.", "",
              "| Sulfur | His contacts | ND1-H only | NE2-H only | either | neither |", "|---|---|---|---|---|---|"]
    for label, kinds, window in classes[:1] + classes[2:]:
        sub = histidines[histidines["kind"].isin(kinds) & (histidines["distance"] <= window) & histidines["h_s"].notna()]
        states = {"ND1-H only": 0, "NE2-H only": 0, "either": 0, "neither": 0}
        for _, grp in sub.groupby(["pdb", "sulfur", "donor_res"]):
            nd1 = bool(grp[grp["donor_atom"] == "HIS ND1"]["like"].any())
            ne2 = bool(grp[grp["donor_atom"] == "HIS NE2"]["like"].any())
            states["either" if nd1 and ne2 else "ND1-H only" if nd1 else "NE2-H only" if ne2 else "neither"] += 1
        lines.append(f"| {label} | {sum(states.values())} | {states['ND1-H only']} | {states['NE2-H only']} | {states['either']} | {states['neither']} |")
    neighbors = everything[everything["offset"].abs() == 1].copy()
    neighbors["like"] = hbond_like(neighbors)
    lines += ["", "## Sequence neighbors (the residues just before and after a sulfur residue)", "",
              "`residue_contacts` leaves these out (`min_separation=2`). Backbone-N donors only; H-bond-like by the same test; "
              "cells are H-bond-like / contacts within 4.0 A. Offset +1 is the amide of the next residue, -1 of the previous one.", "",
              "| Sulfur | offset | contacts <= 4.0 A | H-bond-like | N...S of H-bond-like: min / median / max |", "|---|---|---|---|---|"]
    for kind in KINDS:
        for offset in (-1, 1):
            sub = neighbors[(neighbors["kind"] == kind) & (neighbors["offset"] == offset) & (neighbors["donor"] == "backbone N")
                            & neighbors["h_s"].notna() & (neighbors["distance"] <= 4.0)]
            like = sub[sub["like"]]
            spread = f"{like['distance'].min():.2f} / {like['distance'].median():.2f} / {like['distance'].max():.2f}" if len(like) else "none"
            lines.append(f"| {kind} | {offset:+d} | {len(sub)} | {len(like)} | {spread} |")
    lines += ["", "Caveats: the hydrogen is placed, not observed; the sample is atomic-resolution X-ray structures (small, well-ordered "
              "proteins); metal-bound means a metal ion within 2.7 A (Cys) or 3.3 A (Met) of the sulfur in the same model."]
    return "\n".join(lines) + "\n"


def summarize(cache, out_dir):
    cache, out_dir = Path(cache), Path(out_dir)
    meta = json.loads((cache / "chains.json").read_text())
    jobs = [(cache / f"{entry}.cif", chain) for entry, chain in sorted(meta["chains"].items()) if (cache / f"{entry}.cif").exists()]
    with Pool(8) as pool:
        results = pool.map(_job, jobs, chunksize=8)
    frame = pd.DataFrame([row for result in results for row in result])
    out_dir.mkdir(parents=True, exist_ok=True)
    contacts = frame[(frame["donor"] != "sulfur") & frame["h_s"].notna()].copy()
    contacts["h_bond_like"] = hbond_like(contacts)
    contacts.round(3).to_csv(out_dir / "03_sulfur_acceptor_survey_contacts.csv", index=False)
    (out_dir / "03_sulfur_acceptor_survey.md").write_text(report(frame, len(jobs), meta["max_resolution"]))
    print(f"{len(jobs)} structures, {int((frame['donor'] == 'sulfur').sum())} sulfur atoms -> {out_dir}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    fetch_parser = sub.add_parser("fetch", help="query RCSB for representatives and download mmCIF files")
    fetch_parser.add_argument("--cache", default="survey_cache/sulfur")
    fetch_parser.add_argument("--max-resolution", type=float, default=1.2)
    summarize_parser = sub.add_parser("summarize", help="classify sulfurs and write the report")
    summarize_parser.add_argument("--cache", default="survey_cache/sulfur")
    summarize_parser.add_argument("--out", default="docs/scope_checks")
    args = parser.parse_args(argv)
    if args.command == "fetch":
        fetch(args.cache, args.max_resolution)
    else:
        summarize(args.cache, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
