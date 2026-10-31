#!/usr/bin/env python3
"""Survey type 1 (blue) copper ligand distances in high-resolution PDB structures.

Evidence for the per-metal, per-donor-class coordination cutoffs (TODO item B,
docs/scope_checks/02_cu_type1_survey.md). Two steps, kept separate so the slow download runs once:

    python scripts/survey_cu_type1.py fetch     --cache survey_cache/cu_type1
    python scripts/survey_cu_type1.py summarize --cache survey_cache/cu_type1 --out docs/scope_checks

``fetch`` asks RCSB for X-ray entries that contain a Cu ion at or below ``--max-resolution`` and downloads
their mmCIF files into the cache (about 1 GB for 1.8 A; git-ignored, never commit it). ``summarize`` finds
type 1 sites (one Cys SG within 2.7 A and two His N within 2.5 A of the same Cu), groups entries by UniProt
accession so one protein counts once (median over its entries), and writes a Markdown report and a per-protein CSV.

The axial methionine is identified by position: a Met (or MSE) S within 3.6 A of the Cu that sits four or five
residues after the C-terminal His ligand (His87 and Met92 in plastocyanin). The rule is rough, and the report
says how many Met fall outside it.
"""
import argparse
import json
import sys
import urllib.request
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
FILES_URL = "https://files.rcsb.org/download/{}.cif"

CYS_MAX = 2.7        # Cu-SG window used only to recognise a type 1 site, not a proposed cutoff
HIS_MAX = 2.5
NEAR = 3.6           # report every N/O/S donor inside this distance
AXIAL_OFFSETS = (4, 5)


def _post_json(url, payload, timeout=120):
    request = urllib.request.Request(url, json.dumps(payload).encode(), {"Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)


def search_entries(max_resolution=1.8):
    """PDB ids of X-ray entries with a CU ion at or below ``max_resolution`` angstroms."""
    def term(attribute, operator, value):
        return {"type": "terminal", "service": "text",
                "parameters": {"attribute": attribute, "operator": operator, "value": value}}

    query = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        term("rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id", "exact_match", "CU"),
        term("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        term("rcsb_entry_info.resolution_combined", "less_or_equal", max_resolution)]},
        "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 10000}}}
    result = _post_json(SEARCH_URL, query, timeout=60)
    return sorted(item["identifier"] for item in result["result_set"])


def fetch(cache, max_resolution=1.8, workers=8):
    from concurrent.futures import ThreadPoolExecutor

    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    ids = search_entries(max_resolution)
    todo = [i for i in ids if not (cache / f"{i}.cif").exists()]

    def get(pdb_id):
        with urllib.request.urlopen(FILES_URL.format(pdb_id), timeout=120) as response:
            data = response.read()
        (cache / f"{pdb_id}.cif").write_bytes(data)

    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(get, todo))
    (cache / "entries.json").write_text(json.dumps({"max_resolution": max_resolution, "ids": ids}))
    print(f"{len(ids)} entries at <= {max_resolution} A; downloaded {len(todo)} new into {cache}")


def type1_sites(cif_path):
    """Type 1 copper sites in one mmCIF file: one dict per Cu with its N/O/S donors within ``NEAR``."""
    from Bio.PDB import MMCIFParser, NeighborSearch

    warnings.filterwarnings("ignore")
    model = MMCIFParser(QUIET=True).get_structure("x", str(cif_path))[0]
    atoms = list(model.get_atoms())
    coppers = [a for a in atoms if (a.element or "").upper() == "CU"]
    if not coppers:
        return []
    search = NeighborSearch(atoms)
    pdb_id = Path(cif_path).stem[:4].upper()
    sites = []
    for cu in coppers:
        donors = []
        for atom in search.search(cu.coord, NEAR):
            residue = atom.get_parent()
            if residue.id[0] == "W" or atom.get_parent() is cu.get_parent():
                continue
            if residue.id[0] != " " and residue.resname != "MSE":
                continue
            if (atom.element or "").upper() not in ("N", "O", "S", "SE"):
                continue
            donors.append({"resname": residue.resname, "atom": atom.get_id(), "chain": residue.get_parent().id,
                           "resnum": int(residue.id[1]),
                           "distance": round(float(np.linalg.norm(atom.coord - cu.coord)), 3)})
        cys = [d for d in donors if d["resname"] == "CYS" and d["atom"] == "SG" and d["distance"] <= CYS_MAX]
        his = [d for d in donors if d["resname"] == "HIS" and d["atom"] in ("ND1", "NE2") and d["distance"] <= HIS_MAX]
        if cys and len(his) >= 2:
            sites.append({"pdb": pdb_id, "chain": cu.get_parent().get_parent().id, "donors": donors})
    return sites


def characterise(site):
    """Reduce a site to the distances the survey reports (None where the donor is absent)."""
    donors = site["donors"]
    his = sorted((d for d in donors if d["resname"] == "HIS" and d["atom"] in ("ND1", "NE2") and d["distance"] <= HIS_MAX),
                 key=lambda d: d["resnum"])
    cys = min(d["distance"] for d in donors if d["resname"] == "CYS" and d["atom"] == "SG" and d["distance"] <= CYS_MAX)
    last_his = his[-1]["resnum"]
    met = [d for d in donors if d["resname"] in ("MET", "MSE") and d["atom"] in ("SD", "SE")]
    axial = [d["distance"] for d in met if d["resnum"] - last_his in AXIAL_OFFSETS]
    other = [d["distance"] for d in met if d["resnum"] - last_his not in AXIAL_OFFSETS]
    amide_o = [d["distance"] for d in donors if d["resname"] in ("GLN", "GLU", "ASN") and d["atom"].startswith("O")]
    backbone_o = [d["distance"] for d in donors if d["atom"] == "O"]
    smallest = lambda values: min(values) if values else None
    ordered = sorted(d["distance"] for d in his)[:2]
    return {"pdb": site["pdb"], "chain": site["chain"], "cys_s": cys, "his_n_near": ordered[0], "his_n_far": ordered[1],
            "met_s_axial": smallest(axial), "met_s_other": smallest(other),
            "amide_o": smallest(amide_o), "backbone_o": smallest(backbone_o)}


def _extract(path):
    return [characterise(s) for s in type1_sites(path)]


def protein_groups(pdb_ids):
    """(pdb id, auth chain) -> {'protein': UniProt accession or sequence stub, 'name': description, 'organism': ...}."""
    query = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id polymer_entities{
        rcsb_polymer_entity_container_identifiers{auth_asym_ids reference_sequence_identifiers{database_name database_accession}}
        entity_poly{pdbx_seq_one_letter_code_can} rcsb_entity_source_organism{ncbi_scientific_name}
        rcsb_polymer_entity{pdbx_description}}}}"""
    out = {}
    ids = sorted(pdb_ids)
    for start in range(0, len(ids), 100):
        reply = _post_json(GRAPHQL_URL, {"query": query, "variables": {"ids": ids[start:start + 100]}})
        for entry in reply["data"]["entries"]:
            for entity in entry["polymer_entities"]:
                identifiers = entity["rcsb_polymer_entity_container_identifiers"]
                accessions = [r["database_accession"] for r in identifiers.get("reference_sequence_identifiers") or []
                              if r["database_name"] == "UniProt"]
                protein = accessions[0] if accessions else "seq:" + entity["entity_poly"]["pdbx_seq_one_letter_code_can"][:40]
                organism = (entity["rcsb_entity_source_organism"] or [{}])[0].get("ncbi_scientific_name")
                for chain in identifiers["auth_asym_ids"]:
                    out[(entry["rcsb_id"], chain)] = {"protein": protein, "organism": organism,
                                                       "name": entity["rcsb_polymer_entity"]["pdbx_description"]}
    return out


FAMILIES = (("nitrite reductase", "nitrite reductase"), ("plastocyanin", "plastocyanin"), ("pseudoazurin", "pseudoazurin"),
            ("azurin", "azurin"), ("laccase", "laccase / multicopper oxidase"), ("oxidase", "laccase / multicopper oxidase"),
            ("multi", "laccase / multicopper oxidase"), ("amicyanin", "amicyanin"), ("rusticyanin", "rusticyanin"))


def family(name):
    lowered = (name or "").lower()
    for key, label in FAMILIES:
        if key in lowered:
            return label
    return "other"


def per_protein(rows, groups):
    """One row per protein: the median of each distance over that protein's entries."""
    frame = pd.DataFrame(rows)
    keys = [groups.get((r.pdb, r.chain)) for r in frame.itertuples()]
    frame = frame[[k is not None for k in keys]].copy()
    keys = [k for k in keys if k is not None]
    frame["protein"] = [k["protein"] for k in keys]
    frame["name"] = [k["name"] for k in keys]
    frame["organism"] = [k["organism"] for k in keys]
    distances = ["cys_s", "his_n_near", "his_n_far", "met_s_axial", "met_s_other", "amide_o", "backbone_o"]
    grouped = frame.groupby("protein")
    table = grouped[distances].median()
    table["entries"] = grouped.size()
    table["name"] = grouped["name"].first()
    table["organism"] = grouped["organism"].first()
    table["family"] = table["name"].map(family)
    return table.reset_index()


def _quantiles(values, label):
    values = values.dropna()
    q = np.percentile(values, [0, 5, 50, 95, 100])
    return f"| {label} | {len(values)} | {q[0]:.2f} | {q[1]:.2f} | {q[2]:.2f} | {q[3]:.2f} | {q[4]:.2f} |"


def report(table, n_entries, n_sites, max_resolution):
    lines = [
        "# Type 1 copper ligand distances (PDB survey)", "",
        f"Generated by `scripts/survey_cu_type1.py`. X-ray entries with a Cu ion at <= {max_resolution} A: {n_entries}; "
        f"type 1 sites (Cys SG <= {CYS_MAX} A and two His N <= {HIS_MAX} A): {n_sites}; distinct proteins "
        f"(UniProt accession; median over a protein's entries): {len(table)}.", "",
        "| Distance (A) | proteins | min | p5 | median | p95 | max |", "|---|---|---|---|---|---|---|",
        _quantiles(table["cys_s"], "Cu-S(Cys)"), _quantiles(table["his_n_near"], "Cu-N(His), nearer"),
        _quantiles(table["his_n_far"], "Cu-N(His), farther"), _quantiles(table["met_s_axial"], "Cu-S(Met), His+4/+5"),
        _quantiles(table["amide_o"], "Cu-O(Gln/Glu/Asn), nearest within 3.6"),
        _quantiles(table["backbone_o"], "Cu-O(backbone), nearest within 3.6"), "",
        f"Proteins with a Met at His+4/+5 inside 3.6 A: {table['met_s_axial'].notna().sum()}; with a Met inside 3.6 A at another "
        f"offset only: {(table['met_s_axial'].isna() & table['met_s_other'].notna()).sum()}; with no Met inside 3.6 A: "
        f"{(table['met_s_axial'].isna() & table['met_s_other'].isna()).sum()}.", "",
        "## Axial Met by family", "", "| Family | proteins | min | median | max |", "|---|---|---|---|---|"]
    axial = table.dropna(subset=["met_s_axial"])
    for name, group in sorted(axial.groupby("family"), key=lambda kv: -len(kv[1])):
        values = group["met_s_axial"]
        lines.append(f"| {name} | {len(group)} | {values.min():.2f} | {values.median():.2f} | {values.max():.2f} |")
    histogram, edges = np.histogram(axial["met_s_axial"], bins=np.arange(2.3, 3.7, 0.1))
    lines += ["", "## Axial Met histogram (proteins per 0.1 A)", "", "```"]
    lines += [f"{lo:.1f}-{lo + 0.1:.1f} {'#' * count} {count}" for lo, count in zip(edges, histogram)]
    lines += ["```", "",
              "Caveats: refinement can restrain Cu-S(Met); entries mix redox states and mutants; the His+4/+5 test is a "
              "positional rule and can miss an axial Met with an unusual spacing."]
    return "\n".join(lines) + "\n"


def summarize(cache, out_dir):
    cache, out_dir = Path(cache), Path(out_dir)
    meta = json.loads((cache / "entries.json").read_text())
    files = sorted(cache.glob("*.cif"))
    with Pool(8) as pool:
        results = pool.map(_extract, files, chunksize=8)
    rows = [row for result in results for row in result]
    groups = protein_groups({r["pdb"] for r in rows})
    table = per_protein(rows, groups)
    out_dir.mkdir(parents=True, exist_ok=True)
    table.round(3).to_csv(out_dir / "02_cu_type1_survey_proteins.csv", index=False)
    (out_dir / "02_cu_type1_survey.md").write_text(report(table, len(meta["ids"]), len(rows), meta["max_resolution"]))
    print(f"{len(rows)} type 1 sites, {len(table)} proteins -> {out_dir}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    fetch_parser = sub.add_parser("fetch", help="query RCSB and download mmCIF files")
    fetch_parser.add_argument("--cache", default="survey_cache/cu_type1")
    fetch_parser.add_argument("--max-resolution", type=float, default=1.8)
    summarize_parser = sub.add_parser("summarize", help="extract sites and write the report")
    summarize_parser.add_argument("--cache", default="survey_cache/cu_type1")
    summarize_parser.add_argument("--out", default="docs/scope_checks")
    args = parser.parse_args(argv)
    if args.command == "fetch":
        fetch(args.cache, args.max_resolution)
    else:
        summarize(args.cache, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
