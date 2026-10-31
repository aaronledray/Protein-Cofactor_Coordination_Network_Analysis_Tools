#!/usr/bin/env python3
"""Multi-species plastocyanin check of the conserved-interaction vote (scope check 1, follow-up F).

    python scripts/plastocyanin_vote_check.py fetch --cache survey_cache/plastocyanin
    python scripts/plastocyanin_vote_check.py run   --cache survey_cache/plastocyanin --out docs/scope_checks

``fetch`` lists X-ray plastocyanin entries with a Cu ion at 2.0 A or better (RCSB), takes the first chain of the
plastocyanin entity of each and downloads the files. ``run`` aligns every chain to spinach plastocyanin (1AG6), so
that a position means the same residue everywhere (numbering of 1AG6), then runs the machinery that the heme-copper
oxidase study uses: the anchor audit, residue contacts with the one-vote-per-sequence rule, and cofactor contacts
(with the coordination tiers). One sequence is one UniProt accession, so a protein solved ten times votes once.

The known answer is the type 1 copper site: His37, Cys84 and His87 coordinate the copper in every plastocyanin and
Met92 is a weak axial ligand. The check asks whether the vote recovers that site from the structures without being told,
whether the conserved contacts sit on conserved residues, and where the contacts that are not conserved are.
"""
import argparse
import json
import sys
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from survey_cu_type1 import FILES_URL, GRAPHQL_URL, SEARCH_URL, _post_json  # noqa: E402

from modules import conserved_interactions as ci  # noqa: E402

REFERENCE = "1AG6"
LIGANDS = {37: "His", 84: "Cys", 87: "His", 92: "Met"}
CYANOBACTERIA = {"Q51883", "P46444", "P55020", "P0C178"}          # Leptolyngbya, Nostoc, Synechococcus, Trichormus
THREE_LETTER = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU",
                "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}
THREE_TO_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
                "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def term(attribute, operator, value):
    return {"type": "terminal", "service": "text", "parameters": {"attribute": attribute, "operator": operator, "value": value}}


def search_entries(max_resolution=2.0):
    query = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        term("struct.title", "contains_words", "plastocyanin"), term("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        term("rcsb_entry_info.resolution_combined", "less_or_equal", max_resolution),
        term("rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id", "exact_match", "CU")]},
        "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 1000}}}
    return sorted(item["identifier"] for item in _post_json(SEARCH_URL, query, timeout=120)["result_set"])


def describe_entries(ids):
    """{entry: {chain, uniprot, organism, resolution}} for the plastocyanin entity of each entry."""
    query = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id rcsb_entry_info{resolution_combined} polymer_entities{
        rcsb_polymer_entity_container_identifiers{auth_asym_ids reference_sequence_identifiers{database_name database_accession}}
        rcsb_entity_source_organism{ncbi_scientific_name}}}}"""
    out = {}
    for entry in _post_json(GRAPHQL_URL, {"query": query, "variables": {"ids": list(ids)}})["data"]["entries"]:
        for entity in entry["polymer_entities"]:
            identifiers = entity["rcsb_polymer_entity_container_identifiers"]
            accession = next((r["database_accession"] for r in identifiers.get("reference_sequence_identifiers") or []
                              if r["database_name"] == "UniProt"), None)
            if accession and identifiers["auth_asym_ids"]:
                organism = (entity["rcsb_entity_source_organism"] or [{}])[0].get("ncbi_scientific_name")
                out[entry["rcsb_id"]] = {"chain": identifiers["auth_asym_ids"][0], "uniprot": accession, "organism": organism,
                                         "resolution": entry["rcsb_entry_info"]["resolution_combined"][0]}
                break
    return out


def fetch(cache):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    entries = describe_entries(search_entries())
    for entry in entries:
        target = cache / f"{entry}.cif"
        if not target.exists():
            with urllib.request.urlopen(FILES_URL.format(entry), timeout=120) as response:
                target.write_bytes(response.read())
    (cache / "entries.json").write_text(json.dumps(entries, indent=1))
    print(f"{len(entries)} entries, {len({e['uniprot'] for e in entries.values()})} sequences -> {cache}")


def chain_residues(path, chain_id):
    """[((number, insertion code), one-letter code)] of the standard residues of one chain (first model)."""
    from Bio.PDB import MMCIFParser

    warnings.filterwarnings("ignore")
    model = MMCIFParser(QUIET=True).get_structure("x", str(path))[0]
    return [((int(r.id[1]), r.id[2].strip()), THREE_TO_ONE.get(r.resname.strip().upper(), "X"))
            for r in model[chain_id] if r.id[0] == " "]


def reference_map(reference, residues):
    """{(number, insertion code): reference position} from a global alignment of ``residues`` to ``reference``.

    Both are ``[(key, letter)]``; ``reference`` positions are its residue numbers. Residues aligned to a gap have no entry.
    """
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -10.0, -0.5
    reference_sequence = "".join(letter for _, letter in reference)
    sequence = "".join(letter for _, letter in residues)
    alignment = aligner.align(reference_sequence, sequence)[0]
    mapping = {}
    for (r0, r1), (s0, s1) in zip(*alignment.aligned):
        for offset in range(r1 - r0):
            mapping[residues[s0 + offset][0]] = reference[r0 + offset][0][0]
    return mapping


def build_tables(cache, entries, groups="single"):
    """(atoms, structures, features, residues_by_structure, letters) in reference numbering."""
    cache = Path(cache)
    reference = chain_residues(cache / f"{REFERENCE}.cif", entries[REFERENCE]["chain"])
    atoms, structures, features, letters = [], [], [], {}
    for entry, info in sorted(entries.items()):
        path = cache / f"{entry}.cif"
        residues = chain_residues(path, info["chain"])
        mapping = reference_map(reference, residues)
        table = ci.load_atoms(path, info["chain"], position_map=mapping, structure_id=entry)
        atoms.append(table)
        letters[entry] = {mapping[key]: letter for key, letter in residues if key in mapping}
        metals = ci.load_metals(path, structure_id=entry)
        copper = metals[metals["element"] == "CU"]
        features.append(pd.DataFrame({"structure_id": entry, "feature": "Cu", "atom": "CU", "element": "CU", "kind": "metal",
                                      "x": copper["x"], "y": copper["y"], "z": copper["z"]}).head(1))
        group = "plastocyanin" if groups == "single" else ("cyanobacteria" if info["uniprot"] in CYANOBACTERIA else "eukaryote")
        structures.append((entry, group, info["uniprot"]))
    return (pd.concat(atoms, ignore_index=True), pd.DataFrame(structures, columns=["structure_id", "group", "sequence_id"]),
            pd.concat(features, ignore_index=True), letters)


def all_metals(cache, entries):
    return pd.concat([ci.load_metals(Path(cache) / f"{entry}.cif", structure_id=entry) for entry in sorted(entries)], ignore_index=True)


def organism_of(entries, ids):
    return next((entries[i]["organism"] for i in ids if entries[i]["organism"]), "(not recorded)")


def position_conservation(letters, structures, reference_letters):
    """Per reference position: share of sequences (one per UniProt accession) with the reference residue, and how many have it at all."""
    frame = pd.DataFrame([(sid, p, aa) for sid, mapped in letters.items() for p, aa in mapped.items()], columns=["structure_id", "position", "aa"])
    frame = frame.merge(structures[["structure_id", "sequence_id"]], on="structure_id").drop_duplicates(["sequence_id", "position"])
    total = structures["sequence_id"].nunique()
    rows = []
    for position, group in frame.groupby("position"):
        rows.append({"position": position, "present": len(group), "identical": int((group["aa"] == reference_letters.get(position)).sum()), "total": total})
    return pd.DataFrame(rows)


def run(cache, out_dir):
    cache, out_dir = Path(cache), Path(out_dir)
    entries = json.loads((cache / "entries.json").read_text())
    atoms, structures, features, letters = build_tables(cache, entries)
    reference_letters = letters[REFERENCE]
    n_sequences = structures["sequence_id"].nunique()
    lines = ["# Multi-species plastocyanin check of the conserved-interaction vote", "",
             f"Generated by `scripts/plastocyanin_vote_check.py`. {len(entries)} X-ray structures at 2.0 A or better, {n_sequences} sequences "
             f"(UniProt accessions; a protein solved several times votes once). Positions are the residue numbers of spinach plastocyanin (1AG6) "
             f"after a global alignment of every chain to it.", "",
             "| Sequence | Organism | Structures | Chain length |", "|---|---|---|---|"]
    by_sequence = structures.groupby("sequence_id")["structure_id"].apply(list)
    for sequence_id, ids in sorted(by_sequence.items(), key=lambda kv: organism_of(entries, kv[1])):
        lines.append(f"| {sequence_id} | {organism_of(entries, ids)} | {', '.join(ids)} | {len(letters[ids[0]])} aligned |")

    # ---- audits
    audit = ci.anchor_audit(atoms, features[["structure_id", "x", "y", "z"]], 84, structures)
    backbone = ci.backbone_only_audit(atoms)
    lines += ["", "## Input audits", "", f"- Anchor audit (Cys84 within 8 A of the copper): {int((~audit['ok']).sum())} of {len(audit)} structures flagged.",
              f"- Backbone-only audit: {int(backbone['backbone_only'].sum())} of {len(backbone)} structures flagged."]

    # ---- position conservation
    conservation = position_conservation(letters, structures, reference_letters)
    conservation["fraction_identical"] = (conservation["identical"] / conservation["total"]).round(3)
    invariant = set(conservation.loc[conservation["identical"] == conservation["total"], "position"])
    lines += ["", f"- Residues identical to spinach in all {n_sequences} sequences: {len(invariant)} of {len(reference_letters)} positions.",
              "- Ligand positions: " + ", ".join(f"{name}{pos}: {int(conservation.set_index('position').loc[pos, 'identical'])}/{n_sequences} identical"
                                                for pos, name in LIGANDS.items()) + "."]

    # ---- cofactor contacts and tiers
    cofactor = ci.cofactor_contacts(atoms, features, coordination_model="by_pair")
    cofactor["bond"] = cofactor["coord_tier"] == "bond"
    cofactor["weak"] = cofactor["coord_tier"] == "weak"
    cofactor["ligated"] = cofactor["bond"] | cofactor["weak"]
    cofactor_vote = ci.sequence_vote(cofactor, structures, ["feature", "position"], ["coord", "bond", "weak", "ligated", "contact", "polar"])
    lines += ["", "## Copper contacts (one vote per sequence)", "",
              "A residue is listed when it passes the vote for at least one class. `bond` and `weak` are the coordination tiers "
              "(`docs/scope_checks/02_cu_type1_survey.md`); `ligated` is bond or weak; `coord` is the single 2.8 A cutoff. The tier cutoffs come from a survey of type 1 copper sites that includes these proteins (`docs/scope_checks/02_cu_type1_survey.md`), so recovering their ligands is a consistency check, not an independent test.", "",
              "| Position | Spinach residue | coord | bond | weak | ligated | contact | identical in sequences |", "|---|---|---|---|---|---|---|---|"]
    flags = ["coord", "bond", "weak", "ligated", "contact"]
    cofactor_vote = cofactor_vote.sort_values("position")
    for row in cofactor_vote.itertuples():
        if not any(getattr(row, f"conserved_{f}") for f in flags):
            continue
        identical = int(conservation.set_index("position").loc[row.position, "identical"])
        cells = [f"{getattr(row, f'{f}_plastocyanin'):.2f}{' *' if getattr(row, f'conserved_{f}') else ''}" for f in flags]
        lines.append(f"| {row.position} | {reference_letters.get(row.position, '?')}{row.position} | " + " | ".join(cells) + f" | {identical}/{n_sequences} |")
    lines += ["", "Cells are the fraction of sequences voting positive; `*` marks a conserved class (at least 90% of at least 5 sequences)."]

    # ---- residue contacts
    metals = all_metals(cache, entries)
    contacts = ci.conserved_contacts(atoms, structures, sulfur_acceptor=True, metals=metals)
    conserved = contacts[contacts["conserved_contact"]].copy()
    conserved["identical_a"] = conserved["position_a"].map(conservation.set_index("position")["identical"]).fillna(0).astype(int)
    conserved["identical_b"] = conserved["position_b"].map(conservation.set_index("position")["identical"]).fillna(0).astype(int)
    all_pairs = len(contacts)
    both_invariant = int(((conserved["identical_a"] == n_sequences) & (conserved["identical_b"] == n_sequences)).sum())
    neither = int(((conserved["identical_a"] < n_sequences) & (conserved["identical_b"] < n_sequences)).sum())
    lines += ["", "## Residue-residue contacts", "",
              f"Candidate pairs (in contact in at least one structure): {all_pairs}. Conserved contacts: {len(conserved)}. Of these, both residues "
              f"identical in every sequence: {both_invariant}; neither residue invariant (packing that survives substitution): {neither}.", ""]
    ligand_set = set(LIGANDS)
    touching = conserved[conserved["position_a"].isin(ligand_set) | conserved["position_b"].isin(ligand_set)]
    lines += [f"Conserved contacts that involve a copper ligand residue: {len(touching)}.", ""]
    lines += ["### N-H...S contacts to the ligand sulfur (opt-in `sulfur_acceptor` flags, metals passed)", "",
              "Pairs that contain Cys84 and are conserved in contact; cells are the fraction of sequences voting positive among contact-positive structures "
              "(`*` = conserved).", "", "| Pair | contact | polar | polar_s | polar_s_bb | polar_s_geom |", "|---|---|---|---|---|---|"]
    cys = conserved[(conserved["position_a"] == 84) | (conserved["position_b"] == 84)].sort_values(["position_a", "position_b"])
    for row in cys.itertuples():
        def cell(flag):
            value = getattr(row, f"{flag}_plastocyanin")
            return f"{value:.2f}{' *' if getattr(row, f'conserved_{flag}') else ''}" if pd.notna(value) else "-"
        lines.append(f"| {reference_letters.get(row.position_a, '?')}{row.position_a} - {reference_letters.get(row.position_b, '?')}{row.position_b} | "
                     + " | ".join(cell(f) for f in ("contact", "polar", "polar_s", "polar_s_bb", "polar_s_geom")) + " |")
    lines.append("")
    out_dir.mkdir(parents=True, exist_ok=True)
    keep = ["position_a", "position_b", "contact_plastocyanin", "contact_n_plastocyanin", "closest", "identical_a", "identical_b"]
    conserved.assign(residue_a=conserved["position_a"].map(reference_letters), residue_b=conserved["position_b"].map(reference_letters))[
        ["position_a", "residue_a", "position_b", "residue_b", "contact_plastocyanin", "contact_n_plastocyanin", "identical_a", "identical_b"]
        + [c for c in conserved.columns if c.startswith("conserved_") and c != "conserved_contact"]].to_csv(
        out_dir / "04_plastocyanin_vote_conserved_contacts.csv", index=False)

    # ---- sequence-neighbor contacts (opt-in)
    classes = ci.NEIGHBOR_CHEMISTRY_ALL
    with_neighbors = ci.conserved_contacts(atoms, structures, sulfur_acceptor=True, metals=metals, neighbor_contacts=True, neighbor_chemistry=classes)
    structure_rows = ci.residue_contacts(atoms, neighbor_contacts=True, neighbor_chemistry=classes, metals=metals)
    neighbor_pairs = set(map(tuple, structure_rows.loc[structure_rows["neighbor"], ["position_a", "position_b"]].drop_duplicates().to_numpy()))
    found = with_neighbors[[(a, b) in neighbor_pairs for a, b in zip(with_neighbors["position_a"], with_neighbors["position_b"])]]
    kept = found[found["conserved_contact"]].sort_values("position_a")
    background = pd.read_csv(out_dir / "06_neighbor_chemistry_background.csv").set_index(["resname_a", "resname_b"]) if (out_dir / "06_neighbor_chemistry_background.csv").exists() else None
    lines += ["### Sequence-neighbor contacts (opt-in `neighbor_contacts`, all four classes)", "",
              "Residues one apart in sequence are left out of the table above (`min_separation=2`). `neighbor_contacts` judges them per atom pair "
              "(at least 4 bonds apart through the peptide bond; a polar pair is a hydrogen bond, with the placed hydrogen pointing at the acceptor; "
              "hydrophobic and aromatic pairs involve side-chain atoms beyond C-beta). "
              f"Pairs with such an atom pair in at least one structure: {len(found)}; conserved (in at least 90% of sequences): {len(kept)}. "
              "`background` is the share of that residue-type pair that shows the class in 2,048 atomic-resolution structures "
              "(`docs/scope_checks/06_neighbor_chemistry_survey.md`): a high value means the contact is common whatever the protein.", "",
              "| Pair | contact | classes voted conserved (fraction of sequences) | background |", "|---|---|---|---|"]
    table_rows = []
    for row in kept.itertuples():
        voted = [f"{c} ({getattr(row, f'{c}_plastocyanin'):.2f})" for c in classes if getattr(row, f"conserved_{c}")]
        names = (THREE_LETTER.get(reference_letters.get(row.position_a)), THREE_LETTER.get(reference_letters.get(row.position_b)))
        rates = []
        if background is not None and names in background.index:
            rates = [f"{c} {100 * background.loc[names, f'{c}_rate']:.1f}%" for c in classes if getattr(row, f"conserved_{c}")]
        lines.append(f"| {reference_letters.get(row.position_a, '?')}{row.position_a} - {reference_letters.get(row.position_b, '?')}{row.position_b} | "
                     f"{row.contact_plastocyanin:.2f} | {', '.join(voted) or 'none'} | {', '.join(rates) or '-'} |")
        table_rows.append((row.position_a, reference_letters.get(row.position_a), row.position_b, reference_letters.get(row.position_b),
                           row.contact_plastocyanin, "; ".join(voted), "; ".join(rates)))
    if not len(kept):
        lines.append("| (none conserved) | | | |")
    near = found[(~found["conserved_contact"]) & (found["contact_plastocyanin"] >= 0.5)].sort_values("contact_plastocyanin", ascending=False)
    lines += ["", f"Seen in at least half of the sequences but not conserved ({len(near)}; fraction of sequences with the contact):", ""]
    lines.append(", ".join(f"{reference_letters.get(r.position_a, '?')}{r.position_a}-{reference_letters.get(r.position_b, '?')}{r.position_b} ({r.contact_plastocyanin:.2f})"
                           for r in near.itertuples()) or "none")
    lines.append("")
    pd.DataFrame(table_rows, columns=["position_a", "residue_a", "position_b", "residue_b", "contact_plastocyanin", "conserved_classes", "background_rate"]).to_csv(
        out_dir / "04_plastocyanin_vote_conserved_neighbor_contacts.csv", index=False)

    # ---- two lineages
    atoms2, structures2, features2, _ = build_tables(cache, entries, groups="lineages")
    counts = structures2.groupby("group")["sequence_id"].nunique().to_dict()
    cofactor2 = ci.cofactor_contacts(atoms2, features2, coordination_model="by_pair")
    cofactor2["bond"] = cofactor2["coord_tier"] == "bond"
    vote2 = ci.sequence_vote(cofactor2, structures2, ["feature", "position"], ["bond", "contact"], min_sequences=min(counts.values()))
    shared = vote2[vote2["conserved_bond"]]["position"].tolist()
    contacts2 = ci.conserved_contacts(atoms2, structures2, min_sequences=min(counts.values()))
    lines += ["## Two lineages", "", f"Eukaryotic ({counts.get('eukaryote')} sequences) and cyanobacterial ({counts.get('cyanobacteria')} sequences) plastocyanins "
              f"voted separately, requiring support in both (`min_sequences` lowered to {min(counts.values())} so the smaller group can vote).", "",
              "- Copper bond-tier residues conserved in both lineages: " + (", ".join(f"{reference_letters.get(p, '?')}{p}" for p in shared) or "none") + ".",
              f"- Conserved contacts in both lineages: {int(contacts2['conserved_contact'].sum())} (single group: {len(conserved)}).", ""]
    (out_dir / "04_plastocyanin_vote.md").write_text("\n".join(lines) + "\n")
    print(f"{len(entries)} structures, {n_sequences} sequences; {len(conserved)} conserved contacts -> {out_dir}")
    return {"conserved": conserved, "cofactor_vote": cofactor_vote, "conservation": conservation, "audit": audit, "contacts": contacts}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    fetch_parser = sub.add_parser("fetch")
    fetch_parser.add_argument("--cache", default="survey_cache/plastocyanin")
    run_parser = sub.add_parser("run")
    run_parser.add_argument("--cache", default="survey_cache/plastocyanin")
    run_parser.add_argument("--out", default="docs/scope_checks")
    args = parser.parse_args(argv)
    fetch(args.cache) if args.command == "fetch" else run(args.cache, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
