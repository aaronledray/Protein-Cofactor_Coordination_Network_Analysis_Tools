#!/usr/bin/env python3
"""Shuffle control for the conserved-interaction vote (scope check 2).

The vote asks which residue pairs touch in nearly every sequence of every group. If the position labels that make
structures comparable carry no information, nothing should pass. This script relabels positions within each structure
(the atoms, their coordinates and their sequence order stay as they are, so every structure keeps its own geometry)
and runs the vote again:

* ``full``: every structure's position labels are permuted at random. The conserved set should collapse to zero.
* ``local-w``: labels are permuted only within blocks of ``w`` consecutive positions. Shows how much a small local
  misalignment costs.
* ``sequences-k``: ``k`` randomly chosen sequences are fully scrambled (one permutation per sequence, shared by its
  structures) and the rest are left intact. Shows how many misaligned sequences the 90% rule and the group minimum tolerate.

    python scripts/shuffle_control.py --dataset plastocyanin --replicates 20
    python scripts/shuffle_control.py --dataset hco --replicates 3 --workers 1

Memory: one vote on the heme-copper oxidase set (about 1 million atoms, 249 structures) takes about 50 seconds and its peak memory has not
yet been measured; a first run with 6 workers exhausted a 16 GB laptop. Keep ``--workers`` at 1 until the peak is known.

``plastocyanin`` needs the files fetched by ``plastocyanin_vote_check.py fetch``. ``hco`` reads the heme-copper oxidase
manuscript's position map read-only (``HCO_MANUSCRIPT``, default the sibling checkout) and is skipped when absent.
Writes ``docs/scope_checks/07_shuffle_control.md``.
"""
import argparse
import json
import os
import pickle
import sys
import multiprocessing
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True        # never write bytecode next to the manuscript's code
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from modules import conserved_interactions as ci  # noqa: E402

BASE_SEED = 20261031
DATA = {}


def permute_labels(atoms: pd.DataFrame, rng: np.random.Generator, *, window=None, only=None) -> pd.DataFrame:
    """``atoms`` with the position labels of each structure permuted.

    Whole structure by default; with ``window`` the permutation stays inside blocks of that many consecutive positions
    (in sorted order); ``only`` limits it to a set of structure ids. Coordinates, atom names and ``seq_index`` do not change.
    """
    out = atoms.copy()
    positions = out["position"].to_numpy().copy()
    for sid, index in out.groupby("structure_id", sort=False).indices.items():
        if only is not None and sid not in only:
            continue
        labels = np.unique(positions[index])
        if window is None:
            shuffled = rng.permutation(labels)
        else:
            shuffled = labels.copy()
            for start in range(0, len(labels), window):
                shuffled[start:start + window] = rng.permutation(labels[start:start + window])
        positions[index] = pd.Series(positions[index]).map(dict(zip(labels, shuffled))).to_numpy()
    out["position"] = positions
    return out


def scramble_sequences(atoms: pd.DataFrame, structures: pd.DataFrame, k: int, rng: np.random.Generator):
    """Scramble ``k`` randomly chosen sequences (one permutation per sequence over all its positions); returns (atoms, chosen ids)."""
    if k > structures["sequence_id"].nunique():
        raise ValueError(f"cannot scramble {k} sequences: there are {structures['sequence_id'].nunique()}")
    chosen = list(rng.choice(sorted(structures["sequence_id"].unique()), size=k, replace=False))
    out = atoms.copy()
    positions = out["position"].to_numpy().copy()
    for sequence in chosen:
        members = set(structures.loc[structures["sequence_id"] == sequence, "structure_id"])
        rows = np.nonzero(out["structure_id"].isin(members).to_numpy())[0]
        labels = np.unique(positions[rows])
        mapping = dict(zip(labels, rng.permutation(labels)))
        positions[rows] = pd.Series(positions[rows]).map(mapping).to_numpy()
    out["position"] = positions
    return out, chosen


def conserved_set(atoms, structures, groups, **kwargs):
    votes = ci.conserved_contacts(atoms, structures, groups=groups, **kwargs)
    return set(map(tuple, votes.loc[votes["conserved_contact"], ["position_a", "position_b"]].to_numpy().tolist())), len(votes)


# ----------------------------------------------------------------------------- datasets
def load_plastocyanin(cache):
    import plastocyanin_vote_check as pv

    entries = json.loads((Path(cache) / "entries.json").read_text())
    atoms, structures, features, _ = pv.build_tables(cache, entries)
    return {"name": "plastocyanin (12 sequences, 30 structures)", "atoms": atoms, "structures": structures, "groups": ["plastocyanin"],
            "kwargs": {}, "features": features, "reference": "168 conserved contacts"}


def load_hco(controls):
    controls = Path(controls)
    sys.path.insert(0, str(controls))
    import misaligned                                   # the manuscript's own code, imported read-only

    out = controls / "output"
    auth = misaligned.apply_keyed(pickle.load(open(out / "auth_map.pkl", "rb"))["structs"])
    cache = pickle.load(open(out / "cache.pkl", "rb"))
    seq_pos = {(s["entry"], s["seq"]): {c: v["pos"] for c, v in s["colinfo"].items()} for s in cache["structs"]}
    sid_of = lambda k: f"{k[0]}_{k[1][:8]}"
    rows, meta = [], []
    for k, v in auth.items():
        meta.append((sid_of(k), v["cls"], k[1]))
        for col, res in v["cols"].items():
            for name, elem, xyz, *_ in res["atoms"]:
                rows.append((sid_of(k), col, seq_pos[k][col], res["resname"], name, elem, float(xyz[0]), float(xyz[1]), float(xyz[2])))
    atoms = pd.DataFrame(rows, columns=["structure_id", "position", "seq_index", "resname", "atom", "element", "x", "y", "z"])
    structures = pd.DataFrame(meta, columns=["structure_id", "group", "sequence_id"])
    return {"name": "heme-copper oxidases (54 sequences, 249 structures; groups A, B1, C)", "atoms": atoms, "structures": structures,
            "groups": ["A", "B1", "C"], "kwargs": {}, "features": None, "reference": "490 conserved contacts"}


def plan(replicates, windows, ks):
    tasks = [("full", None, r) for r in range(replicates)]
    tasks += [(f"local-{w}", w, r) for w in windows for r in range(replicates)]
    tasks += [(f"sequences-{k}", k, r) for k in ks for r in range(replicates)]
    return tasks


def run_task(args):
    scheme, parameter, replicate, index = args
    data = DATA["current"]
    rng = np.random.default_rng([BASE_SEED, index])
    if scheme == "full":
        atoms = permute_labels(data["atoms"], rng)
    elif scheme.startswith("local"):
        atoms = permute_labels(data["atoms"], rng, window=parameter)
    else:
        atoms, _ = scramble_sequences(data["atoms"], data["structures"], parameter, rng)
    found, candidates = conserved_set(atoms, data["structures"], data["groups"], **data["kwargs"])
    return scheme, replicate, found, candidates


def run_dataset(data, replicates, windows, ks, workers):
    DATA["current"] = data
    real, candidates = conserved_set(data["atoms"], data["structures"], data["groups"], **data["kwargs"])
    tasks = [(s, p, r, i) for i, (s, p, r) in enumerate(plan(replicates, windows, ks))]
    with multiprocessing.get_context("fork").Pool(workers) as pool:                  # fork: the workers inherit the loaded tables
        results = pool.map(run_task, tasks, chunksize=1)
    summary = []
    for scheme in dict.fromkeys(s for s, *_ in tasks):
        sets = [found for s, _, found, _ in results if s == scheme]
        counts = np.array([len(f) for f in sets])
        shared = np.array([len(f & real) for f in sets])
        summary.append({"scheme": scheme, "replicates": len(sets), "mean": counts.mean(), "min": counts.min(), "max": counts.max(),
                        "mean_shared": shared.mean(), "max_shared": shared.max()})
    return real, candidates, pd.DataFrame(summary)


def write_report(results, out_dir):
    lines = ["# Shuffle control of the conserved-interaction vote", "",
             f"Generated by `scripts/shuffle_control.py` (scope check 2). Position labels are permuted within each structure (atoms, coordinates "
             f"and sequence order unchanged), the vote is run again with the unchanged rules, and the number of conserved contacts is compared with the real one. "
             f"Seeds are fixed; the number of replicates per scheme is in each table.", ""]
    for data, real, candidates, summary in results:
        lines += [f"## {data['name']}", "", f"Real data: {len(real)} conserved contacts among {candidates} candidate pairs ({data['reference']} in the reference).", "",
                  "| Scheme | replicates | conserved contacts: mean | min | max | of those, also in the real set (mean / max) |", "|---|---|---|---|---|---|"]
        for r in summary.itertuples():
            lines.append(f"| {r.scheme} | {r.replicates} | {r.mean:.1f} | {r.min} | {r.max} | {r.mean_shared:.1f} / {r.max_shared} |")
        lines.append("")
    lines += ["Schemes: `full` permutes every structure's labels; `local-w` permutes inside blocks of w consecutive positions; `sequences-k` scrambles k "
              "randomly chosen sequences completely and leaves the others intact. With zero conserved contacts in all n replicates of `full`, the empirical "
              "p value of the real count is below 1/(n + 1). The null tested here is relabelling positions within each structure; it does not test other nulls (unrelated proteins, a different alignment), and the oxidase rows use few replicates, so the sequences-k means depend on which groups the scrambled sequences fall in (a group of exactly five sequences has no tolerance)."]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "07_shuffle_control.md").write_text("\n".join(lines) + "\n")


def copper_task(args):
    scheme, parameter, index = args
    data = DATA["current"]
    rng = np.random.default_rng([BASE_SEED, 1000 + index])
    if scheme == "real":
        atoms = data["atoms"]
    elif scheme == "full":
        atoms = permute_labels(data["atoms"], rng)
    else:
        atoms, _ = scramble_sequences(data["atoms"], data["structures"], parameter, rng)
    contacts = ci.cofactor_contacts(atoms, data["features"], coordination_model="by_pair")
    contacts["bond"] = contacts["coord_tier"] == "bond"
    contacts["ligated"] = contacts["coord_tier"].isin(["bond", "weak"])
    vote = ci.sequence_vote(contacts, data["structures"], ["feature", "position"], ["bond", "ligated"], groups=data["groups"])
    return scheme, set(vote.loc[vote["conserved_bond"], "position"]), set(vote.loc[vote["conserved_ligated"], "position"])


def copper_section(data, replicates, workers):
    """Markdown section: do the copper ligands (known answer: His37, Cys84, His87, Met92) survive the shuffles?"""
    DATA["current"] = data
    tasks = [("real", None, 0)] + [("full", None, r) for r in range(replicates)] + [("sequences-1", 1, 100 + r) for r in range(replicates)]
    tasks += [("sequences-2", 2, 200 + r) for r in range(replicates)]
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        results = pool.map(copper_task, tasks, chunksize=1)
    lines = ["## Copper site under the shuffle (plastocyanin, known answer)", "",
             "Conserved copper bond-tier residues (and bond-or-weak residues) among the sequences, after `cofactor_contacts(coordination_model='by_pair')` "
             "and the one-vote-per-sequence rule. Known answer: His37, Cys84, His87 (bond) and Met92 (bond or weak).", "",
             "| Scheme | replicates | bond: His37, Cys84, His87 all conserved | bond: residues conserved (mean) | bond-or-weak: residues conserved (mean) |", "|---|---|---|---|---|"]
    for scheme in dict.fromkeys(s for s, *_ in tasks):
        rows = [(b, l) for s, b, l in results if s == scheme]
        everything = sum(1 for b, _ in rows if {37, 84, 87} <= b)
        lines.append(f"| {scheme} | {len(rows)} | {everything} of {len(rows)} | {np.mean([len(b) for b, _ in rows]):.1f} | {np.mean([len(l) for _, l in rows]):.1f} |")
    return "\n".join(lines) + "\n"


def append_section(out_dir, section):
    path = Path(out_dir) / "07_shuffle_control.md"
    text = path.read_text() if path.exists() else ""
    marker = "## Copper site under the shuffle"
    if marker in text:
        text = text[:text.index(marker)].rstrip("\n") + "\n"
    path.write_text(text.rstrip("\n") + "\n\n" + section)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dataset", choices=["plastocyanin", "hco", "both"], default="both")
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--hco-replicates", type=int, default=None, help="replicates for the heme-copper oxidase set (default: --replicates); one vote there takes about a minute")
    parser.add_argument("--workers", type=int, default=1, help="parallel processes; each holds a copy of the tables and a vote needs gigabytes on the heme-copper oxidase set, so measure first")
    parser.add_argument("--windows", type=int, nargs="*", default=[3, 5, 10])
    parser.add_argument("--sequences", type=int, nargs="*", default=[1, 2, 3, 5])
    parser.add_argument("--cache", default="survey_cache/plastocyanin")
    parser.add_argument("--copper", action="store_true", help="only add (or refresh) the copper-site section for plastocyanin in the existing report")
    parser.add_argument("--out", default="docs/scope_checks")
    args = parser.parse_args(argv)
    if args.copper:
        append_section(args.out, copper_section(load_plastocyanin(args.cache), args.replicates, args.workers))
        return 0
    datasets = []
    if args.dataset in ("plastocyanin", "both"):
        datasets.append(load_plastocyanin(args.cache))
    controls = Path(os.environ.get("HCO_MANUSCRIPT", ROOT.parent / "CcO_Meta_Analysis_Manuscript/manuscript/controls"))
    if args.dataset in ("hco", "both"):
        if (controls / "output/auth_map.pkl").exists():
            datasets.append(load_hco(controls))
        else:
            print(f"skipping hco: no {controls}/output/auth_map.pkl")
    results = []
    for data in datasets:
        count = args.hco_replicates if (args.hco_replicates and data["groups"] != ["plastocyanin"]) else args.replicates
        real, candidates, summary = run_dataset(data, count, args.windows, args.sequences, args.workers)
        print(data["name"], len(real), "\n", summary.round(1).to_string(index=False))
        results.append((data, real, candidates, summary))
    if results:
        write_report(results, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
