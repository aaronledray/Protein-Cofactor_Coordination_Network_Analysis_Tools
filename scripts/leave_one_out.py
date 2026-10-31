#!/usr/bin/env python3
"""Leave-one-sequence-out stability of the conserved-contact vote (scope check 3).

The vote asks whether at least ``support`` (90%) of the voting sequences of every group show a contact, with at least
``min_sequences`` (5) voting sequences per group. How much does the conserved set depend on any one sequence? Removing a sequence
changes only the counts of its group, not any geometry, so the per-sequence votes are computed once and the support rule is
re-evaluated 54 times for the heme-copper oxidase set:

* ``lost``: conserved contacts that stop passing when the sequence is left out. Leaving a sequence out can only raise the support
  of a group, so a contact is lost only because its group falls below the minimum number of voting sequences.
* ``gained``: contacts that start passing, because the sequence was the single negative that blocked them.
* a margin per conserved contact: how many more negative sequences it could take before it stops passing.

Two scenarios: the rule as it is, and the minimum relaxed by one for the group the sequence is left out of, which separates a real
dependence on a sequence from the group-size rule.

    python scripts/leave_one_out.py --dataset hco          # about two minutes, about 2.5 GB, one process
    python scripts/leave_one_out.py --dataset plastocyanin

Writes ``docs/scope_checks/08_leave_one_out.md`` and ``08_leave_one_out_sequences.csv``.
"""
import argparse
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from modules import conserved_interactions as ci  # noqa: E402

KEYS = ["position_a", "position_b"]


def count_tables(votes: pd.DataFrame, groups):
    """(voting sequences, positive votes) per key and group, as frames with one column per group."""
    grouped = votes.groupby(KEYS + ["group"], sort=False)["vote"]
    n = grouped.count().unstack("group").reindex(columns=list(groups)).fillna(0).astype(int)
    pos = grouped.sum().unstack("group").reindex(columns=list(groups)).fillna(0).astype(int)
    return n, pos


def passes(n: pd.DataFrame, pos: pd.DataFrame, support: float, minimum) -> pd.Series:
    """The support rule of :func:`sequence_vote`: in every group at least ``minimum`` voters and at least ``support`` of them positive.

    ``minimum`` is an int or a dict group -> int.
    """
    ok = pd.Series(True, index=n.index)
    for group in n.columns:
        floor = minimum[group] if isinstance(minimum, dict) else minimum
        ok &= (n[group] >= floor) & (pos[group] / n[group].where(n[group] > 0) >= support)
    return ok


def leave_one_out(votes: pd.DataFrame, groups, *, support: float = 0.90, min_sequences: int = 5) -> pd.DataFrame:
    """One row per sequence: contacts lost and gained when it is left out, with the minimum as it is and relaxed for its group."""
    n, pos = count_tables(votes, groups)
    baseline = passes(n, pos, support, min_sequences)
    by_sequence = {sequence: frame.set_index(KEYS)["vote"] for sequence, frame in votes.groupby("sequence_id")}
    group_of = votes.drop_duplicates("sequence_id").set_index("sequence_id")["group"]
    rows = []
    for sequence, vote in by_sequence.items():
        group = group_of[sequence]
        n_left, pos_left = n.copy(), pos.copy()
        voted = n_left.index.isin(vote.index)
        n_left.loc[voted, group] -= 1
        pos_left.loc[voted, group] -= vote.reindex(n_left.index[voted]).astype(int).to_numpy()
        row = {"sequence_id": sequence, "group": group, "voting_keys": len(vote)}
        for label, minimum in (("as_is", min_sequences), ("relaxed", {g: (min_sequences - 1 if g == group else min_sequences) for g in groups})):
            now = passes(n_left, pos_left, support, minimum)
            row[f"lost_{label}"] = int((baseline & ~now).sum())
            row[f"gained_{label}"] = int((~baseline & now).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def margins(votes: pd.DataFrame, groups, *, support: float = 0.90, min_sequences: int = 5) -> pd.DataFrame:
    """For every conserved key: the number of further negative votes each group could take, and the spare voters above the minimum.

    A group of ``n`` voters passes with at most ``floor(n * (1 - support))`` negatives, so one with n = 5 passes with none.
    """
    n, pos = count_tables(votes, groups)
    conserved = passes(n, pos, support, min_sequences)
    out = pd.DataFrame(index=n.index[conserved])
    for group in groups:
        allowed = np.floor(n[group] * (1 - support) + 1e-9).astype(int)
        out[f"spare_negatives_{group}"] = (allowed - (n[group] - pos[group]))[conserved]
        out[f"spare_voters_{group}"] = (n[group] - min_sequences)[conserved]
    out["margin"] = out[[c for c in out.columns if c.startswith("spare_negatives_")]].min(axis=1)
    return out


def label_of(structures: pd.DataFrame) -> dict:
    """A readable name per sequence id: its first structure id (the sequence ids of the oxidase set are 64-character hashes)."""
    first = structures.sort_values("structure_id").drop_duplicates("sequence_id").set_index("sequence_id")["structure_id"]
    return first.to_dict()


def report(name, votes, groups, structures, support=0.90, min_sequences=5):
    n, pos = count_tables(votes, groups)
    baseline = passes(n, pos, support, min_sequences)
    loo = leave_one_out(votes, groups, support=support, min_sequences=min_sequences)
    margin = margins(votes, groups, support=support, min_sequences=min_sequences)
    sizes = structures.groupby("group")["sequence_id"].nunique()
    lines = [f"## {name}", "",
             f"{int(baseline.sum())} conserved contacts among {len(n)} candidate pairs (support {support:.0%}, at least {min_sequences} voting sequences per group).", "",
             "| Group | sequences | negatives a contact can take and still pass (at 90%) |", "|---|---|---|"]
    for group in groups:
        lines.append(f"| {group} | {sizes[group]} | {math.floor(sizes[group] * (1 - support) + 1e-9)} (of {sizes[group]} if every sequence votes) |")
    lines += ["", "### Leaving each sequence out", "",
              "| Group | sequences | lost, rule as is (mean / max) | lost, minimum relaxed (max) | gained, minimum relaxed (mean / max) |", "|---|---|---|---|---|"]
    for group in groups:
        sub = loo[loo["group"] == group]
        lines.append(f"| {group} | {len(sub)} | {sub['lost_as_is'].mean():.1f} / {sub['lost_as_is'].max()} | {sub['lost_relaxed'].max()} | "
                     f"{sub['gained_relaxed'].mean():.1f} / {sub['gained_relaxed'].max()} |")
    names = label_of(structures)
    top = loo.sort_values("gained_relaxed", ascending=False).head(5)
    lines += ["", "Sequences whose removal frees the most contacts (the sequence was the single negative vote), named by their first structure: " +
              ", ".join(f"{names[r.sequence_id]} (group {r.group}, {r.gained_relaxed})" for r in top.itertuples()) + ".", ""]
    lines += ["### Margin of the conserved contacts", "",
              "Margin = the smallest, over groups, of the further negative votes a contact could take and still pass (0 means one more negative sequence "
              "in some group removes it).", "", "| Margin | conserved contacts | share |", "|---|---|---|"]
    counts = margin["margin"].value_counts().sort_index()
    for value, count in counts.items():
        lines.append(f"| {value} | {count} | {100 * count / len(margin):.1f}% |")
    for group in groups:
        zero = int((margin[f"spare_negatives_{group}"] == 0).sum())
        spare = int((margin[f"spare_voters_{group}"] == 0).sum())
        lines.append(f"\nGroup {group}: {zero} of {len(margin)} conserved contacts have no spare negative vote; {spare} sit exactly at the minimum of voters.")
    return "\n".join(lines) + "\n", loo, margin


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dataset", choices=["plastocyanin", "hco", "both"], default="both")
    parser.add_argument("--cache", default="survey_cache/plastocyanin")
    parser.add_argument("--out", default="docs/scope_checks")
    args = parser.parse_args(argv)
    import shuffle_control as sc

    sections, tables = ["# Leave-one-sequence-out stability of the conserved-contact vote", "",
                        "Generated by `scripts/leave_one_out.py` (scope check 3). The per-sequence contact votes are computed once; the support rule is then "
                        "re-evaluated with each sequence left out. Sequences are identified by the ids of the vote tables.", ""], []
    datasets = []
    if args.dataset in ("plastocyanin", "both"):
        datasets.append(sc.load_plastocyanin(args.cache))
    controls = Path(os.environ.get("HCO_MANUSCRIPT", ROOT.parent / "CcO_Meta_Analysis_Manuscript/manuscript/controls"))
    if args.dataset in ("hco", "both"):
        if (controls / "output/auth_map.pkl").exists():
            datasets.append(sc.load_hco(controls))
        else:
            print("skipping hco: manuscript position map not found")
    for data in datasets:
        observations = ci.contact_observations(data["atoms"])
        votes = ci.sequence_votes(observations, data["structures"], KEYS, "contact")
        del observations
        text, loo, margin = report(data["name"], votes, data["groups"], data["structures"])
        sections.append(text)
        tables.append(loo.assign(dataset=data["name"].split(" (")[0], first_structure=loo["sequence_id"].map(label_of(data["structures"]))))
        print(text)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "08_leave_one_out.md").write_text("\n".join(sections))
    pd.concat(tables, ignore_index=True).to_csv(out / "08_leave_one_out_sequences.csv", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
