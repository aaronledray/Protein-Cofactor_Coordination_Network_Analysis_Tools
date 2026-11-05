"""State-specific sequence votes for residue contacts (experimental, opt-in)."""

from __future__ import annotations

from typing import Optional

import pandas as pd

from . import conserved_interactions as ci


def state_structure_votes(atoms: pd.DataFrame, structures: pd.DataFrame, *, cutoff: float = 4.0,
                          min_separation: int = 2, structure_fraction: float = 0.5,
                          state_column: str = "state", sulfur_acceptor: bool = False,
                          metals: Optional[pd.DataFrame] = None, neighbor_contacts: bool = False) -> pd.DataFrame:
    """Return the structure counts behind every sequence-state contact vote.

    Rows cover each pooled candidate pair, caller-supplied state, group and sequence. A sequence-state with
    no structure in which both positions are mapped is retained with zero eligible structures and a missing
    vote; it is not a negative vote.
    """
    needed = {"structure_id", "group", "sequence_id", state_column}
    missing = needed - set(structures.columns)
    if missing:
        raise ValueError(f"structures lacks {sorted(missing)}")
    if structures["structure_id"].duplicated().any():
        raise ValueError("structures must contain one row per structure_id")
    if structures[state_column].isna().any() or structures[state_column].astype(str).str.strip().eq("").any():
        raise ValueError("every structure needs a non-empty state label")
    mismatch = set(atoms["structure_id"]) ^ set(structures["structure_id"])
    if mismatch:
        raise ValueError(f"atoms and structures must list the same structures; differing: {sorted(mismatch)[:4]}")

    keys = ["position_a", "position_b"]
    columns = [*keys, "state", "group", "sequence_id", "eligible_structures",
               "positive_structures", "structure_fraction", "sequence_vote"]
    rules = [
        {"table": "state_structure_votes", "name": "eligible_structures", "rule": "number of state-specific structures for this sequence where both contact positions are mapped", "cutoff_A": float("nan"), "status": "default"},
        {"table": "state_structure_votes", "name": "positive_structures", "rule": f"eligible structures in which the residue pair is within {cutoff:g} A", "cutoff_A": cutoff, "status": "default"},
        {"table": "state_structure_votes", "name": "structure_fraction", "rule": "positive_structures / eligible_structures; missing when no structures are eligible", "cutoff_A": float("nan"), "status": "default"},
        {"table": "state_structure_votes", "name": "sequence_vote", "rule": f"positive when at least {structure_fraction:.0%} of eligible structures show the contact; missing when none are eligible", "cutoff_A": float("nan"), "status": "default"},
    ]
    observations = ci.contact_observations(atoms, cutoff=cutoff, min_separation=min_separation,
                                           sulfur_acceptor=sulfur_acceptor, metals=metals,
                                           neighbor_contacts=neighbor_contacts)
    candidates = observations[keys].drop_duplicates().sort_values(keys)
    if candidates.empty:
        result = pd.DataFrame(columns=columns)
        result.attrs["rules"] = rules
        return result

    positions_by_structure = atoms.groupby("structure_id")["position"].agg(lambda values: set(values))
    positive_by_structure = {
        key: set(frame["structure_id"].unique())
        for key, frame in observations.loc[observations["contact"]].groupby(keys, sort=False)
    }
    states = sorted(structures[state_column].astype(str).unique())
    rows = []
    for pair in candidates.itertuples(index=False, name=None):
        pair_positive = positive_by_structure.get(pair, set())
        for state in states:
            state_structures = structures.loc[structures[state_column].astype(str) == state]
            for (group, sequence_id), sequence_structures in state_structures.groupby(["group", "sequence_id"], sort=True):
                eligible = []
                for structure_id in sequence_structures["structure_id"]:
                    present = positions_by_structure.get(structure_id, set())
                    if pair[0] in present and pair[1] in present:
                        eligible.append(structure_id)
                positive = pair_positive.intersection(eligible)
                n_eligible, n_positive = len(eligible), len(positive)
                fraction = n_positive / n_eligible if n_eligible else float("nan")
                vote = (n_positive / n_eligible >= structure_fraction) if n_eligible else pd.NA
                rows.append((*pair, state, group, sequence_id, n_eligible, n_positive, fraction, vote))
    result = pd.DataFrame(rows, columns=columns)
    result.attrs["rules"] = rules
    return result


def state_stratified_contacts(atoms: pd.DataFrame, structures: pd.DataFrame, *, cutoff: float = 4.0,
                              min_separation: int = 2, support: float = 0.90, min_sequences: int = 5,
                              structure_fraction: float = 0.5, state_column: str = "state", sulfur_acceptor: bool = False,
                              metals: Optional[pd.DataFrame] = None, neighbor_contacts: bool = False) -> pd.DataFrame:
    """Return one row per pooled candidate pair and state, with state-specific sequence votes.

    Structures are first evaluated together to define the union of contact candidates and the eligible
    (pair, structure) observations. Each state's observations are then voted independently using the
    same within-sequence majority and between-sequence support rules as :func:`ci.conserved_contacts`.
    A pair-state with no eligible sequences remains in the output with zero voters and fails the minimum.
    """
    needed = {"structure_id", "group", "sequence_id", state_column}
    missing = needed - set(structures.columns)
    if missing:
        raise ValueError(f"structures lacks {sorted(missing)}")
    if structures["structure_id"].duplicated().any():
        raise ValueError("structures must contain one row per structure_id")
    mismatch = set(atoms["structure_id"]) ^ set(structures["structure_id"])
    if mismatch:
        raise ValueError(f"atoms and structures must list the same structures; differing: {sorted(mismatch)[:4]}")
    if structures[state_column].isna().any() or structures[state_column].astype(str).str.strip().eq("").any():
        raise ValueError("every structure needs a non-empty state label")
    states = sorted(structures[state_column].astype(str).unique())
    if len(states) < 2:
        raise ValueError("state stratification needs at least two distinct state labels")

    keys = ["position_a", "position_b"]
    groups = sorted(structures["group"].unique())
    observations = ci.contact_observations(atoms, cutoff=cutoff, min_separation=min_separation, sulfur_acceptor=sulfur_acceptor,
                                           metals=metals, neighbor_contacts=neighbor_contacts)
    candidates = observations[keys].drop_duplicates().sort_values(keys).reset_index(drop=True)
    frames = []
    for state in states:
        state_structures = structures.loc[structures[state_column].astype(str) == state, ["structure_id", "group", "sequence_id"]]
        state_ids = set(state_structures["structure_id"])
        state_observations = observations[observations["structure_id"].isin(state_ids)]
        contact = ci.sequence_vote(state_observations, state_structures, keys, ["contact"], support=support,
                                   min_sequences=min_sequences, structure_fraction=structure_fraction)
        chemistry_flags = [*ci.CONTACT_CLASSES, *ci.CONTACT_FLAGS, *(ci.SULFUR_FLAGS if sulfur_acceptor else ())]
        chemistry = ci.sequence_vote(state_observations, state_structures, keys, chemistry_flags, condition="contact",
                                     support=support, min_sequences=min_sequences, structure_fraction=structure_fraction)
        frame = contact.merge(chemistry, on=keys, how="outer")
        votes = ci.sequence_votes(state_observations, state_structures, keys, "contact", structure_fraction=structure_fraction)
        for group in groups:
            column = f"contact_positive_{group}"
            if len(votes):
                positives = votes.loc[votes["group"] == group].groupby(keys)["vote"].sum().rename(column).reset_index()
            else:
                positives = pd.DataFrame(columns=[*keys, column])
            frame = frame.merge(positives, on=keys, how="left")
        if len(state_observations):
            coverage = ci.vote_coverage(votes, groups, support=support, min_sequences=min_sequences)
            frame = frame.merge(coverage, on=keys, how="left")
        frame["state"] = state
        frames.append(frame)

    if len(candidates):
        state_frame = pd.DataFrame({"state": states})
        candidates["_join"] = 1
        state_frame["_join"] = 1
        complete = candidates.merge(state_frame, on="_join").drop(columns="_join")
        result = complete.merge(pd.concat(frames, ignore_index=True), on=[*keys, "state"], how="left")
    else:
        result = pd.DataFrame(columns=[*keys, "state"])

    vote_flags = ["contact", *ci.CONTACT_CLASSES, *ci.CONTACT_FLAGS, *(ci.SULFUR_FLAGS if sulfur_acceptor else ())]
    for flag in vote_flags:
        for group in groups:
            count = f"{flag}_n_{group}"
            fraction = f"{flag}_{group}"
            conserved = f"conserved_{flag}"
            if count not in result:
                result[count] = 0
            else:
                result[count] = result[count].fillna(0).astype(int)
            if fraction not in result:
                result[fraction] = float("nan")
            if conserved not in result:
                result[conserved] = False
            else:
                result[conserved] = result[conserved].astype("boolean").fillna(False).astype(bool)

    for group in groups:
        positive = f"contact_positive_{group}"
        if positive not in result:
            result[positive] = 0
        else:
            result[positive] = result[positive].astype("Int64").fillna(0).astype(int)
        result[f"passes_contact_{group}"] = ((result[f"contact_n_{group}"] >= min_sequences)
                                             & (result[positive] / result[f"contact_n_{group}"].replace(0, float("nan")) >= support))

    for group in groups:
        for column, default in ((f"spare_negatives_{group}", 0), (f"spare_voters_{group}", -min_sequences)):
            if column not in result:
                result[column] = default
            else:
                result[column] = result[column].fillna(default).astype(int)
    for column, default in (("conserved_in", ""), ("held_back_by", ""), ("one_negative_away", ""),
                            ("failing_groups", "+".join(str(g) for g in groups))):
        if column not in result:
            result[column] = default
        else:
            result[column] = result[column].fillna(default).astype(str)

    ordered = [*keys, "state"]
    for flag in vote_flags:
        for group in groups:
            ordered.extend((f"{flag}_n_{group}", *([f"contact_positive_{group}"] if flag == "contact" else []), f"{flag}_{group}"))
            if flag == "contact":
                ordered.append(f"passes_contact_{group}")
        ordered.append(f"conserved_{flag}")
    ordered += [f"spare_negatives_{g}" for g in groups] + [f"spare_voters_{g}" for g in groups]
    ordered += ["conserved_in", "failing_groups", "held_back_by", "one_negative_away"]
    result = result[ordered].sort_values(["position_a", "position_b", "state"], kind="stable").reset_index(drop=True)
    rules = ci.class_rules("conserved_contacts", cutoff=cutoff, min_separation=min_separation, support=support,
                           min_sequences=min_sequences, structure_fraction=structure_fraction, coverage=True,
                           sulfur_acceptor=sulfur_acceptor, neighbor_contacts=neighbor_contacts).to_dict("records")
    for rule in rules:
        rule["table"] = "state_contacts"
    rules.append({"table": "state_contacts", "name": "state", "rule": "the vote is calculated independently within each caller-supplied state; structures are not independent observations", "cutoff_A": float("nan"), "status": "opt-in"})
    rules.append({"table": "state_contacts", "name": "contact_positive_<group>", "rule": "number of positive sequence votes in that state and group", "cutoff_A": float("nan"), "status": "default"})
    rules.append({"table": "state_contacts", "name": "passes_contact_<group>", "rule": f"at least {min_sequences} voting sequences and at least {support:.0%} positive sequence votes in that state and group", "cutoff_A": float("nan"), "status": "default"})
    result.attrs["rules"] = rules
    return result
