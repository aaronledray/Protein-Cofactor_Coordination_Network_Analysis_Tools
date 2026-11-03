"""Contacts observed across the structures of one protein (experimental; ``sscna conserve --single-protein``).

The sequence vote of :mod:`modules.conserved_interactions` needs several sequences: with one sequence there is one vote, and the minimum
of voting sequences per group can never be met. What a single protein offers is its structures (several entries, chains or states).
They are not independent (the same crystal form, laboratory and refinement protocol), so nothing here is evidence of conservation. The
tables report what is seen: in how many of the structures a contact is present, and how close the residues come in each.

``contact_persistence`` has one row per residue pair in contact in at least one structure; ``cofactor_persistence`` the same for the
residues recorded around a cofactor feature. The distance is part of the evidence because within one sequence most of the variation
between structures sits at the contact cutoff: a pair whose median closest distance lies within ``near_margin`` of the cutoff is
flagged ``near_cutoff``, since its presence depends on where the cutoff falls.
"""

from __future__ import annotations

import re
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd

from . import conserved_interactions as ci

NEAR_MARGIN = 0.3            # A: a median closest distance within this of the cutoff is flagged near_cutoff
DISTANCE_REACH = 2.0         # A beyond the cutoff up to which a closest distance is measured; farther is recorded as the limit
KEYS = ["position_a", "position_b"]
CHEMISTRY = (*ci.CONTACT_CLASSES, *ci.CONTACT_FLAGS)


def require_one_sequence(structures: pd.DataFrame) -> None:
    """Raise ``ValueError`` unless every structure belongs to one sequence (the evidence of this mode is not the sequence vote)."""
    found = sorted(structures["sequence_id"].unique())
    if len(found) != 1:
        raise ValueError(f"single-protein mode needs exactly one sequence_id, found {len(found)}: {found[:4]}{' ...' if len(found) > 4 else ''}; "
                         "with several sequences use the sequence vote (without --single-protein)")


def states_of(structures: pd.DataFrame) -> List[str]:
    """The labels of the optional ``state`` column (empty without it); a label must be letters, digits and underscores."""
    if "state" not in structures.columns:
        return []
    labels = sorted(str(s) for s in structures["state"].dropna().unique())
    bad = [s for s in labels if not re.fullmatch(r"[A-Za-z0-9_]+", s)]
    if bad:
        raise ValueError(f"state labels must be letters, digits and underscores, not {bad[:3]}")
    return labels


def _state_columns(frame: pd.DataFrame, structures: pd.DataFrame, keys: Sequence[str], present: str, columns: dict) -> pd.DataFrame:
    """Per state, ``<present>_<state>`` (structures with a row) and the structures showing each flag; ``columns`` maps flag -> column prefix."""
    states = states_of(structures)
    if len(states) < 2:
        return pd.DataFrame()
    labelled = frame.merge(structures[["structure_id", "state"]].drop_duplicates("structure_id"), on="structure_id", how="left")
    out = {}
    for state in states:
        sub = labelled[labelled["state"].astype(str) == state].groupby(list(keys))
        out[f"{present}_{state}"] = sub.size()
        for flag, prefix in columns.items():
            out[f"{prefix}_{state}"] = sub[flag].sum().astype(int)
    return pd.DataFrame(out).fillna(0).astype(int)


def _same_structures(atoms: pd.DataFrame, structures: pd.DataFrame) -> None:
    extra = set(atoms["structure_id"]) ^ set(structures["structure_id"])
    if extra:
        raise ValueError(f"atoms and structures must list the same structures; differing: {sorted(extra)[:4]}")


def contact_persistence(atoms: pd.DataFrame, structures: pd.DataFrame, *, cutoff: float = 4.0, min_separation: int = 2,
                        near_margin: float = NEAR_MARGIN) -> pd.DataFrame:
    """One row per residue pair in contact (rules of :func:`residue_contacts`) in at least one structure, with the evidence of this protein.

    ``structures_present``: structures that map both positions; ``structures_in_contact``: of those, structures where the pair is in contact;
    ``fraction_in_contact``; ``closest_min`` / ``closest_median`` / ``closest_max``: the closest heavy-atom distance of the pair over the
    structures that map both (a distance beyond ``cutoff + 2 A`` is recorded as that limit); ``near_cutoff``: the median is within
    ``near_margin`` of ``cutoff``; ``chemistry``: the contact classes and flags present in every structure where the pair is in contact. With
    two or more states in ``structures["state"]``, ``present_<state>`` and ``in_contact_<state>`` count the same per state.
    """
    require_one_sequence(structures)
    states_of(structures)
    _same_structures(atoms, structures)
    limit = cutoff + DISTANCE_REACH
    obs = ci.contact_observations(atoms, cutoff=cutoff, min_separation=min_separation)
    wide = ci.residue_contacts(atoms, cutoff=limit, min_separation=min_separation, chemistry=False)[["structure_id", *KEYS, "min_distance"]]
    obs = obs.drop(columns=["min_distance"]).merge(wide, on=["structure_id", *KEYS], how="left")
    obs["distance"] = obs["min_distance"].fillna(limit).clip(upper=limit)
    grouped = obs.groupby(KEYS, sort=True)
    table = pd.DataFrame({"structures_present": grouped.size(), "structures_in_contact": grouped["contact"].sum().astype(int)})
    table["fraction_in_contact"] = (table["structures_in_contact"] / table["structures_present"]).round(3)
    table["closest_min"] = grouped["distance"].min().round(3)
    table["closest_median"] = grouped["distance"].median().round(3)
    table["closest_max"] = grouped["distance"].max().round(3)
    table["near_cutoff"] = (table["closest_median"] - cutoff).abs() <= near_margin + 1e-9
    in_contact = obs[obs["contact"]]
    named = [c for c in CHEMISTRY if c in in_contact.columns]
    every = in_contact.groupby(KEYS)[named].agg(lambda s: bool(s.eq(True).all()))
    table["chemistry"] = every.apply(lambda row: ";".join(c for c in named if row[c]), axis=1).reindex(table.index).fillna("")
    names = in_contact.groupby(KEYS)[["resname_a", "resname_b"]].first()
    table = names.join(table, how="right").reset_index()
    table = table[table["structures_in_contact"] > 0]
    extra = _state_columns(obs, structures, KEYS, "present", {"contact": "in_contact"})
    if len(extra.columns):
        table = table.merge(extra.reset_index(), on=KEYS, how="left")
    table = table.reset_index(drop=True)
    rules = ci.class_rules("residue_contacts", cutoff=cutoff, min_separation=min_separation).to_dict("records")
    rules += [
        _rule("contact_persistence", "structures_present", "structures of this protein that map both positions (a structure that lacks one is not counted as a negative)"),
        _rule("contact_persistence", "structures_in_contact", f"of those, structures where the closest heavy-atom pair is within {cutoff:g} A (the rule of `contact`)", cutoff),
        _rule("contact_persistence", "closest_*", f"the closest heavy-atom distance of the pair in each structure that maps both; min, median and max over those structures; a distance "
              f"beyond {limit:g} A (cutoff + {DISTANCE_REACH:g}) is recorded as {limit:g}"),
        _rule("contact_persistence", "near_cutoff", f"the median closest distance is within {near_margin:g} A of the {cutoff:g} A cutoff, so the contact's presence depends on where the cutoff falls",
              near_margin),
        _rule("contact_persistence", "chemistry", "the contact classes and flags that are true in every structure where the pair is in contact"),
        _rule("contact_persistence", "present_<state> / in_contact_<state>", "the same counts within the structures labelled with this state (only with two or more states)", np.nan, "opt-in"),
    ]
    table.attrs["rules"] = rules
    return table


def _rule(table: str, name: str, rule: str, cutoff: float = np.nan, status: str = "default") -> dict:
    return {"table": table, "name": name, "rule": rule, "cutoff_A": cutoff, "status": status}


def cofactor_persistence(cofactor: pd.DataFrame, structures: pd.DataFrame, flags: Sequence[str]) -> pd.DataFrame:
    """One row per (feature, position) in which a class in ``flags`` is true in at least one structure, from :func:`cofactor_contacts`.

    ``structures_recorded``: structures in which the residue is recorded around the feature (within its ``record_distance``);
    ``<flag>_structures``: of those, structures where the class holds; ``distance_min`` / ``distance_median`` / ``distance_max``: the
    closest distance of the residue to the feature over those structures. With two or more states, ``recorded_<state>`` and
    ``<flag>_<state>`` count the same per state. Distances to cofactor classes use several cutoffs, so there is no ``near_cutoff`` here.
    """
    require_one_sequence(structures)
    states_of(structures)
    keys = ["feature", "position"]
    grouped = cofactor.groupby(keys, sort=True)
    table = pd.DataFrame({"resname": grouped["resname"].first(), "structures_recorded": grouped.size()})
    for flag in flags:
        table[f"{flag}_structures"] = grouped[flag].sum().astype(int)
    table["distance_min"] = grouped["min_distance"].min().round(3)
    table["distance_median"] = grouped["min_distance"].median().round(3)
    table["distance_max"] = grouped["min_distance"].max().round(3)
    table = table[(table[[f"{flag}_structures" for flag in flags]] > 0).any(axis=1)].reset_index()
    extra = _state_columns(cofactor, structures, keys, "recorded", {flag: flag for flag in flags})
    if len(extra.columns):
        table = table.merge(extra.reset_index(), on=keys, how="left")
    rules = ci.rules_of(cofactor).to_dict("records") + [
        _rule("cofactor_persistence", "structures_recorded", "structures in which the residue lies within the record distance of the feature (a structure without a row is not counted)"),
        _rule("cofactor_persistence", "<flag>_structures", "of those, structures where the class of `cofactor_contacts` of that name holds"),
        _rule("cofactor_persistence", "distance_*", "the closest distance of the residue to the feature; min, median and max over the recorded structures"),
        _rule("cofactor_persistence", "recorded_<state> / <flag>_<state>", "the same counts within the structures labelled with this state (only with two or more states)", np.nan, "opt-in"),
    ]
    table.attrs["rules"] = rules
    return table
