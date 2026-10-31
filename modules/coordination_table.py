"""Per-metal, per-donor-class coordination cutoffs (experimental, opt-in).

The analysis defaults use one global cutoff for every metal-donor pair (``cofactor_contacts`` 2.8 A,
``analyze_structure`` direct coordination 2.6 A). Bond lengths differ by metal and donor, and some ligands
are real but weak (the axial methionine of a blue copper site), so this table keys a cutoff pair on
(metal element, donor class):

- ``bond_max``: a metal-donor distance at or below this is a ``bond``.
- ``weak_max``: beyond ``bond_max`` and at or below this is ``weak`` (None when the class has no weak range).

A pair that is not in the table is ``default``: it follows the caller's global cutoff, so unlisted
metals and donors behave as before and say so. Backbone amide N is not a donor class here: it stays a
contact and a shell member but is never tagged ``bond`` or ``weak``.

The table grows one metal at a time, from evidence a scope check needs. Every row carries its source.
Current rows: type 1 copper, from the PDB survey in docs/scope_checks/02_cu_type1_survey.md (87
proteins, X-ray, 1.8 A or better). The Cu-S(Met) tiers have no natural break in that survey (axial Met
runs 2.46 to 3.57 A across families); 2.7 and 3.3 A are conventions.
"""

from typing import NamedTuple, Optional, Tuple

import pandas as pd

from .coordination_api import METAL_ION_RESNAMES

TIERS = ("bond", "weak", "default", "none")
_RANK = {"none": 0, "default": 1, "weak": 2, "bond": 3}


class Rule(NamedTuple):
    metal: str                      # element symbol, upper case
    donor_class: str
    donors: Tuple[Tuple[str, str], ...]      # (residue name, atom name)
    bond_max: float
    weak_max: Optional[float]
    source: str


class Call(NamedTuple):
    tier: str                       # one of TIERS
    donor_class: str                # "" when the pair is not in the table
    rule: str                       # human-readable statement of the cutoff that decided the tier


SURVEY = "PDB survey of type 1 Cu sites (docs/scope_checks/02_cu_type1_survey.md)"
RULES: Tuple[Rule, ...] = (
    Rule("CU", "His N", (("HIS", "ND1"), ("HIS", "NE2")), 2.4, None, SURVEY + ": 87 proteins, max 2.37 A"),
    Rule("CU", "Cys S", (("CYS", "SG"),), 2.5, None, SURVEY + ": 87 proteins, max 2.42 A"),
    Rule("CU", "Met S", (("MET", "SD"), ("MSE", "SE")), 2.7, 3.3,
         SURVEY + ": 55 proteins with an axial Met, 2.46 to 3.57 A, no natural break (convention)"),
)


def find_rule(metal: str, resname: str, atom: str) -> Optional[Rule]:
    key = (resname.strip().upper(), atom.strip().upper())
    for rule in RULES:
        if rule.metal == metal.strip().upper() and key in rule.donors:
            return rule
    return None


def coordination_tier(metal: str, resname: str, atom: str, distance: float, *, default_max: float = 2.8) -> Call:
    """Tier of one metal-donor atom pair; ``default_max`` is the global cutoff used for unlisted pairs."""
    rule = find_rule(metal, resname, atom)
    if rule is None:
        if distance <= default_max:
            return Call("default", "", f"not in table; global <= {default_max:g}")
        return Call("none", "", f"not in table; beyond global {default_max:g}")
    text = f"{rule.metal.title()}-{rule.donor_class}: bond <= {rule.bond_max:g}"
    if rule.weak_max is not None:
        text += f", weak <= {rule.weak_max:g}"
    if distance <= rule.bond_max:
        return Call("bond", rule.donor_class, text)
    if rule.weak_max is not None and distance <= rule.weak_max:
        return Call("weak", rule.donor_class, text)
    return Call("none", rule.donor_class, text)


def strongest(calls) -> Call:
    """The best call among (distance, Call) pairs: highest tier, then the shorter distance."""
    best = None
    for distance, call in calls:
        if best is None or (_RANK[call.tier], -distance) > (_RANK[best[1].tier], -best[0]):
            best = (distance, call)
    return best[1] if best else Call("none", "", "")


TIER_COLUMNS = ("coord_tier", "coord_donor_class", "coord_rule")


def annotate_contacts(contacts: pd.DataFrame, *, default_max: float = 2.6) -> pd.DataFrame:
    """Add ``coord_tier``, ``coord_donor_class`` and ``coord_rule`` to the ``contacts`` table of ``analyze_structure``.

    Only cofactor-to-primary-shell rows whose source atom is a metal are judged; the other rows get empty
    strings. ``default_max`` is the cutoff for pairs not in the table (the 2.6 A of ``direct_coordination``).
    The input is not modified and the existing columns, including ``direct_coordination``, are untouched.
    """
    out = contacts.copy()
    tiers = [("", "", "")] * len(out)
    if len(out):
        judged = (out["link_type"] == "cofactor->pcs") & out["src_element"].astype(str).str.upper().isin(METAL_ION_RESNAMES)
        for n, row in enumerate(out.itertuples(index=False)):
            if judged.iloc[n]:
                tiers[n] = coordination_tier(str(row.src_element), str(row.dst_resname), str(row.dst_atom),
                                             float(row.distance_A), default_max=default_max)
    out[list(TIER_COLUMNS)] = pd.DataFrame(tiers, index=out.index, columns=list(TIER_COLUMNS))
    return out
