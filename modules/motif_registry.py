"""Residue-specific motif names and atom membership.

The original chemistry table is still the source of the broad legacy
vocabulary.  This module adds a small, explicit layer on top of it so that
motifs used by the canonical API, comparison code, and interactive viewer
have one resolver and one description of their atom membership.

The public labels intentionally retain compatibility-sensitive names such as
``COO`` and ``heme_por``.  Human-friendly aliases are accepted by
``motif_atoms`` for reports and downstream code.
"""

from __future__ import annotations

from types import MappingProxyType
from typing import Dict, FrozenSet, Mapping, Optional, Set, Tuple

from .chemistry import chemical_moieties


RESIDUE_NAME_ALIASES: Mapping[str, str] = MappingProxyType(
    {
        "HID": "HIS",
        "HIE": "HIS",
        "HIP": "HIS",
    }
)

HEME_PROPIONATE_OXYGENS: FrozenSet[str] = frozenset(
    {"O1A", "O2A", "O1D", "O2D"}
)

_MOTIF_ALIASES: Mapping[str, str] = MappingProxyType(
    {
        "CARBOXYLATE": "COO",
        "COO": "COO",
        "HEME": "heme_por",
        "HEME_PORPHYRIN": "heme_por",
        "HEME_POR": "heme_por",
        "THIOL": "thiol",
        "IMIDAZOLE": "imidazole",
    }
)

# Explicit atom membership for motifs where residue identity matters.  The
# values are merged with the legacy chemistry table below, rather than
# replacing that table wholesale.
_EXPLICIT_MOTIF_ATOMS: Mapping[Tuple[str, str], FrozenSet[str]] = {
    # Histidine protonation variants all share the same imidazole ring.  Ring
    # hydrogens are included when a deposited structure contains them.
    ("HIS", "imidazole"): frozenset(
        {"CG", "ND1", "CD2", "CE1", "NE2", "HD1", "HD2", "HE1", "HE2"}
    ),
    # The carbonyl carbon belongs to the carboxylate motif as well as both O
    # atoms; this preserves the observed OE1/OE2 or OD1/OD2 identity.
    ("ASP", "COO"): frozenset({"CG", "OD1", "OD2"}),
    ("GLU", "COO"): frozenset({"CD", "OE1", "OE2"}),
    ("CYS", "thiol"): frozenset({"SG", "HG"}),
    # Heme propionate oxygens are kept distinct from the porphyrin macrocycle
    # because they can participate in inferred hydrogen-bond contacts.
    ("HEM", "heme_propionate"): HEME_PROPIONATE_OXYGENS,
    # Deposited FeMo-cofactor and homocitrate names vary slightly between
    # dictionaries; include the names present in the reference nitrogenase.
    ("ICS", "femoco_mo"): frozenset({"MO", "MO1"}),
    ("HCA", "hca_c"): frozenset({*(f"C{i}" for i in range(1, 8))}),
    ("HCA", "hca_o"): frozenset({*(f"O{i}" for i in range(1, 8))}),
}

_IRON_SULFUR_RESIDUES = frozenset({"SF4", "FES", "FS4", "CLF", "BCLF"})

_IRON_SULFUR_ATOMS: Mapping[str, Mapping[str, FrozenSet[str]]] = {
    residue: {
        "iron_sulfur_fe": frozenset({*(f"FE{i}" for i in range(1, 5))}),
        "iron_sulfur_s": frozenset({*(f"S{i}" for i in range(1, 5))}),
    }
    for residue in ("SF4", "FES", "FS4")
}
_IRON_SULFUR_ATOMS = {
    **_IRON_SULFUR_ATOMS,
    "CLF": {
        "iron_sulfur_fe": frozenset({*(f"FE{i}" for i in range(1, 9))}),
        "iron_sulfur_s": frozenset({"S1", "S2A", "S3A", "S4A", "S2B", "S3B", "S4B"}),
    },
    "BCLF": {
        "iron_sulfur_fe": frozenset({*(f"FE{i}" for i in range(1, 9))}),
        "iron_sulfur_s": frozenset({"S1", "S2A", "S3A", "S4A", "S2B", "S3B", "S4B"}),
    },
}


def _normalized_residue(residue_name: object) -> str:
    name = str(residue_name or "").strip().upper()
    return RESIDUE_NAME_ALIASES.get(name, name)


def _normalized_motif(motif: object) -> str:
    raw = str(motif or "").strip()
    return _MOTIF_ALIASES.get(raw.upper(), raw)


def _build_membership() -> Dict[str, Dict[str, Set[str]]]:
    membership: Dict[str, Dict[str, Set[str]]] = {}
    for (residue_name, atom_name), motif in chemical_moieties.items():
        residue = str(residue_name).strip().upper()
        atom = str(atom_name).strip().upper()
        if residue == "ANY":
            continue
        membership.setdefault(residue, {}).setdefault(str(motif), set()).add(atom)

    # Make aliases visible in the public registry, not only in the resolver.
    for alias, base in RESIDUE_NAME_ALIASES.items():
        for motif, atoms in membership.get(base, {}).items():
            membership.setdefault(alias, {}).setdefault(motif, set()).update(atoms)

    for (residue, motif), atoms in _EXPLICIT_MOTIF_ATOMS.items():
        membership.setdefault(residue, {}).setdefault(motif, set()).update(atoms)

    # The explicit HEM propionate oxygen classification supersedes the broad
    # legacy ``heme_por`` label for those four coordinating atoms.
    porphyrin_atoms = membership.get("HEM", {}).get("heme_por", set())
    porphyrin_atoms.difference_update(HEME_PROPIONATE_OXYGENS)

    # The common Fe-S cluster residue names are not all present in the old
    # table.  Add their structural classes here so they do not become
    # ``unknown_motif`` in canonical network signatures.
    for residue, motifs in _IRON_SULFUR_ATOMS.items():
        names = membership.setdefault(residue, {})
        for motif, atoms in motifs.items():
            names.setdefault(motif, set()).update(atoms)

    return membership


_MOTIF_MEMBERSHIP = _build_membership()

# Public read-only nested mapping.  A tuple-key mapping would be compact, but
# residue -> motif -> atoms is easier to inspect and consume in applications.
MOTIF_ATOM_MEMBERSHIP: Mapping[str, Mapping[str, FrozenSet[str]]] = MappingProxyType(
    {
        residue: MappingProxyType(
            {motif: frozenset(atoms) for motif, atoms in motifs.items()}
        )
        for residue, motifs in _MOTIF_MEMBERSHIP.items()
    }
)


def _dynamic_iron_sulfur_motif(residue: str, atom: str) -> Optional[str]:
    if residue not in _IRON_SULFUR_RESIDUES:
        return None
    if atom in _MOTIF_MEMBERSHIP.get(residue, {}).get("iron_sulfur_fe", set()):
        return "iron_sulfur_fe"
    if atom in _MOTIF_MEMBERSHIP.get(residue, {}).get("iron_sulfur_s", set()):
        return "iron_sulfur_s"
    return None


def motif_for_atom(residue_name: object, atom_name: object) -> str:
    """Return the canonical motif label for one residue atom.

    Resolution order is explicit residue-specific vocabulary, protonation
    aliases, cluster naming rules, the legacy chemistry table, and finally
    the legacy ``ANY`` atom fallback.  The final value is always a string so
    it can be written directly into a pandas table or HTML hover label.
    """
    residue = _normalized_residue(residue_name)
    atom = str(atom_name or "").strip().upper()

    if residue == "HEM" and atom in HEME_PROPIONATE_OXYGENS:
        return "heme_propionate"

    explicit = _EXPLICIT_MOTIF_ATOMS
    for candidate in (residue, str(residue_name or "").strip().upper()):
        for (candidate_residue, motif), atoms in explicit.items():
            if candidate_residue == candidate and atom in atoms:
                return motif

    cluster_motif = _dynamic_iron_sulfur_motif(residue, atom)
    if cluster_motif is not None:
        return cluster_motif

    label = chemical_moieties.get((residue, atom))
    if label is None:
        raw_residue = str(residue_name or "").strip().upper()
        label = chemical_moieties.get((raw_residue, atom))
    if label is None:
        label = chemical_moieties.get(("ANY", atom))
    return str(label or "unknown_motif")


def motif_atoms(residue_name: object, motif: object) -> FrozenSet[str]:
    """Return atom names belonging to ``motif`` in a residue.

    Unknown residue/motif pairs return an empty set.  Histidine aliases and
    friendly motif names such as ``carboxylate`` are accepted.
    """
    residue = _normalized_residue(residue_name)
    canonical_motif = _normalized_motif(motif)
    direct = MOTIF_ATOM_MEMBERSHIP.get(residue, {}).get(canonical_motif)
    if direct is not None:
        return direct
    raw_residue = str(residue_name or "").strip().upper()
    return MOTIF_ATOM_MEMBERSHIP.get(raw_residue, {}).get(canonical_motif, frozenset())


def motif_vocabulary() -> Mapping[str, Mapping[str, FrozenSet[str]]]:
    """Return the read-only residue-specific motif vocabulary."""
    return MOTIF_ATOM_MEMBERSHIP


__all__ = [
    "HEME_PROPIONATE_OXYGENS",
    "MOTIF_ATOM_MEMBERSHIP",
    "RESIDUE_NAME_ALIASES",
    "motif_atoms",
    "motif_for_atom",
    "motif_vocabulary",
]
