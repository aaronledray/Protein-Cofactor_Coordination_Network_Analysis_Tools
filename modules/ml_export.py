"""Deterministic residue-level label export for downstream ML workflows.

This adapter is separate from the analysis core: it only reshapes the tables
returned by ``analyze_structure()`` and never changes legacy outputs. Each
exported row is one residue in one explicit structure/site/model example, so
labels are never pooled silently across models or cofactor copies.
"""

import re
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple, Union

import pandas as pd

from .cofactor_classes import (
    normalize_cofactor_class_config,
    resolve_cofactor_classes,
    resolve_effective_distance_cutoff,
)
from .assembly_policy import symmetry_context
from .coordination_api import analyze_structure

LABELING_CONTRACT_VERSION = "1.0"

LABEL_COLUMNS = [
    "contract_version", "structure_id", "site_id", "model_id",
    "cofactor_residue_name", "cofactor_residue_number", "cofactor_chain",
    "residue_name", "residue_number", "chain", "insertion_code", "hetero_flag",
    "shell", "shell_depth", "shells_present", "is_primary",
    "direct_coordination", "coordination_role", "motifs", "atoms_involved",
    "minimum_distance_A", "n_contacts",
]

SOLVENT_RESIDUES = frozenset({"HOH", "WAT", "DOD", "H2O"})

_RESIDUE_KEY = ["residue_name", "residue_number", "chain", "insertion_code", "hetero_flag"]
_ROLE_ORDER = ["primary_coordinator", "active_site_component", "network_context"]


def shell_depth(shell: str) -> int:
    """Map a shell label (PCS, SCS, TCS, ShellN) to its integer depth."""
    fixed = {"PCS": 1, "SCS": 2, "TCS": 3}
    if shell in fixed:
        return fixed[shell]
    match = re.fullmatch(r"Shell(\d+)", str(shell))
    if not match:
        raise ValueError(f"Unrecognized shell label: {shell!r}")
    return int(match.group(1))


def labeling_contract(**parameters: Any) -> Dict[str, Any]:
    """Return the contract record to store beside exported labels.

    ``parameters`` should be the keyword arguments passed to
    ``analyze_structure()`` (cutoff, shells, carbon-seed policy, ...). They are
    recorded verbatim so labels can be reproduced.
    """
    return {
        "contract_version": LABELING_CONTRACT_VERSION,
        "label_columns": list(LABEL_COLUMNS),
        "shell_semantics": (
            "Each residue is labeled with the shallowest shell it occupies; "
            "all occupied shells are listed in shells_present."
        ),
        "contacts_vs_links": (
            "contacts is the full atom-pair representation; links is the "
            "nearest-atom edge list kept for legacy compatibility."
        ),
        "analysis_parameters": dict(sorted(parameters.items())),
    }


class MixedCofactorError(ValueError):
    """Raised when one labeling example would mix cofactor families."""


# Defaults that define contract 1.0; overriding any is allowed but recorded.
CONTRACT_DEFAULTS: Dict[str, Any] = {
    "distance_cutoff": 3.6,
    "shells": 3,
    "include_carbon_seeds": False,
    "direct_coordination": True,
    "direct_coordination_cutoff": 2.6,
    "expand_residues": False,
    "site_mode": "per-site",
    "site_model_mode": "per-model",
}
_FORBIDDEN_OVERRIDES = {"site_mode", "first_model_only", "combinatorial_cofactor_cutoff"}


def analyze_labeling_example(
    structure_path,
    cofactor_resname: Union[str, Sequence[str]],
    *,
    model_id: Optional[Any] = None,
    site_id: Optional[str] = None,
    allow_mixed_cofactors: bool = False,
    exclude_solvent: bool = False,
    **options: Any,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Analyze one structure under the frozen labeling contract.

    Policy: every example is one explicit site and model (``per-site`` with
    ``per-model`` boundaries), and exactly one cofactor family per analysis.
    Mixed families would force a single shared cutoff (the largest configured
    one), so they are rejected unless ``allow_mixed_cofactors`` is set; the
    resulting contract then records ``mixed_cofactor_classes``. Run one
    analysis per family when per-family cutoffs matter.
    """
    bad = _FORBIDDEN_OVERRIDES & set(options)
    if bad:
        raise ValueError(f"Not overridable under the labeling contract: {sorted(bad)}")
    names = [cofactor_resname] if isinstance(cofactor_resname, str) else list(cofactor_resname)
    names = [n.strip() for part in names for n in str(part).split(",") if n.strip()]
    class_cutoffs = options.get("cofactor_class_cutoffs")
    classes = resolve_cofactor_classes(names, normalize_cofactor_class_config(class_cutoffs))
    if len(classes) > 1 and not allow_mixed_cofactors:
        raise MixedCofactorError(
            f"Cofactors {names} span families {list(classes)}; analyze each family "
            "separately or pass allow_mixed_cofactors=True."
        )
    params = {**CONTRACT_DEFAULTS, **options}
    tables = analyze_structure(structure_path, names, **params)
    labels = export_residue_labels(
        tables, model_id=model_id, site_id=site_id, exclude_solvent=exclude_solvent
    )
    contract = labeling_contract(**params)
    contract.update(
        cofactor_names=sorted(n.upper() for n in names),
        cofactor_classes=list(classes),
        effective_cutoff_A=resolve_effective_distance_cutoff(
            names, params["distance_cutoff"], class_cutoffs
        ),
        mixed_cofactor_classes=len(classes) > 1,
        symmetry_context=symmetry_context(structure_path),
    )
    if exclude_solvent:  # recorded only when set, so contract 1.0 output is unchanged
        contract["exclude_solvent"] = True
    return labels, contract


def _key(frame: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    return frame[list(columns)].astype(str).agg("\x1f".join, axis=1)


def export_residue_labels(
    tables: Mapping[str, pd.DataFrame],
    *,
    model_id: Optional[Any] = None,
    site_id: Optional[str] = None,
    exclude_solvent: bool = False,
) -> pd.DataFrame:
    """Emit one deterministic label row per structure/site/model/residue.

    ``model_id`` and ``site_id`` select a single training example. When the
    tables span several models or sites and no selector is given, the rows keep
    those identifiers as separate examples; a structure with multiple models
    in one pooled site raises instead, because pooled shell membership would
    then be attributed to each model ambiguously.
    """
    atoms = tables["atoms"]
    contacts = tables["contacts"]
    if atoms.empty:
        return pd.DataFrame(columns=LABEL_COLUMNS)

    if site_id is not None:
        atoms = atoms[atoms["site_id"] == site_id]
        contacts = contacts[contacts["site_id"] == site_id]
        if atoms.empty:
            raise ValueError(f"No atoms found for site_id {site_id!r}")

    models = sorted(atoms["model_id"].dropna().unique().tolist(), key=str)
    if model_id is not None:
        if model_id not in models:
            raise ValueError(f"model_id {model_id!r} not present; found {models}")
        atoms = atoms[atoms["model_id"] == model_id]
        contacts = contacts[
            (contacts["src_model_id"] == model_id) | (contacts["dst_model_id"] == model_id)
        ]
        models = [model_id]
    elif len(models) > 1 and atoms["site_id"].eq("all").any():
        raise ValueError(
            f"Structure has multiple models {models} pooled into one site; "
            "pass model_id, use first_model_only=True, or site_model_mode='per-model'."
        )

    cofactors = (
        atoms[atoms["shell"] == "Cofactor"]
        .groupby(["structure_id", "site_id"], sort=True)
        .agg(
            cofactor_residue_name=("residue_name", lambda v: ";".join(sorted(set(map(str, v))))),
            cofactor_residue_number=("residue_number", lambda v: ";".join(sorted(set(map(str, v))))),
            cofactor_chain=("chain", lambda v: ";".join(sorted(set(map(str, v))))),
        )
        .reset_index()
    )

    shell_atoms = atoms[atoms["shell"] != "Cofactor"].copy()
    if exclude_solvent:
        # Only the exported rows are filtered; waters still take part in shell
        # propagation, so water-mediated residues keep their shell labels.
        shell_atoms = shell_atoms[~shell_atoms["residue_name"].isin(SOLVENT_RESIDUES)]
    if shell_atoms.empty:
        return pd.DataFrame(columns=LABEL_COLUMNS)
    shell_atoms["shell_depth"] = shell_atoms["shell"].map(shell_depth)

    group_cols = ["structure_id", "site_id", "model_id"] + _RESIDUE_KEY
    rows = []
    for key, group in shell_atoms.groupby(group_cols, sort=True, dropna=False):
        depths = sorted(group["shell_depth"].unique())
        top = depths[0]
        top_group = group[group["shell_depth"] == top]
        roles = set(top_group["coordination_role"])
        role = next((r for r in _ROLE_ORDER if r in roles), "")
        motifs = ",".join(sorted({str(m) for m in top_group["motif"] if str(m)}))
        rows.append(
            dict(
                zip(group_cols, key),
                shell=top_group["shell"].iloc[0],
                shell_depth=int(top),
                shells_present=",".join(
                    group.loc[group["shell_depth"] == d, "shell"].iloc[0] for d in depths
                ),
                is_primary=bool(top == 1),
                coordination_role=role,
                motifs=motifs,
                atoms_involved=",".join(sorted(set(map(str, top_group["atom_name"])))),
            )
        )
    labels = pd.DataFrame(rows)

    # Contact-derived evidence, keyed on the destination residue.
    if contacts.empty:
        labels["minimum_distance_A"] = float("nan")
        labels["n_contacts"] = 0
        labels["direct_coordination"] = False
    else:
        dst = contacts.rename(
            columns={
                "dst_resname": "residue_name", "dst_resnum": "residue_number",
                "dst_chain": "chain", "dst_insertion_code": "insertion_code",
                "dst_hetero_flag": "hetero_flag", "dst_model_id": "model_id",
            }
        )
        dst_shell = dst["dst_shell"].map(shell_depth)
        # Use contacts that enter the residue's own (shallowest) shell.
        stats = (
            dst.assign(_depth=dst_shell)
            .groupby(["structure_id", "site_id", "model_id"] + _RESIDUE_KEY + ["_depth"], dropna=False)
            .agg(
                minimum_distance_A=("distance_A", "min"),
                n_contacts=("distance_A", "size"),
                direct_coordination=("direct_coordination", "any"),
            )
            .reset_index()
            .rename(columns={"_depth": "shell_depth"})
        )
        merge_cols = ["structure_id", "site_id", "model_id"] + _RESIDUE_KEY + ["shell_depth"]
        for frame in (labels, stats):
            for col in _RESIDUE_KEY:
                frame[col] = frame[col].astype(str)
            frame["model_id"] = frame["model_id"].astype(str)
        labels = labels.merge(stats, on=merge_cols, how="left")
        labels["n_contacts"] = labels["n_contacts"].fillna(0).astype(int)
        labels["direct_coordination"] = labels["direct_coordination"].fillna(False).astype(bool)

    labels = labels.merge(cofactors, on=["structure_id", "site_id"], how="left")
    labels["contract_version"] = LABELING_CONTRACT_VERSION
    labels["minimum_distance_A"] = labels["minimum_distance_A"].round(4)
    labels = labels[LABEL_COLUMNS]
    sort_cols = ["structure_id", "site_id", "model_id", "shell_depth", "chain",
                 "residue_number", "insertion_code", "residue_name"]
    return labels.sort_values(sort_cols, kind="mergesort", key=_numeric_aware).reset_index(drop=True)


def _numeric_aware(column: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(column, errors="coerce")
    return numeric if numeric.notna().all() else column.astype(str)
