"""Policy-sensitivity comparison of exported residue labels.

Compares label tables produced under alternative analysis parameters against a
baseline so a labeling contract can be frozen with known stability. This module
only reads ``export_residue_labels`` output and never touches legacy outputs.
"""

from typing import Any, Dict, Mapping, Optional, Sequence

import pandas as pd

from .coordination_api import analyze_structure
from .ml_export import export_residue_labels

_KEY = ["site_id", "model_id", "chain", "residue_name", "residue_number", "insertion_code"]

DEFAULT_VARIANTS: Dict[str, Dict[str, Any]] = {
    "cutoff_3.2": {"distance_cutoff": 3.2},
    "cutoff_3.4": {"distance_cutoff": 3.4},
    "cutoff_3.8": {"distance_cutoff": 3.8},
    "cutoff_4.0": {"distance_cutoff": 4.0},
    "carbon_seeds": {"include_carbon_seeds": True},
    "expand_residues": {"expand_residues": True},
    "exclude_ala_sidechain": {"exclude_moieties": ["alanine_sidechain"]},
}

SUMMARY_COLUMNS = [
    "variant", "n_labels", "n_primary", "primary_jaccard", "labeled_jaccard",
    "shell_changed", "gained", "lost",
]


def _keys(frame: pd.DataFrame) -> pd.Series:
    if frame.empty:
        return pd.Series(dtype=str)
    return frame[_KEY].astype(str).agg("|".join, axis=1)


def _jaccard(left: set, right: set) -> float:
    union = left | right
    return 1.0 if not union else len(left & right) / len(union)


def compare_labels(baseline: pd.DataFrame, variant: pd.DataFrame) -> Dict[str, Any]:
    """Summarize how a variant's shell labels differ from the baseline's."""
    base = baseline.assign(_k=_keys(baseline))
    var = variant.assign(_k=_keys(variant))
    base_primary = set(base.loc[base["is_primary"].astype(bool), "_k"]) if len(base) else set()
    var_primary = set(var.loc[var["is_primary"].astype(bool), "_k"]) if len(var) else set()
    base_depth = dict(zip(base["_k"], base["shell_depth"]))
    var_depth = dict(zip(var["_k"], var["shell_depth"]))
    shared = set(base_depth) & set(var_depth)
    return {
        "n_labels": len(variant),
        "n_primary": len(var_primary),
        "primary_jaccard": round(_jaccard(base_primary, var_primary), 3),
        "labeled_jaccard": round(_jaccard(set(base_depth), set(var_depth)), 3),
        "shell_changed": sum(base_depth[k] != var_depth[k] for k in shared),
        "gained": len(set(var_depth) - set(base_depth)),
        "lost": len(set(base_depth) - set(var_depth)),
    }


def sensitivity_table(
    structure_path,
    cofactor_resname,
    *,
    baseline_options: Optional[Mapping[str, Any]] = None,
    variants: Optional[Mapping[str, Mapping[str, Any]]] = None,
    shells: int = 3,
    label_selector: Optional[Mapping[str, Any]] = None,
) -> pd.DataFrame:
    """Run baseline plus variants and return one summary row per variant.

    A variant that raises is reported with ``error`` rather than aborting.
    Only keyword options accepted by ``analyze_structure`` are valid.
    """
    base_options = {"shells": shells, "site_mode": "per-site", **(baseline_options or {})}
    selector = dict(label_selector or {})

    def labels_for(options: Mapping[str, Any]) -> pd.DataFrame:
        return export_residue_labels(
            analyze_structure(structure_path, cofactor_resname, **options), **selector
        )

    baseline = labels_for(base_options)
    rows = [{"variant": "baseline", **compare_labels(baseline, baseline)}]
    for name, overrides in (variants or DEFAULT_VARIANTS).items():
        try:
            row = compare_labels(baseline, labels_for({**base_options, **overrides}))
        except Exception as exc:
            row = {"error": f"{type(exc).__name__}: {exc}"}
        rows.append({"variant": name, **row})
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS + ["error"])
