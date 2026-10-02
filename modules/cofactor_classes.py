"""Named cofactor families and configurable network-distance rules.

The legacy interface accepts ``CLASS=ANGSTROMS`` overrides.  This module
keeps that form while adding named families and structured YAML/JSON rules.
No built-in family changes the default cutoff; an explicit rule is required
so existing analyses remain reproducible.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, FrozenSet, Iterable, Mapping, Optional, Sequence, Tuple, Union

from .structure_processing import METAL_ION_RESNAMES


@dataclass(frozen=True)
class CofactorClassRule:
    """A named cofactor family and its optional shell-distance override."""

    name: str
    residue_names: FrozenSet[str]
    cutoff_A: Optional[float] = None


# These are deliberately conservative starter families.  They identify the
# reference structures supported by the project without changing behavior for
# any structure unless a cutoff is explicitly configured.
_BUILTIN_RESIDUE_NAMES: Dict[str, FrozenSet[str]] = {
    "metal_ion": frozenset(METAL_ION_RESNAMES),
    "heme": frozenset({"HEM"}),
    "iron_sulfur_cluster": frozenset({"SF4", "FES", "FS4"}),
    "metallo_cluster": frozenset({"OEX", "OEC", "ICS", "CLF", "HCA"}),
    "organic_cofactor": frozenset(
        {
            "ATP", "ADP", "AMP", "COA", "COB", "FAD", "FMN", "NAD", "NAP",
            "PLP", "PQQ", "SAM", "THF",
        }
    ),
}

_CLASS_ALIASES = {
    "metal": "metal_ion",
    "metals": "metal_ion",
    "metal-ion": "metal_ion",
    "metal_ions": "metal_ion",
    "heme": "heme",
    "hemes": "heme",
    "iron-sulfur": "iron_sulfur_cluster",
    "iron_sulfur": "iron_sulfur_cluster",
    "iron-sulfur-cluster": "iron_sulfur_cluster",
    "iron_sulfur_clusters": "iron_sulfur_cluster",
    "cluster": "metallo_cluster",
    "clusters": "metallo_cluster",
    "metallo": "metallo_cluster",
    "organic": "organic_cofactor",
    "organics": "organic_cofactor",
}


def canonical_class_name(name: str) -> str:
    """Normalize a family name while retaining readable custom names."""
    normalized = str(name).strip().lower().replace(" ", "_")
    normalized = normalized.replace("-", "_")
    return _CLASS_ALIASES.get(normalized, normalized)


def builtin_cofactor_classes() -> Dict[str, CofactorClassRule]:
    """Return the built-in family registry without mutable shared state."""
    return {
        name: CofactorClassRule(name, residue_names)
        for name, residue_names in _BUILTIN_RESIDUE_NAMES.items()
    }


def _residue_names(value: Any) -> FrozenSet[str]:
    if value is None:
        return frozenset()
    if isinstance(value, str):
        values: Iterable[str] = value.split(",")
    elif isinstance(value, Sequence):
        values = (str(item) for item in value)
    else:
        raise ValueError("cofactor class residues must be a string or a list")
    return frozenset(str(item).strip().upper() for item in values if str(item).strip())


def _cutoff(value: Any, label: str) -> Optional[float]:
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a positive number in Å") from exc
    if result <= 0:
        raise ValueError(f"{label} must be greater than zero")
    return result


def _source_mapping(config: Mapping[str, Any]) -> Mapping[str, Any]:
    for key in ("cofactor_classes", "cofactor_class_rules", "cofactor_class_cutoffs"):
        nested = config.get(key)
        if isinstance(nested, Mapping):
            return nested
    return config


def normalize_cofactor_class_config(
    config: Optional[Mapping[str, Any]],
) -> Dict[str, CofactorClassRule]:
    """Normalize flat cutoff maps or structured family-rule mappings.

    Accepted forms include::

        {"metal": 2.8, "heme": 3.3}
        {"heme": {"residues": ["HEM"], "cutoff_A": 3.3}}
    """
    if not config:
        return {}
    source = _source_mapping(config)
    normalized: Dict[str, CofactorClassRule] = {}
    builtins = builtin_cofactor_classes()

    for raw_name, raw_rule in source.items():
        name = canonical_class_name(str(raw_name))
        if isinstance(raw_rule, CofactorClassRule):
            normalized[name] = CofactorClassRule(
                name=name,
                residue_names=raw_rule.residue_names,
                cutoff_A=raw_rule.cutoff_A,
            )
            continue
        base = builtins.get(name)
        if base is None:
            # A residue name is a useful shorthand for the existing CLI form,
            # e.g. ``HEM=3.3``. Unknown family names must use structured rules
            # with an explicit ``residues`` list.
            residue_name = str(raw_name).strip().upper()
            matching_builtin = next(
                (
                    (builtin_name, rule)
                    for builtin_name, rule in builtins.items()
                    if residue_name in rule.residue_names
                ),
                None,
            )
            if matching_builtin is not None:
                name, base = matching_builtin
            else:
                base = CofactorClassRule(name, frozenset())

        if isinstance(raw_rule, Mapping):
            residue_value = raw_rule.get(
                "residues",
                raw_rule.get("resnames", raw_rule.get("residue_names")),
            )
            residue_names = (
                _residue_names(residue_value)
                if residue_value is not None
                else base.residue_names
            )
            cutoff_value = raw_rule.get(
                "cutoff_A",
                raw_rule.get("cutoff", raw_rule.get("distance_cutoff")),
            )
        else:
            residue_names = base.residue_names
            cutoff_value = raw_rule

        normalized[name] = CofactorClassRule(
            name=name,
            residue_names=residue_names,
            cutoff_A=_cutoff(cutoff_value, f"cofactor class {raw_name!r} cutoff"),
        )
    return normalized


def merge_cofactor_class_configs(
    *configs: Optional[Mapping[str, Any]],
) -> Dict[str, CofactorClassRule]:
    """Merge rules from low to high precedence."""
    merged: Dict[str, CofactorClassRule] = {}
    for config in configs:
        for name, rule in normalize_cofactor_class_config(config).items():
            previous = merged.get(name)
            if previous is not None:
                residue_names = rule.residue_names or previous.residue_names
                cutoff_A = rule.cutoff_A if rule.cutoff_A is not None else previous.cutoff_A
                rule = CofactorClassRule(name, residue_names, cutoff_A)
            merged[name] = rule
    return merged


def resolve_cofactor_classes(
    cofactor_names: Sequence[str],
    config: Optional[Mapping[str, Any]] = None,
) -> Tuple[str, ...]:
    """Return all families represented by the selected cofactor names."""
    names = {str(value).strip().upper() for value in cofactor_names if str(value).strip()}
    rules = builtin_cofactor_classes()
    rules.update(normalize_cofactor_class_config(config))
    matched = {
        name for name, rule in rules.items() if names.intersection(rule.residue_names)
    }

    if not matched:
        return ("organic_cofactor",) if names else tuple()

    covered = set().union(*(rules[name].residue_names for name in matched))
    if names - covered:
        matched.add("organic_cofactor")
    return tuple(sorted(matched))


def resolve_effective_distance_cutoff(
    cofactor_names: Sequence[str],
    fallback: float,
    config: Optional[Mapping[str, Any]] = None,
) -> float:
    """Resolve a shell cutoff while preserving the legacy fallback.

    If several cofactor families are selected, the largest explicitly
    configured cutoff is used. Unconfigured families keep the fallback in the
    maximum so a mixed analysis does not silently lose network atoms.
    """
    fallback_value = _cutoff(fallback, "fallback distance cutoff")
    if fallback_value is None:
        raise ValueError("fallback distance cutoff is required")
    rules = normalize_cofactor_class_config(config)
    if not rules:
        return fallback_value

    classes = list(resolve_cofactor_classes(cofactor_names, rules))
    # ``organic`` was the legacy catch-all for every non-metal cofactor. Keep
    # that behavior when an explicit organic override is supplied, while the
    # named families above remain available for more targeted rules.
    names = {str(value).strip().upper() for value in cofactor_names if str(value).strip()}
    if (
        "organic_cofactor" in rules
        and any(name not in METAL_ION_RESNAMES for name in names)
        and "organic_cofactor" not in classes
    ):
        classes.append("organic_cofactor")
    configured = [rules[name].cutoff_A for name in classes if name in rules]
    configured = [value for value in configured if value is not None]
    if not configured:
        return fallback_value
    unconfigured = [
        name
        for name in classes
        if name not in rules or rules[name].cutoff_A is None
    ]
    # An explicit organic rule is a catch-all for non-metal families, just as
    # it was before named families were introduced.
    if "organic_cofactor" in rules and rules["organic_cofactor"].cutoff_A is not None:
        unconfigured = [name for name in unconfigured if name == "metal_ion"]
    if unconfigured:
        configured.append(fallback_value)
    return max(configured)


def load_cofactor_class_config(path: Union[str, Path]) -> Mapping[str, Any]:
    """Load a YAML or JSON class-rule file."""
    config_path = Path(path).expanduser()
    if not config_path.is_file():
        raise FileNotFoundError(f"Cofactor class configuration not found: {config_path}")
    if config_path.suffix.lower() in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore
        except ImportError as exc:
            raise RuntimeError("PyYAML is required to read YAML cofactor class configuration") from exc
        payload = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    else:
        payload = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise ValueError("Cofactor class configuration must contain a mapping")
    return payload
