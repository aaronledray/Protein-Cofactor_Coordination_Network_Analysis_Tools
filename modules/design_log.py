"""Write the log of a conserved-interaction analysis: tidy tables, a manifest with provenance, and a readable summary (experimental).

The log is a folder. ``manifest.json`` says what the log is, what data and parameters produced it and, for every table, its file, row
count, columns and the rule behind each column (the legends that :mod:`modules.conserved_interactions` attaches to its results);
each table is a tidy CSV; ``structures.csv`` lists the structures with their group and sequence; ``SUMMARY.md`` is the same
information for a person. The log is a description of the structures given to the analysis. It carries no directives: what a
directive may say is a separate decision, and nothing here writes one.

The writer takes tables that the analysis functions already return (``conserved_contacts(..., coverage=True)``, the cofactor vote,
``load_metal_ligands``, the audits, ...); composing them into a pipeline is the caller's.
"""

import datetime
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

import pandas as pd

from .conserved_interactions import format_rules, rules_of

LOG_FORMAT = "cni-log"
LOG_VERSION = 1

LOG_MEANING = (
    "This log describes which contacts, with what chemistry, recur in the sequences and groups of the structure set it was made from. "
    "It does not say that a contact is important, that it stabilizes or enables anything, or that changing it would matter; conservation can "
    "reflect shared ancestry as well as constraint. A contact conserved in a group is conserved in that group's sequences as counted in the "
    "dataset table, and a group with few sequences tolerates few or no exceptions. Hydrogens, where a rule uses them, are placed from the "
    "heavy atoms and are a model, not an observation. Waters and small ligands on metals are reported per structure and are never voted on."
)

KNOWN_TABLES = {
    "conserved_contacts": "residue pairs voted over sequences, with the chemistry that is conserved and the depth of the vote behind each pair",
    "cofactor_contacts": "residues in contact with a cofactor feature, voted over sequences (coordination tiers where computed)",
    "metal_ligands": "non-protein atoms bound to metals, per structure (descriptive, not voted)",
    "water_bridges": "residues bridged to a cofactor oxygen by one modeled water, per structure",
    "audit_anchor": "per structure: is the anchor residue where the alignment says it should be",
    "audit_backbone_only": "per structure: is the model backbone-only",
}


def _json_safe(value: Any) -> Any:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if isinstance(value, Mapping):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(v) for v in value]
    if hasattr(value, "item"):
        return _json_safe(value.item())
    return value


def _commit() -> Optional[str]:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=Path(__file__).resolve().parent, timeout=10)
        return out.stdout.strip() if out.returncode == 0 and out.stdout.strip() else None
    except Exception:
        return None


def _descriptions(rules: pd.DataFrame, columns: Sequence[str]) -> Dict[str, str]:
    """Match each column to the legend row with its name (a ``<placeholder>`` in a legend name matches any suffix)."""
    exact = {r.name: r.rule for r in rules.itertuples()}
    prefixes = [(r.name.split("<")[0], r.rule) for r in rules.itertuples() if "<" in r.name and r.name.split("<")[0]]
    out = {}
    for column in columns:
        if column in exact:
            out[column] = exact[column]
        else:
            for prefix, rule in prefixes:
                if column.startswith(prefix):
                    out[column] = rule
                    break
    return out


def group_depth(structures: pd.DataFrame, *, support: float = 0.90, min_sequences: int = 5) -> pd.DataFrame:
    """Per group: sequences, structures, the negative sequences a contact can take and still pass, and whether the group meets the minimum."""
    rows = []
    for group, frame in structures.groupby("group", sort=True):
        n = int(frame["sequence_id"].nunique())
        rows.append({"group": group, "sequences": n, "structures": int(frame["structure_id"].nunique()),
                     "tolerated_negatives": int(math.floor(n * (1 - support) + 1e-9)), "meets_minimum": n >= min_sequences})
    return pd.DataFrame(rows, columns=["group", "sequences", "structures", "tolerated_negatives", "meets_minimum"])


def write_log(directory, tables: Mapping[str, pd.DataFrame], *, structures: pd.DataFrame, title: str,
              parameters: Optional[Mapping[str, Any]] = None, notes: Sequence[str] = (), created: Optional[str] = None,
              overwrite: bool = False) -> Path:
    """Write the log folder and return its path.

    ``tables`` maps a name (letters, digits and underscores) to a DataFrame; each is written as ``<name>.csv`` and described in the
    manifest with its legend (``frame.attrs["rules"]``) when it has one. ``structures`` has ``structure_id, group, sequence_id`` and is
    written as ``structures.csv``. ``parameters`` (JSON-serializable) records how the tables were made; ``support`` and
    ``min_sequences`` in it set the group-depth statement (defaults 0.90 and 5). An existing log is not overwritten unless ``overwrite``;
    other files in the folder are never touched.
    """
    directory = Path(directory)
    for name in tables:
        if not re.fullmatch(r"[A-Za-z0-9_]+", name) or name in ("manifest", "structures"):
            raise ValueError(f"table name {name!r} must be letters, digits and underscores and not 'manifest' or 'structures'")
    missing = {"structure_id", "group", "sequence_id"} - set(structures.columns)
    if missing:
        raise ValueError(f"structures lacks {sorted(missing)}")
    if (directory / "manifest.json").exists() and not overwrite:
        raise FileExistsError(f"{directory / 'manifest.json'} exists; pass overwrite=True to replace this log")
    directory.mkdir(parents=True, exist_ok=True)
    parameters = dict(parameters or {})
    support, minimum = float(parameters.get("support", 0.90)), int(parameters.get("min_sequences", 5))
    depth = group_depth(structures, support=support, min_sequences=minimum)

    structures[["structure_id", "group", "sequence_id"]].drop_duplicates("structure_id").to_csv(directory / "structures.csv", index=False)
    described = {}
    for name, frame in tables.items():
        frame.to_csv(directory / f"{name}.csv", index=False)
        legend = rules_of(frame)
        meaning = _descriptions(legend, list(frame.columns))
        described[name] = {"file": f"{name}.csv", "rows": int(len(frame)), "what": KNOWN_TABLES.get(name, ""),
                           "columns": [{"name": c, "type": str(frame[c].dtype), **({"rule": meaning[c]} if c in meaning else {})} for c in frame.columns],
                           "rules": _json_safe(legend.to_dict("records"))}
    manifest = {
        "log_format": LOG_FORMAT, "log_version": LOG_VERSION, "title": title, "created": created or datetime.date.today().isoformat(),
        "tool": {"repository_commit": _commit(), "module": "modules.design_log", "pandas": pd.__version__},
        "meaning": LOG_MEANING,
        "dataset": {"structures": int(structures["structure_id"].nunique()), "sequences": int(structures["sequence_id"].nunique()),
                    "groups": _json_safe(depth.set_index("group").to_dict("index")), "structures_file": "structures.csv"},
        "parameters": _json_safe(parameters), "tables": described, "notes": list(notes),
    }
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=False) + "\n")
    (directory / "SUMMARY.md").write_text(_summary(manifest, tables, depth, minimum))
    return directory


def _summary(manifest: Dict[str, Any], tables: Mapping[str, pd.DataFrame], depth: pd.DataFrame, minimum: int) -> str:
    lines = [f"# {manifest['title']}", "", f"Log format {manifest['log_format']} version {manifest['log_version']}, created {manifest['created']}. "
             "Every table, column and rule is described in `manifest.json`.", "", "## What this log means", "", manifest["meaning"], "", "## Dataset", "",
             f"{manifest['dataset']['structures']} structures, {manifest['dataset']['sequences']} sequences.", "",
             "| Group | sequences | structures | negative sequences a contact can take and still pass | meets the minimum of voting sequences |", "|---|---|---|---|---|"]
    for r in depth.itertuples():
        lines.append(f"| {r.group} | {r.sequences} | {r.structures} | {r.tolerated_negatives} | {'yes' if r.meets_minimum else 'no'} |")
    for r in depth.itertuples():
        if r.tolerated_negatives == 0:
            lines += ["", f"Group {r.group} has {r.sequences} sequence(s) and tolerates no negative sequence: a contact conserved in this group is present in every one of them."]
        if not r.meets_minimum:
            lines += ["", f"Group {r.group} has fewer than the minimum of {minimum} voting sequences, so no contact can be conserved in it."]
    lines += ["", "## Tables", "", "| Table | file | rows | what it holds |", "|---|---|---|---|"]
    for name, info in manifest["tables"].items():
        lines.append(f"| {name} | `{info['file']}` | {info['rows']} | {info['what']} |")
    lines.append("")
    contacts = tables.get("conserved_contacts")
    if contacts is not None and "conserved_contact" in contacts.columns:
        lines += ["## Conserved contacts", "", f"{int(contacts['conserved_contact'].sum())} of {len(contacts)} candidate pairs are conserved in every group (the headline)."]
        if "conserved_in" in contacts.columns:
            counts = contacts["conserved_in"].replace("", "(none)").value_counts()
            lines += ["", "| conserved in | pairs |", "|---|---|"] + [f"| {k} | {v} |" for k, v in counts.items()]
        lines.append("")
    ligands = tables.get("metal_ligands")
    if ligands is not None and len(ligands):
        lines += ["## Metal ligands", "", "Non-protein atoms bound to metals, per structure (not voted). Structures with one, by metal and species:", ""]
        counts = ligands.groupby(["metal_element", "ligand_resname"])["structure_id"].nunique().sort_values(ascending=False)
        lines += ["| metal | ligand | structures |", "|---|---|---|"] + [f"| {m} | {l} | {n} |" for (m, l), n in counts.items()] + [""]
    if manifest["notes"]:
        lines += ["## Notes", ""] + [f"- {n}" for n in manifest["notes"]] + [""]
    lines += ["## Parameters", "", "```json", json.dumps(manifest["parameters"], indent=1), "```", ""]
    for name, frame in tables.items():
        legend = rules_of(frame)
        if len(legend):
            lines += [f"## Rules in force: {name}", "", format_rules(legend)]
    return "\n".join(lines) + "\n"
