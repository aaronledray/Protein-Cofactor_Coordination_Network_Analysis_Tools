"""Write the log of a conserved-interaction analysis: tidy tables, a manifest with provenance, and a readable summary (experimental).

The log is a folder. ``manifest.json`` says what the log is, what data and parameters produced it and, for every table, its file, row
count, columns and the rule behind each column (the legends that :mod:`modules.conserved_interactions` attaches to its results);
each table is a tidy CSV; ``structures.csv`` lists the structures with their group and sequence; ``SUMMARY.md`` is the same
information for a person. The log is a description of the structures given to the analysis. Directives are optional and
descriptive (below); without them the log contains none.

Directives are optional and descriptive. ``contact_directives`` and ``cofactor_directives`` restate rows of the vote tables, one directive per row
and (for cofactor residues) class, in three tiers, each with its per-group evidence and a fixed statement; ``write_log(..., directives=...)`` adds
them. A directive makes no claim about function, stability or the effect of changing a position (``DIRECTIVE_POLICY``).

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

import numpy as np
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

TIERS = ("conserved_in_all", "conserved_in_subset", "family_specific")
DIRECTIVE_POLICY = (
    "A directive restates one row of a table in this log together with the evidence behind it: the groups in which the contact or cofactor class is "
    "conserved, and for every group the sequences that vote, the positives and the negative sequences it could still take. It is descriptive "
    "(claim_level 'descriptive'): it states what the structures in this set show, and makes no claim about the effect of keeping or changing the position. "
    "Tiers: conserved_in_all (the headline: every group passes the support rule), conserved_in_subset (some groups pass, or one passes and the others "
    "are not mostly absent), family_specific (exactly one group passes and every other group has voting sequences and shows the contact in fewer than "
    "half of them). Directives are generated from the tables by fixed rules; none is written by hand."
)

MODES = ("vote", "single_protein")
SINGLE_TABLES = {"directives": "one observed directive per contact or cofactor class: a restated table row with its per-structure evidence and a tier"}
SINGLE_MEANING = (
    "This log describes the contacts seen in the structures of one protein: in how many of those structures each residue pair is in contact, and how "
    "close the residues come in each. It is not evidence of conservation: there is one sequence, and structures of one protein are not independent "
    "(they may share a crystal form, a laboratory and a refinement protocol). It does not say that a contact is important, that it stabilizes or "
    "enables anything, or that changing it would matter. A contact whose median closest distance lies near the cutoff is flagged, because its presence "
    "depends on where the cutoff falls. Hydrogens, where a rule uses them, are placed from the heavy atoms and are a model, not an observation. "
    "Waters and small ligands on metals are reported per structure."
)
OBSERVED_TIERS = ("observed_in_all_structures", "observed_in_some_structures", "observed_in_one_structure")
OBSERVED_PHRASES = {"observed_in_all_structures": "observed in all structures", "observed_in_some_structures": "observed in some structures",
                    "observed_in_one_structure": "observed in one structure"}
SINGLE_DIRECTIVE_POLICY = (
    "A directive restates one row of a table in this log together with its evidence: the structures in which the contact or cofactor class is seen, "
    "out of the structures that show both positions, and the closest distance in them. It is an observation (claim_level 'observed') about the structures "
    "of one protein: it is not evidence of conservation and makes no claim about the effect of keeping or changing the position. Tiers: "
    "observed_in_all_structures (in contact in every structure that shows both positions, and in at least two), observed_in_one_structure (in contact in "
    "exactly one structure, whether or not other structures show both positions) and observed_in_some_structures (all other cases). Directives are "
    "generated from the tables by fixed rules; none is written by hand."
)

TIER_PHRASES = {"conserved_in_all": "conserved in all groups", "conserved_in_subset": "conserved in some groups only",
                "family_specific": "conserved in one group and mostly absent from the others"}

KNOWN_TABLES = {
    "directives": "one directive per conserved (or, in single-protein mode, observed) contact or cofactor class: a restated table row with its evidence and a tier",
    "contact_persistence": "residue pairs in contact in at least one structure of one protein: in how many structures, and the closest distance in them (not a vote)",
    "cofactor_persistence": "residues in contact with a cofactor feature across the structures of one protein: structures per class and distances (not a vote)",
    "conserved_contacts": "residue pairs voted over sequences, with the chemistry that is conserved and the depth of the vote behind each pair",
    "state_contacts": "residue-pair votes calculated independently within each caller-supplied state (descriptive; not a test of a state effect)",
    "state_structure_votes": "eligible and positive structure counts behind each residue-pair vote, per state and sequence",
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


def _is_true(value) -> bool:
    """True only for a real boolean True (a NaN from a left merge is not true)."""
    return isinstance(value, (bool, np.bool_)) and bool(value)


def _tier(passing: Sequence[str], groups: Sequence[str], fractions: Mapping[str, float], voters: Mapping[str, int], specific_below: float) -> str:
    if len(passing) == len(groups):
        return "conserved_in_all"
    if len(passing) == 1 and all(voters[g] >= 1 and fractions[g] < specific_below for g in groups if g not in passing):
        return "family_specific"
    return "conserved_in_subset"


def _directives(vote: pd.DataFrame, structures: pd.DataFrame, *, keys: Sequence[str], flag: str, directive: str, source_table: str,
                subject, identifier, chemistry=None, support: float, min_sequences: int, specific_below: float) -> pd.DataFrame:
    """One directive per row of ``vote`` that passes the support rule in at least one group (columns ``<flag>_n_<g>`` and ``<flag>_<g>``)."""
    groups = sorted(structures["group"].unique())
    rows = []
    for record in vote.to_dict("records"):
        voters = {g: int(record.get(f"{flag}_n_{g}", 0) or 0) for g in groups}
        fractions = {g: float(record.get(f"{flag}_{g}", float("nan"))) for g in groups}
        positives = {g: int(round(fractions[g] * voters[g])) if voters[g] else 0 for g in groups}
        passing = [g for g in groups if voters[g] >= min_sequences and positives[g] / voters[g] >= support]
        if not passing:
            continue
        tier = _tier(passing, groups, {g: (positives[g] / voters[g] if voters[g] else float("nan")) for g in groups}, voters, specific_below)
        counts = "; ".join(f"{g} {positives[g]}/{voters[g]}" for g in groups)
        subject_text = subject(record)
        row = {"directive_id": f"{directive}:{identifier(record)}", "directive": directive, "source_table": source_table,
               **{k: record[k] for k in keys}, "subject": subject_text, "tier": tier, "conserved_in": "+".join(passing),
               "failing_groups": "+".join(g for g in groups if g not in passing)}
        for g in groups:
            allowed = int(math.floor(voters[g] * (1 - support) + 1e-9))
            row[f"voters_{g}"] = voters[g]
            row[f"positives_{g}"] = positives[g]
            row[f"spare_negatives_{g}"] = allowed - (voters[g] - positives[g])
        row["chemistry"] = chemistry(record) if chemistry else ""
        silent = [g for g in groups if voters[g] == 0]
        note = f"; no voting sequences in {' and '.join(silent)}" if silent and tier != "conserved_in_all" else ""
        row["statement"] = f"{subject_text}: {TIER_PHRASES[tier]} (conserved in {' and '.join(passing)}; sequences with it / sequences voting: {counts}{note})."
        row["claim_level"] = "descriptive"
        rows.append(row)
    return pd.DataFrame(rows)


def contact_directives(contacts: pd.DataFrame, structures: pd.DataFrame, *, support: float = 0.90, min_sequences: int = 5,
                       specific_below: float = 0.5) -> pd.DataFrame:
    """Descriptive ``conserved_contact`` directives from the output of ``conserved_contacts``: one per pair conserved in at least one group.

    The chemistry column lists the contact classes and flags that are conserved for the pair.
    """
    from .conserved_interactions import CONTACT_CLASSES, CONTACT_FLAGS, SULFUR_FLAGS

    names = [c for c in (*CONTACT_CLASSES, *CONTACT_FLAGS, *SULFUR_FLAGS) if f"conserved_{c}" in contacts.columns]
    return _directives(contacts, structures, keys=["position_a", "position_b"], flag="contact", directive="conserved_contact", source_table="conserved_contacts",
                       subject=lambda r: f"contact {r['position_a']}-{r['position_b']}", identifier=lambda r: f"{r['position_a']}-{r['position_b']}", chemistry=lambda r: ";".join(c for c in names if _is_true(r.get(f"conserved_{c}"))),
                       support=support, min_sequences=min_sequences, specific_below=specific_below)


def cofactor_directives(vote: pd.DataFrame, structures: pd.DataFrame, flags: Sequence[str], *, support: float = 0.90, min_sequences: int = 5,
                        specific_below: float = 0.5) -> pd.DataFrame:
    """Descriptive ``conserved_cofactor_contact`` directives from a ``sequence_vote`` over (feature, position): one per class in ``flags``."""
    frames = [_directives(vote, structures, keys=["feature", "position"], flag=flag, directive="conserved_cofactor_contact", source_table="cofactor_contacts",
                          subject=lambda r, flag=flag: f"{r['feature']} position {r['position']} class {flag}",
                          identifier=lambda r, flag=flag: f"{r['feature']}:{r['position']}:{flag}",
                          support=support, min_sequences=min_sequences, specific_below=specific_below).assign(cofactor_class=flag) for flag in flags]
    return pd.concat([f for f in frames if len(f)], ignore_index=True) if any(len(f) for f in frames) else pd.DataFrame()


def observed_tier(in_contact: int, present: int) -> str:
    """Tier of a contact seen in ``in_contact`` of the ``present`` structures that show both positions."""
    if in_contact == 1:
        return "observed_in_one_structure"
    return "observed_in_all_structures" if in_contact == present else "observed_in_some_structures"


def _observed_rows(table: pd.DataFrame, *, directive: str, source_table: str, keys: Sequence[str], subject, identifier, present: str, seen: str,
                   low: str, mid: str, high: str, total: int, near=None, chemistry=None) -> list:
    rows = []
    for r in table.to_dict("records"):
        k, n = int(r[seen]), int(r[present])
        if k < 1:
            continue
        tier = observed_tier(k, n)
        text = subject(r)
        note = ("; within the near-cutoff margin of the cutoff" if near and _is_true(r.get(near)) else "") + ("; one structure, so nothing is shown about persistence" if total == 1 else "")
        rows.append({"directive_id": f"{directive}:{identifier(r)}", "directive": directive, "source_table": source_table, **{c: r[c] for c in keys}, "subject": text,
                     "tier": tier, "structures_present": n, "structures_in_contact": k, "closest_min": r[low], "closest_median": r[mid], "closest_max": r[high],
                     "near_cutoff": bool(r[near]) if near else "", "chemistry": chemistry(r) if chemistry else "",
                     "statement": f"{text}: {OBSERVED_PHRASES[tier]} (in contact in {k} of {n} structures with both positions present; closest distance median {r[mid]:.2f} A, "
                                  f"range {r[low]:.2f} to {r[high]:.2f} A{note}).",
                     "claim_level": "observed"})
    return rows


def observed_contact_directives(persistence: pd.DataFrame, structures: pd.DataFrame) -> pd.DataFrame:
    """``observed_contact`` directives from :func:`modules.single_protein.contact_persistence`: one per pair, with its per-structure evidence."""
    rows = _observed_rows(persistence, directive="observed_contact", source_table="contact_persistence", keys=["position_a", "position_b"],
                          subject=lambda r: f"contact {r['position_a']}-{r['position_b']}", identifier=lambda r: f"{r['position_a']}-{r['position_b']}",
                          present="structures_present", seen="structures_in_contact", low="closest_min", mid="closest_median", high="closest_max",
                          total=int(structures["structure_id"].nunique()), near="near_cutoff", chemistry=lambda r: r.get("chemistry", ""))
    return pd.DataFrame(rows)


def observed_cofactor_directives(persistence: pd.DataFrame, structures: pd.DataFrame, flags: Sequence[str]) -> pd.DataFrame:
    """``observed_cofactor_contact`` directives from :func:`modules.single_protein.cofactor_persistence`: one per (feature, position, class) seen."""
    rows = []
    for flag in flags:
        rows += [dict(r, cofactor_class=flag) for r in _observed_rows(
            persistence, directive="observed_cofactor_contact", source_table="cofactor_persistence", keys=["feature", "position"],
            subject=lambda r, flag=flag: f"{r['feature']} position {r['position']} class {flag}", identifier=lambda r, flag=flag: f"{r['feature']}:{r['position']}:{flag}",
            present="structures_recorded", seen=f"{flag}_structures", low="distance_min", mid="distance_median", high="distance_max",
            total=int(structures["structure_id"].nunique()))]
    return pd.DataFrame(rows)


def write_log(directory, tables: Mapping[str, pd.DataFrame], *, structures: pd.DataFrame, title: str,
              parameters: Optional[Mapping[str, Any]] = None, notes: Sequence[str] = (), created: Optional[str] = None,
              overwrite: bool = False, directives: Optional[pd.DataFrame] = None, mode: str = "vote") -> Path:
    """Write the log folder and return its path.

    ``tables`` maps a name (letters, digits and underscores) to a DataFrame; each is written as ``<name>.csv`` and described in the
    manifest with its legend (``frame.attrs["rules"]``) when it has one. ``structures`` has ``structure_id, group, sequence_id`` and is
    written as ``structures.csv``. ``parameters`` (JSON-serializable) records how the tables were made; ``support`` and
    ``min_sequences`` in it set the group-depth statement (defaults 0.90 and 5). An existing log is not overwritten unless ``overwrite``;
    other files in the folder are never touched. ``directives`` (from :func:`contact_directives` and :func:`cofactor_directives`, concatenated) is written as
    ``directives.csv`` with ``DIRECTIVE_POLICY`` in the manifest and summary; without it the log contains no directive.
    ``mode="single_protein"`` (the structures of one protein; see :mod:`modules.single_protein`) writes the log with ``SINGLE_MEANING`` and accepts only
    directives from :func:`observed_contact_directives` and :func:`observed_cofactor_directives` (tiers ``OBSERVED_TIERS``, claim level ``observed``);
    in the default mode only the descriptive directives are accepted.
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, not {mode!r}")
    single = mode == "single_protein"
    if single and structures["sequence_id"].nunique() != 1:
        raise ValueError("a single-protein log needs exactly one sequence_id")
    directory = Path(directory)
    for name in tables:
        if not re.fullmatch(r"[A-Za-z0-9_]+", name) or name in ("manifest", "structures"):
            raise ValueError(f"table name {name!r} must be letters, digits and underscores and not 'manifest' or 'structures'")
    missing = {"structure_id", "group", "sequence_id"} - set(structures.columns)
    if missing:
        raise ValueError(f"structures lacks {sorted(missing)}")
    if (directory / "manifest.json").exists() and not overwrite:
        raise FileExistsError(f"{directory / 'manifest.json'} exists; pass overwrite=True to replace this log")
    if directives is not None:
        missing_columns = {"directive_id", "tier", "statement", "claim_level"} - set(directives.columns)
        tiers, level = (OBSERVED_TIERS, "observed") if single else (TIERS, "descriptive")
        if missing_columns or set(directives["tier"]) - set(tiers) or set(directives["claim_level"]) != {level}:
            raise ValueError(f"directives must have {sorted(missing_columns) or f'tiers from {tiers} and claim_level {level}'} in mode {mode}")
        tables = {**tables, "directives": directives}
    directory.mkdir(parents=True, exist_ok=True)
    parameters = dict(parameters or {})
    support, minimum = float(parameters.get("support", 0.90)), int(parameters.get("min_sequences", 5))
    depth = group_depth(structures, support=support, min_sequences=minimum)

    kept = ["structure_id", "group", "sequence_id"] + (["state"] if "state" in structures.columns else [])
    structures[kept].drop_duplicates("structure_id").to_csv(directory / "structures.csv", index=False)
    described = {}
    for name, frame in tables.items():
        frame.to_csv(directory / f"{name}.csv", index=False)
        legend = rules_of(frame)
        meaning = _descriptions(legend, list(frame.columns))
        described[name] = {"file": f"{name}.csv", "rows": int(len(frame)), "what": (SINGLE_TABLES if single else {}).get(name) or KNOWN_TABLES.get(name, ""),
                           "columns": [{"name": c, "type": str(frame[c].dtype), **({"rule": meaning[c]} if c in meaning else {})} for c in frame.columns],
                           "rules": _json_safe(legend.to_dict("records"))}
    manifest = {
        "log_format": LOG_FORMAT, "log_version": LOG_VERSION, "title": title, "created": created or datetime.date.today().isoformat(),
        "tool": {"repository_commit": _commit(), "module": "modules.design_log", "pandas": pd.__version__},
        "mode": mode, "meaning": SINGLE_MEANING if single else LOG_MEANING,
        "dataset": {"structures": int(structures["structure_id"].nunique()), "sequences": int(structures["sequence_id"].nunique()),
                    **({"states": _json_safe(structures.drop_duplicates("structure_id")["state"].astype(str).value_counts().to_dict())} if single and "state" in structures.columns
                       else {} if single else {"groups": _json_safe(depth.set_index("group").to_dict("index"))}), "structures_file": "structures.csv"},
        "parameters": _json_safe(parameters), "tables": described, "notes": list(notes),
    }
    if directives is not None:
        manifest["directive_policy"] = SINGLE_DIRECTIVE_POLICY if single else DIRECTIVE_POLICY
        manifest["directive_tiers"] = _json_safe(directives["tier"].value_counts().to_dict())
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=False) + "\n")
    (directory / "SUMMARY.md").write_text(_summary(manifest, tables, depth, minimum))
    return directory


def _summary(manifest: Dict[str, Any], tables: Mapping[str, pd.DataFrame], depth: pd.DataFrame, minimum: int) -> str:
    if manifest.get("mode") == "single_protein":
        return _summary_single(manifest, tables)
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
    state_contacts = tables.get("state_contacts")
    if state_contacts is not None and len(state_contacts):
        support = float(manifest["parameters"].get("support", 0.90))
        minimum = int(manifest["parameters"].get("min_sequences", 5))
        lines += ["## State-stratified contacts", "", "Each state uses the same sequence-vote rule independently. These counts are descriptive and do not test for a state effect.", "",
                  "| State | Group | candidate pairs | below minimum coverage | passing the group rule |", "|---|---|---:|---:|---:|"]
        for state in sorted(state_contacts["state"].astype(str).unique()):
            subset = state_contacts[state_contacts["state"].astype(str) == state]
            for group in sorted(manifest["dataset"]["groups"]):
                n = subset[f"contact_n_{group}"]
                passing = subset[f"passes_contact_{group}"]
                lines.append(f"| {state} | {group} | {len(subset)} | {int((n < minimum).sum())} | {int(passing.sum())} |")
        lines.append("")
    ligands = tables.get("metal_ligands")
    if ligands is not None and len(ligands):
        lines += ["## Metal ligands", "", "Non-protein atoms bound to metals, per structure (not voted). Structures with one, by metal and species:", ""]
        counts = ligands.groupby(["metal_element", "ligand_resname"])["structure_id"].nunique().sort_values(ascending=False)
        lines += ["| metal | ligand | structures |", "|---|---|---|"] + [f"| {m} | {l} | {n} |" for (m, l), n in counts.items()] + [""]
    directives = tables.get("directives")
    if directives is not None:
        lines += ["## Directives", "", manifest["directive_policy"], "", "| tier | directives |", "|---|---|"]
        lines += [f"| {tier} | {int((directives['tier'] == tier).sum())} |" for tier in TIERS]
        lines += ["", "| kind | directives |", "|---|---|"] + [f"| {k} | {v} |" for k, v in directives["directive"].value_counts().items()] + [""]
    if manifest["notes"]:
        lines += ["## Notes", ""] + [f"- {n}" for n in manifest["notes"]] + [""]
    lines += ["## Parameters", "", "```json", json.dumps(manifest["parameters"], indent=1), "```", ""]
    for name, frame in tables.items():
        legend = rules_of(frame)
        if len(legend):
            lines += [f"## Rules in force: {name}", "", format_rules(legend)]
    return "\n".join(lines) + "\n"


def _summary_single(manifest: Dict[str, Any], tables: Mapping[str, pd.DataFrame]) -> str:
    dataset = manifest["dataset"]
    lines = [f"# {manifest['title']}", "", f"Log format {manifest['log_format']} version {manifest['log_version']}, mode single_protein, created {manifest['created']}. "
             "Every table, column and rule is described in `manifest.json`.", "", "## What this log means", "", manifest["meaning"], "", "## Dataset", "",
             f"{dataset['structures']} structure(s) of one protein (one sequence id)." + (" States: " + ", ".join(f"{k} {v}" for k, v in dataset["states"].items()) + "." if dataset.get("states") else "")]
    if dataset["structures"] == 1:
        lines += ["", "There is one structure, so nothing can be said about whether a contact persists across structures."]
    lines += ["", "## Tables", "", "| Table | file | rows | what it holds |", "|---|---|---|---|"]
    for name, info in manifest["tables"].items():
        lines.append(f"| {name} | `{info['file']}` | {info['rows']} | {info['what']} |")
    lines.append("")
    persistence = tables.get("contact_persistence")
    if persistence is not None and len(persistence):
        every = int((persistence["structures_in_contact"] == persistence["structures_present"]).sum())
        lines += ["## Contacts", "", f"{len(persistence)} residue pairs are in contact in at least one structure; {every} in every structure that shows both positions; "
                  f"{int(persistence['near_cutoff'].sum())} have a median closest distance near the cutoff.", ""]
    ligands = tables.get("metal_ligands")
    if ligands is not None and len(ligands):
        counts = ligands.groupby(["metal_element", "ligand_resname"])["structure_id"].nunique().sort_values(ascending=False)
        lines += ["## Metal ligands", "", "Non-protein atoms bound to metals, per structure. Structures with one, by metal and species:", "",
                  "| metal | ligand | structures |", "|---|---|---|"] + [f"| {m} | {l} | {n} |" for (m, l), n in counts.items()] + [""]
    directives = tables.get("directives")
    if directives is not None:
        lines += ["## Directives", "", manifest["directive_policy"], "", "| tier | directives |", "|---|---|"]
        lines += [f"| {tier} | {int((directives['tier'] == tier).sum())} |" for tier in OBSERVED_TIERS]
        lines += ["", "| kind | directives |", "|---|---|"] + [f"| {k} | {v} |" for k, v in directives["directive"].value_counts().items()] + [""]
    if manifest["notes"]:
        lines += ["## Notes", ""] + [f"- {n}" for n in manifest["notes"]] + [""]
    lines += ["## Parameters", "", "```json", json.dumps(manifest["parameters"], indent=1), "```", ""]
    for name, frame in tables.items():
        legend = rules_of(frame)
        if len(legend):
            lines += [f"## Rules in force: {name}", "", format_rules(legend)]
    return "\n".join(lines) + "\n"
