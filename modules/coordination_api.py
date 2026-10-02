"""Importable, headless coordination-network analysis APIs.

The legacy SSCNA script remains responsible for its existing files and plots.
This module provides a side-effect-free single-structure API and a batch API
for downstream workflows such as model-label generation.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import pandas as pd
import numpy as np

from .io_utils import unpack_pdb_file
from .structure_processing import (
    COORD_LINK_COLUMNS,
    _coord_link_rows,
    _shell_label,
    identify_coordination_shells,
    identify_coordination_network,
    METAL_ION_RESNAMES,
)
from .motif_registry import motif_for_atom
from .cofactor_classes import resolve_effective_distance_cutoff


PathLike = Union[str, Path]

ATOM_COLUMNS = [
    "structure_id", "site_id", "shell", "residue_name", "residue_number",
    "chain", "insertion_code", "hetero_flag", "model_id", "atom_name", "element",
    "motif", "coordination_role", "x", "y", "z",
]

RESIDUE_COLUMNS = [
    "structure_id", "site_id", "cofactor_residue_name", "cofactor_residue_number",
    "cofactor_chain", "cofactor_insertion_code", "shell", "residue_name",
    "residue_number", "chain", "insertion_code", "hetero_flag",
    "atoms_involved", "minimum_distance_A",
]

# Opt-in residue schema (``include_model_id=True``): one row per residue per
# model, with the model column placed after the residue identity fields.
RESIDUE_COLUMNS_WITH_MODEL = (
    RESIDUE_COLUMNS[:12] + ["model_id"] + RESIDUE_COLUMNS[12:]
)

LINK_IDENTITY_COLUMNS = [
    "src_insertion_code", "src_hetero_flag",
    "dst_insertion_code", "dst_hetero_flag",
    "src_element", "dst_element",
    "direct_coordination",
]
LINK_COLUMNS = ["structure_id"] + COORD_LINK_COLUMNS + LINK_IDENTITY_COLUMNS

CONTACT_COLUMNS = [
    "structure_id", "site_id", "src_shell", "dst_shell", "link_type",
    "contact_role", "direct_coordination",
    "src_resname", "src_resnum", "src_chain", "src_insertion_code",
    "src_hetero_flag", "src_model_id", "src_atom", "src_element", "src_motif",
    "dst_resname", "dst_resnum", "dst_chain", "dst_insertion_code",
    "dst_hetero_flag", "dst_model_id", "dst_atom", "dst_element", "dst_motif",
    "distance_A",
]

_INFERRED_HBOND_CUTOFF = 3.2
_PRIMARY_MOTIF_TARGET_ELEMENTS = frozenset({"N", "O", "S"})


def _motif_label(atom: Dict[str, Any]) -> str:
    """Return the shared residue-specific motif for an atom."""
    residue_name = str(atom.get("residue", "")).upper()
    atom_name = str(atom.get("name", "")).upper()
    return motif_for_atom(residue_name, atom_name)


def _annotate_link_motifs(link_rows: Sequence[Dict[str, Any]]) -> None:
    """Populate motif labels in API link rows without changing legacy files."""
    for row in link_rows:
        row["src_moiety"] = _motif_label(
            {"residue": row.get("src_resname", ""), "name": row.get("src_atom", "")}
        )
        row["dst_moiety"] = _motif_label(
            {"residue": row.get("dst_resname", ""), "name": row.get("dst_atom", "")}
        )


def _is_metal_atom(atom: Dict[str, Any]) -> bool:
    element = str(atom.get("element", "")).upper()
    residue_name = str(atom.get("residue", "")).upper()
    return element in METAL_ION_RESNAMES or residue_name in METAL_ION_RESNAMES


def _is_primary_motif_contact(
    source: Dict[str, Any],
    destination: Dict[str, Any],
    source_shell_number: int,
    distance: float,
) -> bool:
    """Identify primary cofactor-motif contacts without calling them direct.

    Heme propionate oxygens are cofactor motifs that can form primary
    hydrogen-bonding/polar contacts with protein or water heteroatoms.  These
    contacts are rooted in the visualized network, while ``direct_coordination``
    remains reserved for validated metal-originating links such as HEM FE to a
    histidine nitrogen.
    """
    if source_shell_number != 0 or distance > _INFERRED_HBOND_CUTOFF:
        return False
    if _motif_label(source) != "heme_propionate":
        return False
    destination_element = str(destination.get("element", "")).upper()
    return destination_element in _PRIMARY_MOTIF_TARGET_ELEMENTS


def _motif_contact_rows(
    structure_id: str,
    site_id: str,
    cofactor_atoms: Sequence[Dict[str, Any]],
    shell_atoms: Dict[int, Sequence[Dict[str, Any]]],
    *,
    distance_cutoff: float,
    direct_coordination_cutoff: float,
) -> pd.DataFrame:
    """Return exact atom-pair contacts between adjacent coordination shells.

    This is intentionally additive to the legacy nearest-atom links. It keeps
    every qualifying pair from the selected shell atom sets, labels both atoms
    with their chemical motifs, and prevents contacts from being formed across
    different structure models.
    """
    transitions: List[Tuple[int, str, Sequence[Dict[str, Any]], int, str, Sequence[Dict[str, Any]]]] = []
    for shell_number in sorted(shell_atoms):
        source_shell_number = shell_number - 1
        source_label = "Cofactor" if source_shell_number == 0 else _shell_label(source_shell_number)
        source_atoms = cofactor_atoms if source_shell_number == 0 else shell_atoms.get(source_shell_number, [])
        transitions.append(
            (
                source_shell_number,
                source_label,
                source_atoms,
                shell_number,
                _shell_label(shell_number),
                shell_atoms[shell_number],
            )
        )

    rows: List[Dict[str, Any]] = []
    for src_number, src_label, sources, dst_number, dst_label, destinations in transitions:
        link_type = f"{src_label.lower()}->{dst_label.lower()}"
        for source in sources:
            source_model = source.get("model_id")
            source_coordinates = np.asarray(source["coordinates"], dtype=float)
            for destination in destinations:
                destination_model = destination.get("model_id")
                if source_model is not None and destination_model is not None and source_model != destination_model:
                    continue
                distance = float(np.linalg.norm(source_coordinates - np.asarray(destination["coordinates"], dtype=float)))
                if distance > distance_cutoff:
                    continue
                direct = bool(
                    src_number == 0
                    and _is_metal_atom(source)
                    and distance <= direct_coordination_cutoff
                )
                primary_motif_contact = _is_primary_motif_contact(
                    source,
                    destination,
                    src_number,
                    distance,
                )
                rows.append(
                    {
                        "structure_id": structure_id,
                        "site_id": site_id,
                        "src_shell": src_label,
                        "dst_shell": dst_label,
                        "link_type": link_type,
                        "contact_role": (
                            "direct_coordination"
                            if direct
                            else "primary_motif_contact" if primary_motif_contact else "shell_contact"
                        ),
                        "direct_coordination": direct,
                        "src_resname": source.get("residue", ""),
                        "src_resnum": source.get("residue_number", ""),
                        "src_chain": source.get("chain", ""),
                        "src_insertion_code": source.get("insertion_code", ""),
                        "src_hetero_flag": source.get("hetero_flag", ""),
                        "src_model_id": source_model,
                        "src_atom": source.get("name", ""),
                        "src_element": source.get("element", ""),
                        "src_motif": _motif_label(source),
                        "dst_resname": destination.get("residue", ""),
                        "dst_resnum": destination.get("residue_number", ""),
                        "dst_chain": destination.get("chain", ""),
                        "dst_insertion_code": destination.get("insertion_code", ""),
                        "dst_hetero_flag": destination.get("hetero_flag", ""),
                        "dst_model_id": destination_model,
                        "dst_atom": destination.get("name", ""),
                        "dst_element": destination.get("element", ""),
                        "dst_motif": _motif_label(destination),
                        "distance_A": distance,
                    }
                )

    rows.sort(
        key=lambda row: (
            row["structure_id"], row["site_id"], row["src_shell"], row["dst_shell"],
            str(row["src_resname"]), str(row["src_resnum"]), str(row["src_chain"]),
            str(row["src_atom"]), str(row["dst_resname"]), str(row["dst_resnum"]),
            str(row["dst_chain"]), str(row["dst_atom"]), row["distance_A"],
        )
    )
    return pd.DataFrame(rows, columns=CONTACT_COLUMNS)


def _as_names(names: Union[str, Sequence[str]]) -> List[str]:
    if isinstance(names, str):
        return [part.strip() for part in names.split(",") if part.strip()]
    return [str(name).strip() for name in names if str(name).strip()]


def _cofactor_site_groups(
    structure,
    cofactor_names: Sequence[str],
    cofactor_names2: Sequence[str],
    *,
    combinatorial: bool,
    cutoff: float,
    first_model_only: bool,
    model_mode: str = "pooled",
) -> List[set]:
    """Return cofactor residue sites, optionally merged into nearby clusters.

    ``pooled`` preserves the historical site identity, which is residue name,
    residue number, chain, insertion code, and hetero flag across all models.
    ``per-model`` adds the model ID to the site identity. Combinatorial
    merging only compares cofactors within the same model in either mode, so
    a cross-model proximity cannot create a false linked site.
    """
    if model_mode not in {"pooled", "per-model"}:
        raise ValueError("model_mode must be 'pooled' or 'per-model'")
    models = list(structure)
    if first_model_only:
        models = models[:1]
    names = {name.upper() for name in (*cofactor_names, *cofactor_names2)}
    atoms_by_site: Dict[Tuple[Any, ...], Dict[Any, List[np.ndarray]]] = {}
    for model in models:
        for chain in model:
            for residue in chain:
                if residue.get_resname().upper() not in names:
                    continue
                residue_id = residue.get_id()
                identity = (
                    residue.get_resname(),
                    residue_id[1],
                    chain.id,
                    str(residue_id[2] or "").strip(),
                    str(residue_id[0] or "").strip(),
                )
                key = (model.id, *identity) if model_mode == "per-model" else identity
                atoms_by_site.setdefault(key, {}).setdefault(model.id, []).extend(
                    np.asarray(atom.coord, dtype=float) for atom in residue
                )

    keys = sorted(atoms_by_site, key=lambda value: tuple(str(part) for part in value))
    if not combinatorial:
        return [{key} for key in keys]

    parent = list(range(len(keys)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for left in range(len(keys)):
        for right in range(left + 1, len(keys)):
            left_by_model = atoms_by_site[keys[left]]
            right_by_model = atoms_by_site[keys[right]]
            for model_id in set(left_by_model).intersection(right_by_model):
                distances = np.linalg.norm(
                    np.asarray(left_by_model[model_id])[:, None, :]
                    - np.asarray(right_by_model[model_id])[None, :, :],
                    axis=2,
                )
                if float(distances.min()) <= cutoff:
                    union(left, right)
                    break

    grouped: Dict[int, set] = {}
    for index, key in enumerate(keys):
        grouped.setdefault(find(index), set()).add(key)
    return [grouped[root] for root in sorted(grouped)]


def _effective_distance_cutoff(
    cofactor_names: Sequence[str],
    fallback: float,
    class_cutoffs: Optional[Dict[str, Any]],
) -> float:
    return resolve_effective_distance_cutoff(cofactor_names, fallback, class_cutoffs)


def _site_label(index: int, site_keys: set, model_mode: str) -> str:
    """Return a stable, human-readable site identifier."""
    if model_mode != "per-model":
        return f"site_{index}"
    model_ids = sorted(
        {key[0] for key in site_keys if len(key) == 6},
        key=str,
    )
    if len(model_ids) == 1:
        return f"site_{index}_model_{model_ids[0]}"
    return f"site_{index}_models_{'_'.join(str(value) for value in model_ids)}"


def _annotate_direct_coordination(
    rows: List[Dict[str, Any]],
    *,
    enabled: bool,
    cutoff: float,
) -> None:
    for row in rows:
        row["direct_coordination"] = bool(
            enabled
            and row.get("link_type") == "cofactor->pcs"
            and str(row.get("src_element", "")).upper() in METAL_ION_RESNAMES
            and float(row.get("distance_A", float("inf"))) <= cutoff
        )


def _structure_id(path: Path) -> str:
    name = path.name
    for suffix in (".gz", ".mmcif", ".cif", ".pdb"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return name


def _atom_identity(atom: Dict[str, Any]) -> Tuple[str, Any, str, str, str]:
    return (
        str(atom.get("residue", "")),
        atom.get("residue_number", ""),
        str(atom.get("chain", "")),
        str(atom.get("insertion_code", "") or "").strip(),
        str(atom.get("hetero_flag", "") or "").strip(),
    )


def _atom_rows(
    structure_id: str,
    site_id: str,
    shell_atoms: Sequence[Tuple[str, Sequence[Dict[str, Any]]]],
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for shell, atoms in shell_atoms:
        for atom in atoms:
            coordinates = list(atom["coordinates"])
            rows.append(
                {
                    "structure_id": structure_id,
                    "site_id": site_id,
                    "shell": shell,
                    "residue_name": atom.get("residue", ""),
                    "residue_number": atom.get("residue_number", ""),
                    "chain": atom.get("chain", ""),
                    "insertion_code": atom.get("insertion_code", ""),
                    "hetero_flag": atom.get("hetero_flag", ""),
                    "model_id": atom.get("model_id", ""),
                    "atom_name": atom.get("name", ""),
                    "element": atom.get("element", ""),
                    "motif": _motif_label(atom),
                    "coordination_role": "active_site_component",
                    "x": float(coordinates[0]),
                    "y": float(coordinates[1]),
                    "z": float(coordinates[2]),
                }
            )
    return pd.DataFrame(rows, columns=ATOM_COLUMNS)


def _normalized_atom_key(
    *,
    model_id: Any,
    residue_name: Any,
    residue_number: Any,
    chain: Any,
    insertion_code: Any,
    hetero_flag: Any,
    atom_name: Any,
) -> Tuple[str, str, str, str, str, str, str]:
    def normalize(value: Any) -> str:
        return "" if value is None else str(value).strip()

    return (
        normalize(model_id), normalize(residue_name), normalize(residue_number),
        normalize(chain), normalize(insertion_code), normalize(hetero_flag),
        normalize(atom_name),
    )


def _contact_atom_key(row: Dict[str, Any], prefix: str) -> Tuple[str, str, str, str, str, str, str]:
    return _normalized_atom_key(
        model_id=row.get(f"{prefix}_model_id"),
        residue_name=row.get(f"{prefix}_resname"),
        residue_number=row.get(f"{prefix}_resnum"),
        chain=row.get(f"{prefix}_chain"),
        insertion_code=row.get(f"{prefix}_insertion_code"),
        hetero_flag=row.get(f"{prefix}_hetero_flag"),
        atom_name=row.get(f"{prefix}_atom"),
    )


def _annotate_atom_roles(atoms: pd.DataFrame, contacts: pd.DataFrame) -> pd.DataFrame:
    """Separate geometric shell membership from validated network roles."""
    primary_keys = set()
    network_keys = set()
    for row in contacts.to_dict("records"):
        contact_role = str(row.get("contact_role", "")).strip().lower()
        direct_value = row.get("direct_coordination", False)
        is_direct = direct_value is True or str(direct_value).strip().lower() in {"true", "1", "yes"}
        is_primary = is_direct or contact_role == "primary_motif_contact"
        link_type = str(row.get("link_type", "")).strip().lower()
        is_network = link_type in {"pcs->scs", "scs->tcs"}
        for prefix in ("src", "dst"):
            key = _contact_atom_key(row, prefix)
            if is_primary:
                primary_keys.add(key)
            elif is_network:
                network_keys.add(key)

    annotated = atoms.copy()
    roles = []
    for row in annotated.to_dict("records"):
        key = _normalized_atom_key(
            model_id=row.get("model_id"),
            residue_name=row.get("residue_name"),
            residue_number=row.get("residue_number"),
            chain=row.get("chain"),
            insertion_code=row.get("insertion_code"),
            hetero_flag=row.get("hetero_flag"),
            atom_name=row.get("atom_name"),
        )
        if key in primary_keys:
            roles.append("primary_coordinator")
        elif key in network_keys:
            roles.append("network_context")
        else:
            roles.append("active_site_component")
    annotated["coordination_role"] = roles
    return annotated


def _residue_rows(
    structure_id: str,
    site_id: str,
    cofactor_atoms: Sequence[Dict[str, Any]],
    shell_atoms: Sequence[Tuple[str, Sequence[Dict[str, Any]]]],
    link_rows: Sequence[Dict[str, Any]],
    include_model_id: bool = False,
    contacts: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    cofactor_ids = list(dict.fromkeys(_atom_identity(atom) for atom in cofactor_atoms))
    cofactor_names = ";".join(str(value[0]) for value in cofactor_ids)
    cofactor_numbers = ";".join(str(value[1]) for value in cofactor_ids)
    cofactor_chains = ";".join(str(value[2]) for value in cofactor_ids)
    cofactor_insertions = ";".join(str(value[3]) for value in cofactor_ids)

    distance_by_residue: Dict[Tuple[str, Any, str, str, str], float] = {}
    for link in link_rows:
        key = (
            str(link.get("dst_resname", "")),
            link.get("dst_resnum", ""),
            str(link.get("dst_chain", "")),
            str(link.get("dst_insertion_code", "") or "").strip(),
            str(link.get("dst_hetero_flag", "") or "").strip(),
        )
        distance = float(link["distance_A"])
        previous = distance_by_residue.get(key)
        if previous is None or distance < previous:
            distance_by_residue[key] = distance

    model_distances: Dict[Tuple, float] = {}
    if include_model_id and contacts is not None and not contacts.empty:
        for row in contacts.itertuples():
            key = (
                row.dst_shell,
                str(row.dst_resname),
                row.dst_resnum,
                str(row.dst_chain),
                str(row.dst_insertion_code or "").strip(),
                str(row.dst_hetero_flag or "").strip(),
                row.dst_model_id,
            )
            distance = round(float(row.distance_A), 3)  # links round to 3 dp
            if key not in model_distances or distance < model_distances[key]:
                model_distances[key] = distance

    rows: List[Dict[str, Any]] = []
    for shell, atoms in shell_atoms:
        grouped: Dict[Tuple, List[Dict[str, Any]]] = {}
        for atom in atoms:
            key = _atom_identity(atom)
            if include_model_id:
                key = key + (atom.get("model_id", ""),)
            grouped.setdefault(key, []).append(atom)
        for key, group in grouped.items():
            residue_name, residue_number, chain = key[:3]
            first = group[0]
            model_distance = (
                model_distances.get((shell,) + key)
                if include_model_id
                else None
            )
            rows.append(
                {
                    "structure_id": structure_id,
                    "site_id": site_id,
                    "cofactor_residue_name": cofactor_names,
                    "cofactor_residue_number": cofactor_numbers,
                    "cofactor_chain": cofactor_chains,
                    "cofactor_insertion_code": cofactor_insertions,
                    "shell": shell,
                    "residue_name": residue_name,
                    "residue_number": residue_number,
                    "chain": chain,
                    "insertion_code": first.get("insertion_code", ""),
                    "hetero_flag": first.get("hetero_flag", ""),
                    "atoms_involved": ",".join(
                        sorted({str(atom.get("name", "")) for atom in group})
                    ),
                    **({"model_id": first.get("model_id", "")} if include_model_id else {}),
                    "minimum_distance_A": (
                        None
                        if shell == "Cofactor"
                        else model_distance
                        if model_distance is not None
                        else distance_by_residue.get(
                            (
                                str(residue_name),
                                residue_number,
                                str(chain),
                                str(first.get("insertion_code", "") or "").strip(),
                                str(first.get("hetero_flag", "") or "").strip(),
                            )
                        )
                    ),
                }
            )
    return pd.DataFrame(
        rows, columns=RESIDUE_COLUMNS_WITH_MODEL if include_model_id else RESIDUE_COLUMNS
    )


def analyze_structure(
    structure_path: PathLike,
    cofactor_resname: Union[str, Sequence[str]],
    *,
    distance_cutoff: float = 3.6,
    expand_residues: bool = False,
    combinatorial: bool = False,
    combinatorial_cofactor_cutoff: float = 20.0,
    cofactor_resname2: Optional[Union[str, Sequence[str]]] = None,
    exclude_moieties: Optional[Sequence[str]] = None,
    first_model_only: bool = False,
    shells: int = 2,
    site_mode: str = "union",
    site_model_mode: str = "pooled",
    include_carbon_seeds: bool = False,
    direct_coordination: bool = False,
    direct_coordination_cutoff: float = 2.6,
    cofactor_class_cutoffs: Optional[Dict[str, Any]] = None,
    include_model_id: bool = False,
) -> Dict[str, pd.DataFrame]:
    """Analyze one structure without creating files or plots.

    Returns a dictionary with three pandas DataFrames:

    ``residues``
        One row per residue and shell. ``site_id`` is ``all`` in union mode or
        a stable ``site_N`` identifier in per-site mode.
    ``atoms``
        Atom-level members of the cofactor, PCS, and SCS sets.
    ``links``
        The same nearest-atom links used by the legacy ``Coord_Links.csv``.
    ``contacts``
        Exact atom-pair contacts between adjacent shells, with chemical motif
        labels and direct-coordination annotations.

    ``site_model_mode`` controls per-site boundaries: ``pooled`` retains one
    residue identity across all scanned models, while ``per-model`` keeps
    otherwise identical cofactor residues separate by model. The default is
    pooled for legacy compatibility; ``first_model_only`` remains the way to
    restrict the analysis to one model.

    ``include_model_id`` is opt-in: the ``residues`` table then gains a
    ``model_id`` column and has one row per residue per model, with minimum
    distances taken from that model's contacts. The default schema is
    unchanged.
    """
    path = Path(structure_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Structure file not found: {path}")

    structure, _ = unpack_pdb_file(str(path))
    cofactor_names = _as_names(cofactor_resname)
    cofactor_names2 = _as_names(cofactor_resname2 or [])
    exclusions = list(exclude_moieties or [])
    effective_cutoff = _effective_distance_cutoff(
        cofactor_names,
        distance_cutoff,
        cofactor_class_cutoffs,
    )
    if site_mode not in {"union", "per-site"}:
        raise ValueError("site_mode must be 'union' or 'per-site'")
    if site_model_mode not in {"pooled", "per-model"}:
        raise ValueError("site_model_mode must be 'pooled' or 'per-model'")
    if shells < 1:
        raise ValueError("shells must be at least 1")

    structure_id = _structure_id(path)

    def analyze_site(
        site_id: str,
        site_keys: Optional[set] = None,
    ) -> Dict[str, pd.DataFrame]:
        if site_keys is not None:
            cofactor_atoms, shell_sets, link_rows = identify_coordination_shells(
                structure=structure,
                cofactor_resname=cofactor_names,
                distance_cutoff=effective_cutoff,
                expand_residues=expand_residues,
                combinatorial_mode=combinatorial,
                combinatorial_cofactor_cutoff=combinatorial_cofactor_cutoff,
                cofactor_resname2=cofactor_names2 or None,
                exclude_moieties=exclusions,
                output_dir=None,
                output_prefix="",
                write_coord_links=False,
                include_link_identity=True,
                first_model_only=first_model_only,
                shells=shells,
                cofactor_site_keys=site_keys,
                include_carbon_seeds=include_carbon_seeds,
            )
        elif shells == 2:
            cofactor_atoms, pcs_atoms, scs_atoms = identify_coordination_network(
                structure=structure,
                cofactor_resname=cofactor_names,
                distance_cutoff=effective_cutoff,
                expand_residues=expand_residues,
                combinatorial_mode=combinatorial,
                combinatorial_cofactor_cutoff=combinatorial_cofactor_cutoff,
                cofactor_resname2=cofactor_names2 or None,
                exclude_moieties=exclusions,
                output_dir=None,
                output_prefix="",
                write_coord_links=False,
                first_model_only=first_model_only,
                include_carbon_seeds=include_carbon_seeds,
            )
            shell_sets = {1: pcs_atoms, 2: scs_atoms}
            link_rows = _coord_link_rows(
                cofactor_atoms,
                pcs_atoms,
                scs_atoms,
                include_identity=True,
            )
        else:
            cofactor_atoms, shell_sets, link_rows = identify_coordination_shells(
                structure=structure,
                cofactor_resname=cofactor_names,
                distance_cutoff=effective_cutoff,
                expand_residues=expand_residues,
                combinatorial_mode=combinatorial,
                combinatorial_cofactor_cutoff=combinatorial_cofactor_cutoff,
                cofactor_resname2=cofactor_names2 or None,
                exclude_moieties=exclusions,
                output_dir=None,
                output_prefix="",
                write_coord_links=False,
                include_link_identity=True,
                first_model_only=first_model_only,
                shells=shells,
                include_carbon_seeds=include_carbon_seeds,
            )

        _annotate_direct_coordination(
            link_rows,
            enabled=direct_coordination,
            cutoff=direct_coordination_cutoff,
        )
        _annotate_link_motifs(link_rows)
        shell_tables = [("Cofactor", cofactor_atoms)] + [
            (_shell_label(number), shell_sets[number])
            for number in sorted(shell_sets)
        ]
        link_columns = LINK_COLUMNS if site_id == "all" else ["structure_id", "site_id"] + COORD_LINK_COLUMNS + LINK_IDENTITY_COLUMNS
        links = pd.DataFrame(
            [
                ({"structure_id": structure_id, **({} if site_id == "all" else {"site_id": site_id}), **row})
                for row in link_rows
            ],
            columns=link_columns,
        )
        atoms = _atom_rows(structure_id, site_id, shell_tables)
        contacts = _motif_contact_rows(
            structure_id,
            site_id,
            cofactor_atoms,
            shell_sets,
            distance_cutoff=effective_cutoff,
            direct_coordination_cutoff=direct_coordination_cutoff,
        )
        residues = _residue_rows(
            structure_id,
            site_id,
            cofactor_atoms,
            shell_tables,
            link_rows,
            include_model_id=include_model_id,
            contacts=contacts,
        )
        atoms = _annotate_atom_roles(atoms, contacts)
        return {"residues": residues, "atoms": atoms, "links": links, "contacts": contacts}

    if site_mode == "union":
        return analyze_site("all")

    site_groups = _cofactor_site_groups(
        structure,
        cofactor_names,
        cofactor_names2,
        combinatorial=combinatorial,
        cutoff=combinatorial_cofactor_cutoff,
        first_model_only=first_model_only,
        model_mode=site_model_mode,
    )
    site_results = [
        analyze_site(_site_label(index, group, site_model_mode), group)
        for index, group in enumerate(site_groups, start=1)
    ]
    if not site_results:
        return {
            "residues": pd.DataFrame(
                columns=RESIDUE_COLUMNS_WITH_MODEL if include_model_id else RESIDUE_COLUMNS
            ),
            "atoms": pd.DataFrame(columns=ATOM_COLUMNS),
            "links": pd.DataFrame(columns=["structure_id", "site_id"] + COORD_LINK_COLUMNS + LINK_IDENTITY_COLUMNS),
            "contacts": pd.DataFrame(columns=CONTACT_COLUMNS),
        }
    return {
        key: pd.concat([result[key] for result in site_results], ignore_index=True)
        for key in ("residues", "atoms", "links", "contacts")
    }


def discover_structure_paths(inputs: Union[PathLike, Iterable[PathLike]]) -> List[Path]:
    """Expand files/directories into a sorted list of PDB/mmCIF paths."""
    if isinstance(inputs, (str, Path)):
        values: Iterable[PathLike] = [inputs]
    else:
        values = inputs

    allowed_suffixes = {".pdb", ".cif", ".mmcif", ".gz"}
    paths: List[Path] = []
    for value in values:
        path = Path(value).expanduser()
        if path.is_dir():
            paths.extend(
                candidate
                for candidate in path.rglob("*")
                if candidate.is_file()
                and candidate.suffix.lower() in allowed_suffixes
                and (
                    candidate.suffix.lower() != ".gz"
                    or candidate.with_suffix("").suffix.lower() in {".pdb", ".cif", ".mmcif"}
                )
            )
        elif path.is_file():
            paths.append(path)
        else:
            # Keep explicit missing inputs in the work list. The batch worker
            # records the error for this structure while allowing other files
            # to complete.
            paths.append(path)

    return sorted(set(path.resolve() for path in paths), key=str)


def _batch_worker(payload: Tuple[str, Dict[str, Any]]) -> Tuple[str, Optional[Dict[str, pd.DataFrame]], Optional[str]]:
    path, options = payload
    try:
        return path, analyze_structure(path, **options), None
    except Exception as exc:  # per-structure isolation is part of the API
        return path, None, f"{type(exc).__name__}: {exc}"


def batch_analyze(
    inputs: Union[PathLike, Iterable[PathLike]],
    cofactor_resname: Union[str, Sequence[str]],
    *,
    workers: int = 1,
    **options: Any,
) -> Dict[str, pd.DataFrame]:
    """Analyze many structures and combine tables without aborting on errors."""
    if workers < 1:
        raise ValueError("workers must be at least 1")

    paths = discover_structure_paths(inputs)
    worker_options = {"cofactor_resname": cofactor_resname, **options}
    payloads = [(str(path), worker_options) for path in paths]
    results: List[Tuple[str, Optional[Dict[str, pd.DataFrame]], Optional[str]]] = []

    if workers == 1:
        results = [_batch_worker(payload) for payload in payloads]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(_batch_worker, payload) for payload in payloads]
            for future in as_completed(futures):
                results.append(future.result())
        results.sort(key=lambda result: result[0])

    residue_tables = [result[1]["residues"] for result in results if result[1] is not None]
    atom_tables = [result[1]["atoms"] for result in results if result[1] is not None]
    link_tables = [result[1]["links"] for result in results if result[1] is not None]
    contact_tables = [result[1]["contacts"] for result in results if result[1] is not None]
    errors = [
        {
            "structure_id": _structure_id(Path(path)),
            "structure_path": path,
            "error": error,
        }
        for path, result, error in results
        if result is None and error is not None
    ]

    return {
        "residues": pd.concat(residue_tables, ignore_index=True)
        if residue_tables
        else pd.DataFrame(
            columns=RESIDUE_COLUMNS_WITH_MODEL
            if options.get("include_model_id")
            else RESIDUE_COLUMNS
        ),
        "atoms": pd.concat(atom_tables, ignore_index=True)
        if atom_tables
        else pd.DataFrame(columns=ATOM_COLUMNS),
        "links": pd.concat(link_tables, ignore_index=True)
        if link_tables
        else pd.DataFrame(columns=LINK_COLUMNS),
        "contacts": pd.concat(contact_tables, ignore_index=True)
        if contact_tables
        else pd.DataFrame(columns=CONTACT_COLUMNS),
        "errors": pd.DataFrame(errors, columns=["structure_id", "structure_path", "error"]),
    }
