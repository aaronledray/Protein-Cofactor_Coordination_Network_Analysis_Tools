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


PathLike = Union[str, Path]

ATOM_COLUMNS = [
    "structure_id", "site_id", "shell", "residue_name", "residue_number",
    "chain", "insertion_code", "hetero_flag", "atom_name", "element",
    "x", "y", "z",
]

RESIDUE_COLUMNS = [
    "structure_id", "site_id", "cofactor_residue_name", "cofactor_residue_number",
    "cofactor_chain", "cofactor_insertion_code", "shell", "residue_name",
    "residue_number", "chain", "insertion_code", "hetero_flag",
    "atoms_involved", "minimum_distance_A",
]

LINK_IDENTITY_COLUMNS = [
    "src_insertion_code", "src_hetero_flag",
    "dst_insertion_code", "dst_hetero_flag",
    "src_element", "dst_element",
    "direct_coordination",
]
LINK_COLUMNS = ["structure_id"] + COORD_LINK_COLUMNS + LINK_IDENTITY_COLUMNS


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
) -> List[set]:
    """Return cofactor residue sites, optionally merged into nearby clusters."""
    models = list(structure)
    if first_model_only:
        models = models[:1]
    names = {name.upper() for name in (*cofactor_names, *cofactor_names2)}
    atoms_by_site: Dict[Tuple[str, Any, str, str, str], List[np.ndarray]] = {}
    for model in models:
        for chain in model:
            for residue in chain:
                if residue.get_resname().upper() not in names:
                    continue
                residue_id = residue.get_id()
                key = (
                    residue.get_resname(),
                    residue_id[1],
                    chain.id,
                    str(residue_id[2] or "").strip(),
                    str(residue_id[0] or "").strip(),
                )
                atoms_by_site.setdefault(key, []).extend(
                    np.asarray(atom.coord, dtype=float) for atom in residue
                )

    keys = sorted(atoms_by_site, key=str)
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
            if not atoms_by_site[keys[left]] or not atoms_by_site[keys[right]]:
                continue
            distances = np.linalg.norm(
                np.asarray(atoms_by_site[keys[left]])[:, None, :]
                - np.asarray(atoms_by_site[keys[right]])[None, :, :],
                axis=2,
            )
            if float(distances.min()) <= cutoff:
                union(left, right)

    grouped: Dict[int, set] = {}
    for index, key in enumerate(keys):
        grouped.setdefault(find(index), set()).add(key)
    return [grouped[root] for root in sorted(grouped)]


def _effective_distance_cutoff(
    cofactor_names: Sequence[str],
    fallback: float,
    class_cutoffs: Optional[Dict[str, float]],
) -> float:
    if not class_cutoffs:
        return fallback
    cofactor_class = (
        "metal"
        if cofactor_names and all(name.upper() in METAL_ION_RESNAMES for name in cofactor_names)
        else "organic"
    )
    return float(class_cutoffs.get(cofactor_class, fallback))


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
                    "atom_name": atom.get("name", ""),
                    "element": atom.get("element", ""),
                    "x": float(coordinates[0]),
                    "y": float(coordinates[1]),
                    "z": float(coordinates[2]),
                }
            )
    return pd.DataFrame(rows, columns=ATOM_COLUMNS)


def _residue_rows(
    structure_id: str,
    site_id: str,
    cofactor_atoms: Sequence[Dict[str, Any]],
    shell_atoms: Sequence[Tuple[str, Sequence[Dict[str, Any]]]],
    link_rows: Sequence[Dict[str, Any]],
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

    rows: List[Dict[str, Any]] = []
    for shell, atoms in shell_atoms:
        grouped: Dict[Tuple[str, Any, str], List[Dict[str, Any]]] = {}
        for atom in atoms:
            grouped.setdefault(_atom_identity(atom), []).append(atom)
        for (residue_name, residue_number, chain, insertion_code, hetero_flag), group in grouped.items():
            first = group[0]
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
                    "minimum_distance_A": (
                        None
                        if shell == "Cofactor"
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
    return pd.DataFrame(rows, columns=RESIDUE_COLUMNS)


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
    include_carbon_seeds: bool = False,
    direct_coordination: bool = False,
    direct_coordination_cutoff: float = 2.6,
    cofactor_class_cutoffs: Optional[Dict[str, float]] = None,
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
        residues = _residue_rows(
            structure_id,
            site_id,
            cofactor_atoms,
            shell_tables,
            link_rows,
        )
        return {"residues": residues, "atoms": atoms, "links": links}

    if site_mode == "union":
        return analyze_site("all")

    site_groups = _cofactor_site_groups(
        structure,
        cofactor_names,
        cofactor_names2,
        combinatorial=combinatorial,
        cutoff=combinatorial_cofactor_cutoff,
        first_model_only=first_model_only,
    )
    site_results = [
        analyze_site(f"site_{index}", group)
        for index, group in enumerate(site_groups, start=1)
    ]
    if not site_results:
        return {
            "residues": pd.DataFrame(columns=RESIDUE_COLUMNS),
            "atoms": pd.DataFrame(columns=ATOM_COLUMNS),
            "links": pd.DataFrame(columns=["structure_id", "site_id"] + COORD_LINK_COLUMNS + LINK_IDENTITY_COLUMNS),
        }
    return {
        key: pd.concat([result[key] for result in site_results], ignore_index=True)
        for key in ("residues", "atoms", "links")
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
        else pd.DataFrame(columns=RESIDUE_COLUMNS),
        "atoms": pd.concat(atom_tables, ignore_index=True)
        if atom_tables
        else pd.DataFrame(columns=ATOM_COLUMNS),
        "links": pd.concat(link_tables, ignore_index=True)
        if link_tables
        else pd.DataFrame(columns=LINK_COLUMNS),
        "errors": pd.DataFrame(errors, columns=["structure_id", "structure_path", "error"]),
    }
