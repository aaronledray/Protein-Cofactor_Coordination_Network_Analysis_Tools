"""Assembly and symmetry policy helpers for labeling workflows.

The analysis always works on the deposited asymmetric unit: biological-assembly
(BIOMT / ``_pdbx_struct_assembly``) transforms and crystal-symmetry mates are
never applied, so contacts across symmetry-related chains are not seen. Copies
that the file itself contains (for example non-crystallographic duplicates)
are separate sites. These helpers record that context and group chemically
equivalent sites so train/test splits can keep copies together.
"""

import hashlib
from pathlib import Path
from typing import Any, Dict, Union

import pandas as pd

ASSEMBLY_POLICY = "asymmetric_unit_only"
_IDENTITY = (1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)


def _pdb_context(path: Path) -> Dict[str, Any]:
    space_group = ""
    biomolecules = set()
    operators: Dict[Any, Dict[int, tuple]] = {}
    with open(path, errors="replace") as handle:
        for line in handle:
            if line.startswith("CRYST1"):
                space_group = line[55:66].strip()
            elif line.startswith("REMARK 350 BIOMOLECULE:"):
                biomolecules.add(line.split(":")[1].strip())
                current = line.split(":")[1].strip()
            elif line.startswith("REMARK 350") and "BIOMT" in line:
                parts = line.split()
                try:
                    row, serial = int(parts[2][5:]), int(parts[3])
                    values = tuple(float(v) for v in parts[4:8])
                except (IndexError, ValueError):
                    continue
                operators.setdefault((current, serial), {})[row] = values
    non_identity = sum(
        1
        for rows in operators.values()
        if len(rows) == 3 and tuple(v for r in (1, 2, 3) for v in rows[r]) != _IDENTITY
    )
    return {
        "space_group": space_group,
        "n_assemblies": len(biomolecules),
        "non_identity_operators": non_identity,
    }


def _cif_context(path: Path) -> Dict[str, Any]:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    data = MMCIF2Dict(str(path))

    def as_list(key):
        value = data.get(key, [])
        return [value] if isinstance(value, str) else list(value)

    types = as_list("_pdbx_struct_oper_list.type")
    group = as_list("_symmetry.space_group_name_H-M")
    return {
        "space_group": group[0].strip() if group else "",
        "n_assemblies": len(as_list("_pdbx_struct_assembly.id")),
        "non_identity_operators": sum(1 for t in types if "identity" not in t.lower()),
    }


def symmetry_context(structure_path: Union[str, Path]) -> Dict[str, Any]:
    """Describe the symmetry information a file carries but analysis ignores."""
    path = Path(structure_path)
    is_cif = path.suffix.lower() in {".cif", ".mmcif"}
    context = _cif_context(path) if is_cif else _pdb_context(path)
    context["assembly_policy"] = ASSEMBLY_POLICY
    context["has_unapplied_transforms"] = context["non_identity_operators"] > 0
    return context


def equivalent_site_groups(labels: pd.DataFrame, level: str = "coordination") -> pd.DataFrame:
    """Group sites with identical labeled chemistry, for leakage-aware splits.

    ``level="coordination"`` (default) signs a site by the sorted residue names
    of its primary coordinators, ignoring chain, numbering, structure, and
    site IDs. It is deliberately coarse: it groups symmetry/NCS copies and
    near-duplicate entries, but also distinct sites with the same coordination
    chemistry (for example two Cys4 Fe-S clusters). That errs toward keeping
    related sites together. ``level="full"`` also uses every shell's residue
    names, motifs, atoms, and flags, which is strict enough that real NCS
    copies often differ through small contact noise; use it only to detect
    near-exact duplicates. Grouping never proves symmetry relatedness.
    """
    if level not in {"coordination", "full"}:
        raise ValueError("level must be 'coordination' or 'full'")
    columns = ["structure_id", "site_id", "model_id", "level", "signature",
               "equivalence_group", "group_size"]
    if labels.empty:
        return pd.DataFrame(columns=columns)

    def parts_for(group: pd.DataFrame):
        if level == "coordination":
            coordinators = group[group["coordination_role"] == "primary_coordinator"]
            return [str(name) for name in coordinators["residue_name"]]
        return [
            "|".join(map(str, (r.shell_depth, r.residue_name, r.motifs, r.atoms_involved,
                               bool(r.is_primary), bool(r.direct_coordination))))
            for r in group.itertuples()
        ]

    rows = []
    for key, group in labels.groupby(["structure_id", "site_id", "model_id"], sort=True, dropna=False):
        digest = hashlib.sha1("\n".join(sorted(parts_for(group))).encode()).hexdigest()[:12]
        rows.append({"structure_id": key[0], "site_id": key[1], "model_id": key[2],
                     "level": level, "signature": digest,
                     "equivalence_group": f"eq_{digest[:8]}"})
    frame = pd.DataFrame(rows)
    frame["group_size"] = frame.groupby("signature")["signature"].transform("size")
    return frame[columns]
