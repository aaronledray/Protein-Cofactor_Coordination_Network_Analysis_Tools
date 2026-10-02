"""Profile-conservation overlays for the cohesive 3D viewer."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Optional

from .comparison_runner import analysis_from_path
from .network_alignment import reference_position_map
from .network_comparison import build_residue_conservation_map
from .io_utils import unpack_pdb_file
from .moieties import bond_lookup
from .plotting import plot_interactive_cohesive_network


def write_profile_conservation_viewer(
    reference_path: Path,
    cofactor_resname: str,
    profile_result: Mapping[str, Any],
    output_filename: Path,
    *,
    analysis_options: Optional[Mapping[str, Any]] = None,
    compact_html: bool = False,
) -> Path:
    """Render the first profile reference with per-residue family support.

    The profile's first reference is used as the coordinate/template view. Its
    residues are colored by support across the profile, while all existing
    cohesive viewer layers and contact controls remain available.
    """
    reference_records = profile_result.get("references", [])
    template_record = profile_result.get("template") or (reference_records[0] if reference_records else None)
    if template_record is None:
        raise ValueError("Profile result has no reference or template signature")
    selected_reference_path = Path(template_record.get("structure_path") or reference_path)
    signature = template_record["signature"]
    analysis = analysis_from_path(
        selected_reference_path,
        cofactor_resname,
        analysis_options=analysis_options,
    )
    residue_map = None
    if any("position" in feature for feature in signature.get("residues", [])):
        residue_map = reference_position_map(analysis)
    conservation_by_residue = build_residue_conservation_map(
        analysis,
        signature,
        profile_result["profile"],
        residue_map=residue_map,
    )

    structure, _ = unpack_pdb_file(str(selected_reference_path))
    shell_groups = {}
    for row in analysis["atoms"].to_dict("records"):
        shell_groups.setdefault(row["shell"], []).append(row)
    output_filename = Path(output_filename)
    output_filename.parent.mkdir(parents=True, exist_ok=True)
    plot_interactive_cohesive_network(
        structure=structure,
        cofactor_atoms=shell_groups.get("Cofactor", []),
        pcs_atoms=shell_groups.get("PCS", []),
        scs_atoms=shell_groups.get("SCS", []),
        focused_atoms_by_shell=shell_groups,
        contacts=analysis["contacts"],
        bond_lookup_table=bond_lookup,
        pdb_name=selected_reference_path.name,
        cofactor_resname=cofactor_resname,
        output_filename=str(output_filename),
        first_model_only=bool((analysis_options or {}).get("first_model_only", False)),
        compact_html=compact_html,
        conservation_by_residue=conservation_by_residue,
        conservation_label="Family conservation",
    )
    return output_filename


__all__ = ["write_profile_conservation_viewer"]
