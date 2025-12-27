# modules/analysis.py

"""
Analysis utilities for Cofactor_Fingerprint_tool.

This module collects functions for:
- distance cutoff assessment
- CA-only and CA→CB vector evaluation
- query structure evaluation
- frequency breakdowns by residue/moiety
- high-level evaluation pipelines
"""

# modules/analysis.py




import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .structure_io import unpack_pdb_file
from .structure_processing import (
    make_residue_centroid_sphere,
    extract_query_box,
    identify_coordination_network,     # <-- add this
    identify_residues_of_interest,
    generate_coordination_csv_with_moieties,
    find_cofactor_atoms,
)
from .residue_matching import Residue_Match_At_Template_Position
from .plotting import plot_evaluation_results  # if you have this; else remove




__all__ = [
    "unpack_pdb_file",
    "make_residue_centroid_sphere",
    "extract_query_box",
    "identify_coordination_network",         # <-- add this
    "identify_residues_of_interest",
    "generate_coordination_csv_with_moieties",
    "find_cofactor_atoms",
    "Residue_Match_At_Template_Position",
    "plot_evaluation_results",
    # "evaluate_CA_only",
    # "evaluate_CA_CB_vectors",
    # "process_query_file",
    # "evaluate_query_files",
]



# --------------------------------------------------------------------
#  Distance Cutoff Assessments
# --------------------------------------------------------------------
def assess_distance_cutoff(template_residue_sphere, pdb_files, 
                           distance_min=0.0, distance_max=4.0, step=0.2):
    """Evaluate how different distance cutoffs affect residue matching."""
    cutoffs = np.arange(distance_min, distance_max + step, step)
    avg_matches, std_matches = [], []

    for cutoff in cutoffs:
        match_counts = []
        for query_pdb_file in pdb_files:
            _, query_atoms = unpack_pdb_file(query_pdb_file)
            query_residue_sphere = make_residue_centroid_sphere(query_atoms)
            results = Residue_Match_At_Template_Position(
                template_residue_sphere, query_residue_sphere, distance_cutoff=cutoff
            )
            valid_matches = sum(1 for r in results if r["Query Matched Residue Name"] != "None")
            match_counts.append(valid_matches)

        avg_matches.append(np.mean(match_counts))
        std_matches.append(np.std(match_counts))

    plt.errorbar(cutoffs, avg_matches, yerr=std_matches,
                 fmt='-', color='black', ecolor='gray', capsize=5)
    plt.xlabel("Distance Cutoff (Å)")
    plt.ylabel("Average Matches")
    plt.title("Effect of Distance Cutoff")
    plt.show()

    return cutoffs, avg_matches, std_matches


# --------------------------------------------------------------------
#  CA-only evaluation
# --------------------------------------------------------------------
def evaluate_CA_only(template_atoms, query_box, distance_cutoff):
    """Evaluate matches using only CA atoms."""
    template_CA = [a for a in template_atoms if a['name'] == 'CA']
    query_CA = [a for a in query_box if a['name'] == 'CA']

    matches, distances, matched_residues = 0, [], []

    for t_atom in template_CA:
        t_coord = np.array(t_atom['coordinates'], dtype=float)
        dists = [np.linalg.norm(t_coord - np.array(q['coordinates'], dtype=float)) for q in query_CA]
        if dists and min(dists) <= distance_cutoff:
            matches += 1
            matched_residues.append(str(t_atom['residue_number']))
            distances.append(min(dists))

    return {
        "Total_Template": len(template_CA),
        "Matches": matches,
        "Misses": len(template_CA) - matches,
        "Avg_Distance": np.mean(distances) if distances else np.nan,
        "Matched_Residues": matched_residues
    }


# --------------------------------------------------------------------
#  CA→CB vector evaluation
# --------------------------------------------------------------------
def evaluate_CA_CB_vectors(template_atoms, query_box, distance_cutoff):
    """Evaluate matches using CA atoms with CB vector orientation & length."""
    # (you pasted this one earlier — drop it here unchanged!)
    ...


# --------------------------------------------------------------------
#  Query file processing
# --------------------------------------------------------------------
def process_query_file(pdb_file, template_atoms_of_interest, distance_cutoff, mode="CA_only"):
    """Process a single query structure for evaluation."""
    _, query_atoms = unpack_pdb_file(pdb_file)
    from .structure_utils import extract_query_box  # lazy import
    query_box, _ = extract_query_box(template_atoms_of_interest, query_atoms, distance_cutoff)

    if mode == "CA_only":
        result = evaluate_CA_only(template_atoms_of_interest, query_box, distance_cutoff)
    elif mode == "CA_CB_vectors":
        result = evaluate_CA_CB_vectors(template_atoms_of_interest, query_box, distance_cutoff)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    result["Query_File"] = os.path.basename(pdb_file)
    return result


def evaluate_query_files(input_dir, template_atoms_of_interest, distance_cutoff, mode="CA_only"):
    """Loop over all query PDBs in a folder and evaluate them."""
    pdb_files = glob.glob(os.path.join(input_dir, "*.pdb"))
    return [process_query_file(f, template_atoms_of_interest, distance_cutoff, mode) for f in pdb_files]


# --------------------------------------------------------------------
#  High-level pipeline
# --------------------------------------------------------------------
def main_evaluation_pipeline(mode="CA_only", distance_cutoff=3.0):
    """Run full evaluation pipeline."""
    aligned_dir = "0_rep10_aligned/0_name_edited"
    results = evaluate_query_files(aligned_dir, template_atoms_of_interest=[],  # TODO: pass properly
                                   distance_cutoff=distance_cutoff, mode=mode)
    df = pd.DataFrame(results)
    print(df)
    out_csv = f"0_evaluation_results_{mode}.csv"
    plot_evaluation_results(df, output_csv=out_csv, mode=mode)
    return df



import csv
import math
from typing import List, Dict, Tuple, Optional

AtomDict = Dict[str, object]
