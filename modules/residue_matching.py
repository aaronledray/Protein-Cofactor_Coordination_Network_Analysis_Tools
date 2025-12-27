import numpy as np
from collections import defaultdict

def Residue_Match_At_Template_Position(template_residue_sphere, query_residue_sphere, distance_cutoff=3.6):
    results = []

    for template_residue in template_residue_sphere:
        template_coords = np.array(template_residue["coordinates"])
        template_residue_number = template_residue["residue_number"]
        template_residue_name = template_residue["residue_name"]

        closest_distance = float("inf")
        closest_residue_name = "None"
        closest_residue_chain = "None"
        closest_residue_moiety = "None"

        for query_residue in query_residue_sphere:
            query_coords = np.array(query_residue["coordinates"])
            distance = np.linalg.norm(template_coords - query_coords)

            if distance < closest_distance and distance <= distance_cutoff:
                closest_distance = distance
                closest_residue_name = query_residue["residue_name"]
                closest_residue_chain = query_residue.get("chain", "None")
                closest_residue_moiety = query_residue.get("moiety", "Unknown")

        results.append({
            "Template Residue Number": template_residue_number,
            "Template Residue Name": template_residue_name,
            "Query Matched Residue Name": closest_residue_name,
            "Query Chain": closest_residue_chain,
            "Distance": closest_distance if closest_residue_name != "None" else "N/A",
            "Query Matched Moiety": closest_residue_moiety,
        })

    return results

def make_residue_centroid_sphere(template_coord_residues):
    residue_groups = defaultdict(list)

    for atom in template_coord_residues:
        residue_number = atom.get("residue_number")
        residue_name = atom.get("residue")
        if residue_number is None or residue_name is None:
            print(f"[WARNING] Skipping atom with missing residue info: {atom}")
            continue
        residue_groups[residue_number].append(atom)

    cofactor_residue_sphere = []
    for residue_number, atoms in residue_groups.items():
        if not atoms:
            print(f"[WARNING] Empty atom group for residue {residue_number}.")
            continue

        coords = np.array([atom["coordinates"] for atom in atoms])
        centroid = coords.mean(axis=0)
        residue_name = atoms[0].get("residue", "Unknown")
        cofactor_residue_sphere.append({
            "residue_number": residue_number,
            "residue_name": residue_name,
            "coordinates": centroid
        })

    return cofactor_residue_sphere
