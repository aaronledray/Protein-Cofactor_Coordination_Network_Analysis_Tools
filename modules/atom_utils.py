# coordination/atom_utils.py

from collections import defaultdict
from Bio.PDB import Residue
import plotly.graph_objects as go

def extract_coordinates(residue, atom_names):
    """Extract specified atom coordinates from a Biopython Residue."""
    atom_coords = []
    for atom in residue:
        if atom.get_name() in atom_names:
            atom_coords.append((atom.get_name(), atom.get_coord()))
    return atom_coords

def get_residue_atoms(residue):
    """Return all atoms in a residue as a list of dictionaries with atom details."""
    residue_number = residue.get_id()[1]
    residue_name = residue.get_resname()
    chain_id = residue.get_full_id()[2]
    return [
        {
            'coordinates': atom.coord,
            'name': atom.get_name(),
            'element': atom.element,
            'residue': residue_name,
            'residue_number': residue_number,
            'chain': chain_id
        }
        for atom in residue
    ]

def get_residue_bonds(residue_atoms, bond_lookup):
    """
    Return bonds (pairs of atom coordinates) for a single residue using the bond lookup table.
    """
    if not isinstance(residue_atoms, list) or not residue_atoms:
        print("[ERROR] Invalid residue_atoms structure.")
        return []

    if not isinstance(residue_atoms[0], dict):
        print("[ERROR] residue_atoms does not contain dictionaries.")
        return []

    residue_name = residue_atoms[0].get("residue")
    if residue_name in ["HOH", "NA", "WAT"]:
        return []

    if residue_name not in bond_lookup:
        if residue_name not in ["HOH", "NA", "WAT"]:
            print(f"[WARNING] Residue {residue_name} not found in bond lookup table.")
        return []

    try:
        atom_dict = {atom["name"]: atom["coordinates"] for atom in residue_atoms}
    except KeyError as e:
        print(f"[ERROR] Missing key in atom dictionary: {e}")
        return []

    bonds = []
    for bond in bond_lookup[residue_name]:
        if bond[0] in atom_dict and bond[1] in atom_dict:
            bonds.append((atom_dict[bond[0]], atom_dict[bond[1]]))
    return bonds

def add_bonds_to_plotly(fig, bonds):
    """Add bonds as lines to a Plotly 3D figure."""
    for bond in bonds:
        x_coords, y_coords, z_coords = zip(*bond)
        fig.add_trace(go.Scatter3d(
            x=x_coords, y=y_coords, z=z_coords,
            mode='lines',
            line=dict(color='black', width=2),
            name='Bond'
        ))

def add_bonds_to_matplotlib(ax, bonds):
    """Add bonds as lines to a Matplotlib 3D axis."""
    for bond in bonds:
        x_coords, y_coords, z_coords = zip(*bond)
        ax.plot(x_coords, y_coords, z_coords, color='black', linewidth=1)

def generate_residue_bonds(residues, bond_lookup):
    """
    Generate bonds for a list of atoms grouped into residues using bond_lookup.
    """
    bonds = []

    if not isinstance(residues, list):
        print("[ERROR] residue_atoms must be a list.")
        return bonds

    for atom in residues:
        if not isinstance(atom, dict) or 'residue' not in atom or 'residue_number' not in atom:
            print(f"[ERROR] Invalid atom entry: {atom}")
            return bonds

    residue_groups = defaultdict(list)
    for atom in residues:
        key = (atom['residue'], atom['residue_number'])
        residue_groups[key].append(atom)

    for (residue_name, _), atoms in residue_groups.items():
        if residue_name in bond_lookup:
            atom_dict = {atom['name']: atom['coordinates'] for atom in atoms}
            for bond in bond_lookup[residue_name]:
                if bond[0] in atom_dict and bond[1] in atom_dict:
                    bonds.append((atom_dict[bond[0]], atom_dict[bond[1]]))
        else:
            print(f"[WARNING] Residue {residue_name} not in bond lookup table.")

    return bonds

def generate_all_bonds(structure, bond_lookup):
    """Generate all bonds for a full Biopython structure object."""
    all_bonds = []
    for residue in structure.get_residues():
        residue_name = residue.get_resname()
        if residue_name is None:
            print("[WARNING] Residue name is None. Skipping.")
            continue

        atoms = [{'name': atom.get_name(),
                  'element': atom.element,
                  'coordinates': atom.coord,
                  'residue': residue_name}
                 for atom in residue]

        if not atoms:
            print(f"[WARNING] Residue {residue_name} has no atoms. Skipping.")
            continue

        all_bonds.extend(get_residue_bonds(atoms, bond_lookup))

    return all_bonds
