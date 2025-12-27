from Bio.PDB import PDBParser

# Handling of pdb file
def unpack_pdb_file(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)

    atoms_data = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_data = {
                        'name': atom.get_name(),
                        'element': atom.element,
                        'residue': residue.get_resname(),
                        'residue_number': residue.get_id()[1],
                        'coordinates': atom.get_coord(),
                        'chain': chain.id
                    }
                    atoms_data.append(atom_data)

    return structure, atoms_data
