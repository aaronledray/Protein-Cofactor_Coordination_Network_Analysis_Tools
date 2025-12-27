# coordination/constants.py

"""
This module defines constants and lookup tables used throughout the
coordination fingerprint analysis project.
"""

import warnings
from Bio import BiopythonWarning



# Suppress warnings globally (e.g., PDBConstructionWarning)
warnings.simplefilter('ignore', BiopythonWarning)



# Atom color scheme for consistent plotting and 3D rendering
atom_type_colors = {
    "C": "black",
    "N": "blue",
    "O": "red",
    "S": "yellow",
    "FE": "orange",
    "MN": "purple",
    "CA": "green",
    "CU": "goldenrod"
}















bond_lookup = {
    
    
    
    
    # New entries for non-standard residues
    "CU": [],           # Residue CU: no bonds defined (or define as needed)
    "OXY": [("O1", "O2")],  # Residue OXY: bond between O1 and O2



    # Standard amino acids
    "ALA": [("N", "CA"), ("CA", "C"), ("C", "O")],
    "GLY": [("N", "CA"), ("CA", "C"), ("C", "O")],
    "ASP": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2"), ("CA", "C"), ("C", "O")],
    "GLU": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2"), ("CA", "C"), ("C", "O")],
    "HIS": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"), ("CD2", "NE2"), ("CA", "C"), ("C", "O"), ("CE1", "NE2")],
    "HIE": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"), ("CD2", "NE2"), ("CA", "C"), ("C", "O"), ("CE1", "NE2")],
    "HD1": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"), ("CD2", "NE2"), ("CA", "C"), ("C", "O"), ("CE1", "NE2")],
    "CYS": [("N", "CA"), ("CA", "CB"), ("CB", "SG"), ("CA", "C"), ("C", "O")],
    "LYS": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ"), ("CA", "C"), ("C", "O")],
    "ARG": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"), ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2"), ("CA", "C"), ("C", "O")],
    "SER": [("N", "CA"), ("CA", "CB"), ("CB", "OG"), ("CA", "C"), ("C", "O")],
    "THR": [("N", "CA"), ("CA", "CB"), ("CB", "OG1"), ("CB", "CG2"), ("CA", "C"), ("C", "O")],
    "TYR": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CZ", "OH"), ("CA", "C"), ("C", "O"), ("CZ", "CE2")],
    "PHE": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CA", "C"), ("C", "O")],
    "TRP": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CD1", "NE1"), ("CG", "CD2"), ("CD2", "CE2"), ("CD2", "CE3"), ("CE2", "CZ2"), ("CE3", "CZ3"), ("CZ2", "CH2"), ("CA", "C"), ("C", "O")],
    "VAL": [("N", "CA"), ("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CA", "C"), ("C", "O")],
    "LEU": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CA", "C"), ("C", "O")],
    "ILE": [("N", "CA"), ("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1"), ("CA", "C"), ("C", "O")],
    "PRO": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "N")],
    "MET": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "SD"), ("SD", "CE")],
    "ASN": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2"), ("CA", "C"), ("C", "O")],
    "GLN": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2"), ("CA", "C"), ("C", "O")],

    #SO3 (sulf)
    "SO3": [("S", "O1"), ("S", "O2"), ("S", "O3")],
    
    
    
    # HEM (HemeB) S
    "HEM": [
            ("CHA", "C1A"), ("C1A", "C2A"), ("CMA", "CDA"), ("CDA", "CEA"),
            ("CEA", "O3A"), ("CEA", "O4A"), ("C2A", "C3A"), ("CAA", "CBA"), ("C3A", "C4A"),
            ("CBA", "CCA"), ("CCA", "O1A"), ("CCA", "O2A"), ("C4A", "NA"), ("C1A", "NA"),
            ("CHB", "C1B"), ("C1B", "C2B"), ("C2B", "CMB"), ("CMB", "CDB"), ("CDB", "CEB"),
            ("CEB", "O3B"), ("CEB", "O4B"), ("C2B", "C3B"), ("C3B", "CAB"), ("CAB", "CBB"),
            ("CBB", "CCB"), ("CCB", "O1B"), ("CCB", "O2B"), ("C4B", "NB"), ("C1B", "NB"), ("C3B", "C4B"),
            ("CHC", "C1C"), ("C1C", "C2C"), ("C2C", "CDC"), ("CDC", "CEC"), ("CEC", "O3C"),
            ("CEC", "O4C"), ("C2C", "C3C"), ("C3C", "CAC"), ("CAC", "CBC"), ("CBC", "CCC"),
            ("CCC", "O1C"), ("CCC", "O2C"), ("C4C", "NC"), ("C1C", "NC"), ("C3C", "C4C"),
            ("CHD", "C1D"), ("C1D", "C2D"), ("C2D", "CAD"), ("CAD", "CBD"), ("CBD", "CCD"),
            ("CCD", "O1D"), ("CCD", "O2D"), ("C2D", "C3D"), ("C3D", "CDD"), ("CDD", "CED"),
            ("CED", "O3D"), ("CED", "O4D"), ("C4D", "ND"), ("C1D", "ND"), ("C3D", "C4D"),
            ("FE", "NA"), ("FE", "NB"), ("FE", "NC"), ("FE", "ND"),  # FE-NA bonds / Central iron and surrounding rings
            ("CHA", "C4D"),  # CHA to C4D
            ("CHD", "C4C"),  # CHD to C4C
            ("CHC", "C4B"),  # CHC to C4B
            ("CHB", "C4A"),   # CHB to C4A
            ("CMD", "C2D"),
            ("C2A", "CAA"),
            ("CMC", "C2C"),
            ("CBD", "CGD"),
            ("CGD", "O2D"),
            ("CGD", "O1D"),
            ("CBA", "CGA"),
            ("C3A", "CMA"),
            ("O2A", "CGA"),
            ("CGA", "O1A")],
          
    # HEM (HemeB) S
    "HEA": [
            ("CHA", "C1A"), ("C1A", "C2A"), ("CMA", "CDA"), ("CDA", "CEA"),
            ("CEA", "O3A"), ("CEA", "O4A"), ("C2A", "C3A"), ("CAA", "CBA"), ("C3A", "C4A"),
            ("CBA", "CCA"), ("CCA", "O1A"), ("CCA", "O2A"), ("C4A", "NA"), ("C1A", "NA"),
            ("CHB", "C1B"), ("C1B", "C2B"), ("C2B", "CMB"), ("CMB", "CDB"), ("CDB", "CEB"),
            ("CEB", "O3B"), ("CEB", "O4B"), ("C2B", "C3B"), ("C3B", "CAB"), ("CAB", "CBB"),
            ("CBB", "CCB"), ("CCB", "O1B"), ("CCB", "O2B"), ("C4B", "NB"), ("C1B", "NB"), ("C3B", "C4B"),
            ("CHC", "C1C"), ("C1C", "C2C"), ("C2C", "CDC"), ("CDC", "CEC"), ("CEC", "O3C"),
            ("CEC", "O4C"), ("C2C", "C3C"), ("C3C", "CAC"), ("CAC", "CBC"), ("CBC", "CCC"),
            ("CCC", "O1C"), ("CCC", "O2C"), ("C4C", "NC"), ("C1C", "NC"), ("C3C", "C4C"),
            ("CHD", "C1D"), ("C1D", "C2D"), ("C2D", "CAD"), ("CAD", "CBD"), ("CBD", "CCD"),
            ("CCD", "O1D"), ("CCD", "O2D"), ("C2D", "C3D"), ("C3D", "CDD"), ("CDD", "CED"),
            ("CED", "O3D"), ("CED", "O4D"), ("C4D", "ND"), ("C1D", "ND"), ("C3D", "C4D"),
            ("FE", "NA"), ("FE", "NB"), ("FE", "NC"), ("FE", "ND"),  # FE-NA bonds / Central iron and surrounding rings
            ("CHA", "C4D"),  # CHA to C4D
            ("CHD", "C4C"),  # CHD to C4C
            ("CHC", "C4B"),  # CHC to C4B
            ("CHB", "C4A"),   # CHB to C4A
            ("CMD", "C2D"),
            ("C2A", "CAA"),
            ("CMC", "C2C"),
            ("CBD", "CGD"),
            ("CGD", "O2D"),
            ("CGD", "O1D"),
            ("CBA", "CGA"),
            ("C3A", "CMA"),
            ("O2A", "CGA"),
            ("CGA", "O1A"), ("C3B", "C11"), ("C11", "C12"), ("C11","O11"), ("CMA","OMA")],
    
    
    
    
    
    
    # HEM (HemeB) S
    "HM1": [
            ("CHA", "C1A"), ("C1A", "C2A"), ("CMA", "CDA"), ("CDA", "CEA"),
            ("CEA", "O3A"), ("CEA", "O4A"), ("C2A", "C3A"), ("CAA", "CBA"), ("C3A", "C4A"),
            ("CBA", "CCA"), ("CCA", "O1A"), ("CCA", "O2A"), ("C4A", "NA"), ("C1A", "NA"),
            ("CHB", "C1B"), ("C1B", "C2B"), ("C2B", "CMB"), ("CMB", "CDB"), ("CDB", "CEB"),
            ("CEB", "O3B"), ("CEB", "O4B"), ("C2B", "C3B"), ("C3B", "CAB"), ("CAB", "CBB"),
            ("CBB", "CCB"), ("CCB", "O1B"), ("CCB", "O2B"), ("C4B", "NB"), ("C1B", "NB"), ("C3B", "C4B"),
            ("CHC", "C1C"), ("C1C", "C2C"), ("C2C", "CDC"), ("CDC", "CEC"), ("CEC", "O3C"),
            ("CEC", "O4C"), ("C2C", "C3C"), ("C3C", "CAC"), ("CAC", "CBC"), ("CBC", "CCC"),
            ("CCC", "O1C"), ("CCC", "O2C"), ("C4C", "NC"), ("C1C", "NC"), ("C3C", "C4C"),
            ("CHD", "C1D"), ("C1D", "C2D"), ("C2D", "CAD"), ("CAD", "CBD"), ("CBD", "CCD"),
            ("CCD", "O1D"), ("CCD", "O2D"), ("C2D", "C3D"), ("C3D", "CDD"), ("CDD", "CED"),
            ("CED", "O3D"), ("CED", "O4D"), ("C4D", "ND"), ("C1D", "ND"), ("C3D", "C4D"),
            ("FE", "NA"), ("FE", "NB"), ("FE", "NC"), ("FE", "ND"),  # FE-NA bonds / Central iron and surrounding rings
            ("CHA", "C4D"),  # CHA to C4D
            ("CHD", "C4C"),  # CHD to C4C
            ("CHC", "C4B"),  # CHC to C4B
            ("CHB", "C4A"),   # CHB to C4A
            ("CMD", "C2D"),
            ("C2A", "CAA"),
            ("CMC", "C2C"),
            ("CBD", "CGD"),
            ("CGD", "O2D"),
            ("CGD", "O1D"),
            ("CBA", "CGA"),
            ("C3A", "CMA"),
            ("O2A", "CGA"),
            ("CGA", "O1A"), ("C3B", "C11"), ("C11", "C12"), ("C11","O11"), ("CMA","OMA")],
          
    
        
    
    # HEM (HemeB) S
    "HM2": [
            ("CHA", "C1A"), ("C1A", "C2A"), ("CMA", "CDA"), ("CDA", "CEA"),
            ("CEA", "O3A"), ("CEA", "O4A"), ("C2A", "C3A"), ("CAA", "CBA"), ("C3A", "C4A"),
            ("CBA", "CCA"), ("CCA", "O1A"), ("CCA", "O2A"), ("C4A", "NA"), ("C1A", "NA"),
            ("CHB", "C1B"), ("C1B", "C2B"), ("C2B", "CMB"), ("CMB", "CDB"), ("CDB", "CEB"),
            ("CEB", "O3B"), ("CEB", "O4B"), ("C2B", "C3B"), ("C3B", "CAB"), ("CAB", "CBB"),
            ("CBB", "CCB"), ("CCB", "O1B"), ("CCB", "O2B"), ("C4B", "NB"), ("C1B", "NB"), ("C3B", "C4B"),
            ("CHC", "C1C"), ("C1C", "C2C"), ("C2C", "CDC"), ("CDC", "CEC"), ("CEC", "O3C"),
            ("CEC", "O4C"), ("C2C", "C3C"), ("C3C", "CAC"), ("CAC", "CBC"), ("CBC", "CCC"),
            ("CCC", "O1C"), ("CCC", "O2C"), ("C4C", "NC"), ("C1C", "NC"), ("C3C", "C4C"),
            ("CHD", "C1D"), ("C1D", "C2D"), ("C2D", "CAD"), ("CAD", "CBD"), ("CBD", "CCD"),
            ("CCD", "O1D"), ("CCD", "O2D"), ("C2D", "C3D"), ("C3D", "CDD"), ("CDD", "CED"),
            ("CED", "O3D"), ("CED", "O4D"), ("C4D", "ND"), ("C1D", "ND"), ("C3D", "C4D"),
            ("FE", "NA"), ("FE", "NB"), ("FE", "NC"), ("FE", "ND"),  # FE-NA bonds / Central iron and surrounding rings
            ("CHA", "C4D"),  # CHA to C4D
            ("CHD", "C4C"),  # CHD to C4C
            ("CHC", "C4B"),  # CHC to C4B
            ("CHB", "C4A"),   # CHB to C4A
            ("CMD", "C2D"),
            ("C2A", "CAA"),
            ("CMC", "C2C"),
            ("CBD", "CGD"),
            ("CGD", "O2D"),
            ("CGD", "O1D"),
            ("CBA", "CGA"),
            ("C3A", "CMA"),
            ("O2A", "CGA"),
            ("CGA", "O1A"), ("C3B", "C11"), ("C11", "C12"), ("C11","O11"), ("CMA","OMA")],
          
    
    
    
    
    # SRM (siroheme) S
    "SRM": [
            ("CHA", "C1A"), ("C1A", "C2A"), ("C2A", "CMA"), ("CMA", "CDA"), ("CDA", "CEA"),
            ("CEA", "O3A"), ("CEA", "O4A"), ("C2A", "C3A"), ("C3A", "CAA"), ("CAA", "CBA"), ("C3A", "C4A"),
            ("CBA", "CCA"), ("CCA", "O1A"), ("CCA", "O2A"), ("C4A", "NA"), ("C1A", "NA"),
            ("CHB", "C1B"), ("C1B", "C2B"), ("C2B", "CMB"), ("CMB", "CDB"), ("CDB", "CEB"),
            ("CEB", "O3B"), ("CEB", "O4B"), ("C2B", "C3B"), ("C3B", "CAB"), ("CAB", "CBB"),
            ("CBB", "CCB"), ("CCB", "O1B"), ("CCB", "O2B"), ("C4B", "NB"), ("C1B", "NB"), ("C3B", "C4B"),
            ("CHC", "C1C"), ("C1C", "C2C"), ("C2C", "CDC"), ("CDC", "CEC"), ("CEC", "O3C"),
            ("CEC", "O4C"), ("C2C", "C3C"), ("C3C", "CAC"), ("CAC", "CBC"), ("CBC", "CCC"),
            ("CCC", "O1C"), ("CCC", "O2C"), ("C4C", "NC"), ("C1C", "NC"), ("C3C", "C4C"),
            ("CHD", "C1D"), ("C1D", "C2D"), ("C2D", "CAD"), ("CAD", "CBD"), ("CBD", "CCD"),
            ("CCD", "O1D"), ("CCD", "O2D"), ("C2D", "C3D"), ("C3D", "CDD"), ("CDD", "CED"),
            ("CED", "O3D"), ("CED", "O4D"), ("C4D", "ND"), ("C1D", "ND"), ("C3D", "C4D"),
            ("FE", "NA"), ("FE", "NB"), ("FE", "NC"), ("FE", "ND"),  # FE-NA bonds / Central iron and surrounding rings
            ("CHA", "C4D"),  # CHA to C4D
            ("CHD", "C4C"),  # CHD to C4C
            ("CHC", "C4B"),  # CHC to C4B
            ("CHB", "C4A")],   # CHB to C4A
            
    
    # SF4 (4Fe4S) and SO3 (sulf)
    "SF4": [("FE1", "S1"), ("FE2", "S1"), ("FE1", "S2"), ("FE3", "S2"), ("FE2", "S3"), ("FE4", "S3"), ("FE3", "S4"),
            ("FE4", "S1"), 
            ("FE4", "S2"), 
            ("FE1", "S3"), 
            ("FE2", "S4"), 
            ("FE3", "S1")],
            
    
    
    
    # Oxygen evolving cluster:
    "OEX": [("MN4", "O4"),
            ("MN4", "O5"),
            ("MN3", "O2"),
            ("MN3", "O5"),
            ("MN3", "O3"),
            ("MN3", "O4"),
            ("MN2", "O3"),
            ("MN2", "O2"),
            ("MN2", "O1"),
            ("MN1", "O5"),
            ("MN1", "O3"),
            ("MN1", "O1"),
            ("CA1", "O1"),
            ("CA1", "O2"),
            ("CA1", "O5")],
    
    
    "OY1": [("O1", "O2")],

    
}



