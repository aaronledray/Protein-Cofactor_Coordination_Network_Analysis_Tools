"""
chemistry.py
------------
Chemical moiety lookup table and helper functions for grouping atoms
by chemical moiety (sidechains, cofactors, backbone, water, etc).
"""

from typing import List, Dict, Optional
from Bio.PDB import Residue, Atom

# ---------------------------------------------------------
# Moieties lookup table: keys are (RESNAME, ATOMNAME), values are moiety labels.
# ---------------------------------------------------------







# chemical_moieties: Dict[tuple, str] = {
    




#     # MISSING BACKBONES:







#     # OXYGEN:
#     ("OXY", "O1"): "oxy_ligand",
#     ("OXY", "O2"): "oxy_ligand",

#     ("HEA", "C16"): "heme",
#     ("HEA", "C17"): "heme",
#     ("HEA", "C18"): "heme",
#     ("HEA", "C19"): "heme",
#     ("HEA", "C20"): "heme",
#     ("HEA", "C21"): "heme",
#     ("HEA", "C22"): "heme",
#     ("HEA", "C23"): "heme",
#     ("HEA", "C24"): "heme",
#     ("HEA", "C25"): "heme",
#     ("HEA", "C27"): "heme",

    
#     ("HIS", "CA"): "backbone",  # CA is backbone for all residues
#     ("HIS", "CB"): "histidine_sidechain",
    
    
#     ("HEA", "C13"): "heme",
#     ("HEA", "C14"): "heme",
#     ("HEA", "C15"): "heme",
#     ("HEA", "C26"): "heme",
#     ("HEA", "O1"): "heme",
#     ("HEA", "O2"): "heme",
#     ("HEA", "CAA"): "heme",
#     ("HEA", "CAC"): "heme",
#     ("HEA", "CAD"): "heme",
#     ("HEA", "NA"): "heme",
#     ("HEA", "CBA"): "heme",
#     ("HEA", "CBC"): "heme",
#     ("HEA", "CBD"): "heme",
#     ("HEA", "NB"): "heme",
#     ("HEA", "CGA"): "heme",
#     ("HEA", "CGD"): "heme",
#     ("HEA", "ND"): "heme",
#     ("HEA", "CHA"): "heme",
#     ("HEA", "CHB"): "heme",
#     ("HEA", "CHC"): "heme",
#     ("HEA", "CHD"): "heme",
#     ("HEA", "CMA"): "heme",
#     ("HEA", "CMB"): "heme",
#     ("HEA", "CMC"): "heme",
#     ("HEA", "CMD"): "heme",
#     ("HEA", "OMA"): "heme",
#     ("HEA", "C1A"): "heme",
#     ("HEA", "C1B"): "heme",
#     ("HEA", "C1C"): "heme",
#     ("HEA", "C1D"): "heme",
#     ("HEA", "O1A"): "heme",
#     ("HEA", "O1D"): "heme",
#     ("HEA", "C2A"): "heme",
#     ("HEA", "C2B"): "heme",
#     ("HEA", "C2C"): "heme",
#     ("HEA", "C2D"): "heme",
#     ("HEA", "O2A"): "heme",
#     ("HEA", "O2D"): "heme",
#     ("HEA", "C3A"): "heme",
#     ("HEA", "C3B"): "heme",
#     ("HEA", "C3C"): "heme",
#     ("HEA", "C3D"): "heme",
#     ("HEA", "C4A"): "heme",
#     ("HEA", "C4B"): "heme",
#     ("HEA", "C4C"): "heme",
#     ("HEA", "C4D"): "heme",
#     ("HEA", "C11"): "heme",
#     ("HEA", "O11"): "heme",
#     ("HEA", "C12"): "heme",
#     ("HEA", "NC"): "heme",
#     ("HEA", "FE"): "heme",
    
    
    
    
    
#     # --- CU Copper ---
#     ("CU", "CU"): "copper_center",
#     ("CU1", "CU"): "copper_center",

#     ("TRP", "CZ3"): "tryptophan_sidechain",


        
#     # --- New entries for Glutamine (GLN) ---
#     ("GLN", "H"):    "backbone",
#     ("GLN", "CA"):   "backbone",
#     ("GLN", "HA"):   "backbone",
#     ("GLN", "CB"):   "glutamine_sidechain",
#     ("GLN", "HB2"):  "glutamine_sidechain",
#     ("GLN", "HB3"):  "glutamine_sidechain",
#     ("GLN", "CG"):   "glutamine_sidechain",
#     ("GLN", "HG2"):  "glutamine_sidechain",
#     ("GLN", "HG3"):  "glutamine_sidechain",
#     ("GLN", "CD"):   "glutamine_sidechain",
#     ("GLN", "HE21"): "glutamine_sidechain",
#     ("GLN", "HE22"): "glutamine_sidechain",
    
    
    

#     # --- Additional entries for Glutamine (GLN) ---
#     ("GLN", "C"):  "backbone",
    
#     # --- Additional entries for Proline (PRO) ---
#     ("PRO", "HG3"): "proline_sidechain",  # Missing HG3 in PRO
    
#     # --- Additional entries for Tryptophan (TRP) ---
#     ("TRP", "H"):    "backbone",
#     ("TRP", "CA"):   "backbone",
#     ("TRP", "HA"):   "backbone",
#     ("TRP", "C"):    "backbone",
#     ("TRP", "CB"):   "tryptophan_sidechain",
#     ("TRP", "HB2"):  "tryptophan_sidechain",
#     ("TRP", "HB3"):  "tryptophan_sidechain",
#     ("TRP", "CG"):   "tryptophan_sidechain",
#     ("TRP", "CD1"):  "tryptophan_sidechain",
#     ("TRP", "HD1"):  "tryptophan_sidechain",
#     ("TRP", "NE1"):  "tryptophan_sidechain",
#     ("TRP", "HE1"):  "tryptophan_sidechain",
#     ("TRP", "CE2"):  "tryptophan_sidechain",
#     ("TRP", "CZ2"):  "tryptophan_sidechain",
#     ("TRP", "HZ2"):  "tryptophan_sidechain",
#     ("TRP", "CE3"):  "tryptophan_sidechain",
#     ("TRP", "HE3"):  "tryptophan_sidechain",
#     ("TRP", "CD2"):  "tryptophan_sidechain",
    
#     # --- Additional entries for Asparagine (ASN) ---
#     ("ASN", "H"):     "backbone",
#     ("ASN", "CB"):    "asparagine_sidechain",
#     ("ASN", "HB2"):   "asparagine_sidechain",
#     ("ASN", "HB3"):   "asparagine_sidechain",
#     ("ASN", "CG"):    "asparagine_sidechain",
#     ("ASN", "HD21"):  "asparagine_sidechain",
#     ("ASN", "HD22"):  "asparagine_sidechain",
    
    

        
#     # --- Additional entries for Glycine (GLY) ---
#     ("GLY", "HA2"): "backbone",
#     ("GLY", "HA3"): "backbone",
    
#     # --- Additional entries for Asparagine (ASN) ---
#     ("ASN", "C"):  "backbone",
#     ("ASN", "HA"): "backbone",
#     ("ASN", "CA"): "backbone",
    
#     # --- Additional entries for Proline (PRO) ---
#     ("PRO", "C"):   "backbone",
#     ("PRO", "CA"):  "backbone",
#     ("PRO", "HA"):  "backbone",
#     ("PRO", "CB"):  "proline_sidechain",
#     ("PRO", "HB2"): "proline_sidechain",
#     ("PRO", "HB3"): "proline_sidechain",
#     ("PRO", "CG"):  "proline_sidechain",
#     ("PRO", "HG2"): "proline_sidechain",
#     ("PRO", "HD2"): "proline_sidechain",
#     ("PRO", "HD3"): "proline_sidechain",
#     ("PRO", "CD"):  "proline_sidechain",
    
#     # --- Additional entries for Tryptophan (TRP) ---
#     ("TRP", "HZ3"): "tryptophan_sidechain",
#     ("TRP", "HH2"): "tryptophan_sidechain",
#     ("TRP", "CH2"): "tryptophan_sidechain",
#     ("TRP", "CZ3"): "tryptophan_sidechain",
    
    

#     # --- New entries for Alanine (ALA) ---
#     ("ALA", "H"):    "backbone",
#     ("ALA", "CA"):   "backbone",
#     ("ALA", "HA"):   "backbone",
#     ("ALA", "C"):    "backbone",
#     ("ALA", "CB"):   "alanine_sidechain",
#     ("ALA", "HB1"):  "alanine_sidechain",
#     ("ALA", "HB2"):  "alanine_sidechain",
#     ("ALA", "HB3"):  "alanine_sidechain",
    
#     # --- New entries for Glycine (GLY) ---
#     ("GLY", "H"):    "backbone",
#     ("GLY", "CA"):   "backbone",
#     ("GLY", "HA"):   "backbone",
#     ("GLY", "C"):    "backbone",
    
#     # --- New entries for Tyrosine (TYR) ---
#     ("TYR", "H"):    "backbone",
#     ("TYR", "CA"):   "backbone",
#     ("TYR", "HA"):   "backbone",
#     ("TYR", "C"):    "backbone",
#     ("TYR", "CB"):   "tyrosine_sidechain",
#     ("TYR", "HB2"):  "tyrosine_sidechain",
#     ("TYR", "HB3"):  "tyrosine_sidechain",
#     ("TYR", "CG"):   "tyrosine_sidechain",
#     ("TYR", "CD1"):  "tyrosine_sidechain",
#     ("TYR", "HD1"):  "tyrosine_sidechain",
#     ("TYR", "CE1"):  "tyrosine_sidechain",
#     ("TYR", "HE1"):  "tyrosine_sidechain",
#     ("TYR", "CZ"):   "tyrosine_sidechain",
#     ("TYR", "HH"):   "tyrosine_sidechain",
#     ("TYR", "CE2"):  "tyrosine_sidechain",
#     ("TYR", "HE2"):  "tyrosine_sidechain",
#     ("TYR", "CD2"):  "tyrosine_sidechain",
#     ("TYR", "HD2"):  "tyrosine_sidechain",
    
#     # --- New entries for Methionine (MET) ---
#     ("MET", "H"):    "backbone",
#     ("MET", "CA"):   "backbone",
#     ("MET", "HA"):   "backbone",
#     ("MET", "C"):    "backbone",
#     ("MET", "CB"):   "methionine_sidechain",
#     ("MET", "HB2"):  "methionine_sidechain",
#     ("MET", "HB3"):  "methionine_sidechain",
#     ("MET", "CG"):   "methionine_sidechain",
#     ("MET", "HG2"):  "methionine_sidechain",
#     ("MET", "HG3"):  "methionine_sidechain",
#     ("MET", "CE"):   "methionine_sidechain",
#     ("MET", "HE1"):  "methionine_sidechain",
#     ("MET", "HE2"):  "methionine_sidechain",
#     ("MET", "HE3"):  "methionine_sidechain",
    
    
    

#     # --- Additional entries for Isoleucine (ILE) ---
#     ("ILE", "H"):    "backbone",
#     ("ILE", "CA"):   "backbone",
#     ("ILE", "HA"):   "backbone",
#     ("ILE", "C"):    "backbone",
#     ("ILE", "CB"):   "isoleucine_sidechain",
#     ("ILE", "HB"):   "isoleucine_sidechain",
#     ("ILE", "CG2"):  "isoleucine_sidechain",
#     ("ILE", "HG21"): "isoleucine_sidechain",
#     ("ILE", "HG22"): "isoleucine_sidechain",
#     ("ILE", "HG23"): "isoleucine_sidechain",
    
#     # --- Additional entries for Leucine (LEU) ---
#     ("LEU", "CG"):   "leucine_sidechain",
#     ("LEU", "HG"):   "leucine_sidechain",
#     ("LEU", "CD1"):  "leucine_sidechain",
#     ("LEU", "HD11"): "leucine_sidechain",
#     ("LEU", "HD12"): "leucine_sidechain",
#     ("LEU", "HD13"): "leucine_sidechain",
#     ("LEU", "CD2"):  "leucine_sidechain",
#     ("LEU", "HD22"): "leucine_sidechain",
#     ("LEU", "HD23"): "leucine_sidechain",
#     # (Backbone atoms for LEU were added earlier.)
    
#     # --- Additional entries for Lysine (LYS) ---
#     ("LYS", "CA"):   "backbone",
#     ("LYS", "HA"):   "backbone",
#     ("LYS", "CB"):   "lysine_sidechain",
#     ("LYS", "HB3"):  "lysine_sidechain",
#     ("LYS", "CG"):   "lysine_sidechain",
#     ("LYS", "HG2"):  "lysine_sidechain",
#     ("LYS", "HG3"):  "lysine_sidechain",
#     ("LYS", "CD"):   "lysine_sidechain",
#     ("LYS", "HD2"):  "lysine_sidechain",
#     ("LYS", "HD3"):  "lysine_sidechain",
#     ("LYS", "CE"):   "lysine_sidechain",
#     ("LYS", "HE2"):  "lysine_sidechain",
#     ("LYS", "HE3"):  "lysine_sidechain",
#     ("LYS", "HZ1"):  "lysine_sidechain",
#     ("LYS", "HZ2"):  "lysine_sidechain",
#     ("LYS", "HZ3"):  "lysine_sidechain",
#     ("LYS", "C"):    "backbone",
#     # (LYS backbone H was added earlier.)
    
#     # --- Additional entries for Phenylalanine (PHE) ---
#     ("PHE", "H"):    "backbone",
#     ("PHE", "CA"):   "backbone",
#     ("PHE", "HA"):   "backbone",
#     ("PHE", "CB"):   "phenyl_sidechain",
#     ("PHE", "HB2"):  "phenyl_sidechain",
#     ("PHE", "HB3"):  "phenyl_sidechain",
#     ("PHE", "CG"):   "phenyl_sidechain",
#     ("PHE", "CE1"):  "phenyl_sidechain",
#     ("PHE", "HE1"):  "phenyl_sidechain",
#     ("PHE", "CD2"):  "phenyl_sidechain",
#     ("PHE", "HD2"):  "phenyl_sidechain",
#     ("PHE", "C"):    "backbone",
#     # (Other PHE entries such as HZ, CD1, CZ, CE2, HE2 were added earlier.)
    
#     # --- Additional entries for residue HD1 (if this residue exists) ---
#     ("HD1", "H"):    "backbone",
#     ("HD1", "CA"):   "backbone",
#     ("HD1", "HA"):   "backbone",
#     ("HD1", "C"):    "backbone",
#     ("HD1", "CB"):   "hd1_sidechain",
#     ("HD1", "CG"):   "hd1_sidechain",
#     ("HD1", "ND1"):  "hd1_sidechain",
#     ("HD1", "HD1"):  "hd1_sidechain",
#     ("HD1", "CE1"):  "hd1_sidechain",
#     ("HD1", "HE1"):  "hd1_sidechain",
#     ("HD1", "NE2"):  "hd1_sidechain",
#     ("HD1", "CD2"):  "hd1_sidechain",

        
#     # --- New entries for Phenylalanine (PHE) ---
#     ("PHE", "HZ"):  "phenyl_sidechain",   # Aromatic ring hydrogen
#     ("PHE", "CD1"): "phenyl_sidechain",
#     ("PHE", "HD1"): "phenyl_sidechain",
#     ("PHE", "CZ"):  "phenyl_sidechain",
#     ("PHE", "CE2"): "phenyl_sidechain",
#     ("PHE", "HE2"): "phenyl_sidechain",
    
#     # --- New entries for Isoleucine (ILE) ---
#     ("ILE", "HD12"): "isoleucine_sidechain",
#     ("ILE", "CD1"):  "isoleucine_sidechain",
#     ("ILE", "HD13"): "isoleucine_sidechain",
#     ("ILE", "HD11"): "isoleucine_sidechain",
#     ("ILE", "HG12"): "isoleucine_sidechain",
#     ("ILE", "CG1"):  "isoleucine_sidechain",
#     ("ILE", "HG13"): "isoleucine_sidechain",
    
#     # --- New entries for Leucine (LEU) ---
#     ("LEU", "HD21"): "leucine_sidechain",
#     ("LEU", "H"):    "backbone",
#     ("LEU", "CA"):   "backbone",
#     ("LEU", "HA"):   "backbone",
#     ("LEU", "CB"):   "leucine_sidechain",
#     ("LEU", "HB2"):  "leucine_sidechain",
#     ("LEU", "HB3"):  "leucine_sidechain",
#     ("LEU", "C"):    "backbone",
    
#     # --- New entries for Lysine (LYS) ---
#     ("LYS", "H"):    "backbone",
#     ("LYS", "HB2"):  "lysine_sidechain",
    
#     # --- New entries for residue HD1 (if this is a valid residue type) ---
#     ("HD1", "C"):    "backbone",
#     ("HD1", "CA"):   "backbone",
#     ("HD1", "HA"):   "backbone",
#     ("HD1", "HB2"):  "histidine_HD1_sidechain",
#     ("HD1", "HB3"):  "histidine_HD1_sidechain",
#     ("HD1", "HD2"):  "histidine_HD1_sidechain",
    
    


    
 
#     # --- Additional entries for Histidine (HIE) ---
#     ("HIE", "HB2"): "histidine_HIE_sidechain",
#     ("HIE", "HB3"): "histidine_HIE_sidechain",
#     ("HIE", "HE2"): "histidine_HIE_sidechain",
#     ("HIE", "HD2"): "histidine_HIE_sidechain",
    
#     # --- Additional entries for Valine (VAL) ---
#     ("VAL", "CA"):  "backbone",
#     ("VAL", "HA"):  "backbone",
#     ("VAL", "CB"):  "valine_sidechain",
#     ("VAL", "HB"):  "valine_sidechain",
#     ("VAL", "CG1"): "valine_sidechain",
#     ("VAL", "HG11"): "valine_sidechain",
#     ("VAL", "HG12"): "valine_sidechain",
#     ("VAL", "HG13"): "valine_sidechain",
#     ("VAL", "CG2"): "valine_sidechain",
#     ("VAL", "HG21"): "valine_sidechain",
#     ("VAL", "HG22"): "valine_sidechain",
#     ("VAL", "C"):   "backbone",

        
    
#     # --- New entries for Glutamate (GLU) ---
#     ("GLU", "H"):      "backbone",
#     ("GLU", "CA"):     "backbone",
#     ("GLU", "HA"):     "backbone",
#     ("GLU", "C"):      "backbone",
#     ("GLU", "CB"):     "glutamate_sidechain",
#     ("GLU", "HB2"):    "glutamate_sidechain",
#     ("GLU", "HB3"):    "glutamate_sidechain",
#     ("GLU", "CG"):     "glutamate_sidechain",
#     ("GLU", "HG2"):    "glutamate_sidechain",
#     ("GLU", "HG3"):    "glutamate_sidechain",
    
#     # --- New entries for Histidine in the HIE form (epsilon-protonated) ---
#     ("HIE", "H"):      "backbone",
#     ("HIE", "CA"):     "backbone",
#     ("HIE", "HA"):     "backbone",
#     ("HIE", "C"):      "backbone",
#     ("HIE", "ND1"):    "imidazole",
#     ("HIE", "NE2"):    "imidazole",
#     ("HIE", "CG"):     "imidazole",
#     ("HIE", "CD2"):    "imidazole",
#     ("HIE", "CE1"):    "imidazole",
#     ("HIE", "HE1"):    "imidazole",
#     ("HIE", "CB"):     "histidine_sidechain",
    
#     # --- New entries for Valine (VAL) ---
#     ("VAL", "H"):      "backbone",
#     ("VAL", "HG23"):   "valine_sidechain",

        
#     # --- New entries for Threonine (THR) ---
#     ("THR", "H"):      "backbone",
#     ("THR", "CA"):     "backbone",
#     ("THR", "HA"):     "backbone",
#     ("THR", "C"):      "backbone",
#     ("THR", "CB"):     "threonine_sidechain",
#     ("THR", "HB"):     "threonine_sidechain",
#     ("THR", "CG2"):    "threonine_sidechain",
#     ("THR", "HG21"):   "threonine_sidechain",
#     ("THR", "HG22"):   "threonine_sidechain",
#     ("THR", "HG23"):   "threonine_sidechain",
#     ("THR", "HG1"):    "threonine_sidechain",
    
#     # --- New entries for Serine (SER) ---
#     ("SER", "H"):      "backbone",
#     ("SER", "CA"):     "backbone",
#     ("SER", "HA"):     "backbone",
#     ("SER", "C"):      "backbone",
#     ("SER", "CB"):     "serine_sidechain",
#     ("SER", "HB2"):    "serine_sidechain",
#     ("SER", "HB3"):    "serine_sidechain",
#     ("SER", "HG"):     "serine_sidechain",
    
    
    

#     # COO moiety (Asp, Glu): includes the two oxygens and the bridging carbon
#     ("ASP", "OD1"): "COO",
#     ("ASP", "OD2"): "COO",
#     ("ASP", "CG"): "COO",  # Bridging carbon
#     ("GLU", "OE1"): "COO",
#     ("GLU", "OE2"): "COO",
#     ("GLU", "CD"): "COO",  # Bridging carbon
    
#     # Imidazole moiety (His): includes all 5 atoms in the imidazole ring
#     ("HIS", "ND1"): "imidazole",
#     ("HIS", "NE2"): "imidazole",
#     ("HIS", "CG"): "imidazole",
#     ("HIS", "CD2"): "imidazole",
#     ("HIS", "CE1"): "imidazole",
    
#     # C-terminal moiety (C_terminus): similar to COO, includes two oxygens and the bridging carbon
#     ("ANY", "OXT"): "c_terminus",  # Oxygen 1
#     #("ANY", "C"): "c_terminus",    # Bridging carbon
    
#     # Serine (Ser)
#     ("SER", "OG"): "hydroxyl",

#     # Threonine (Thr)
#     ("THR", "OG1"): "hydroxyl",

#     # Tyrosine (Tyr)
#     ("TYR", "OH"): "phenol",

#     # Cysteine (Cys)
#     ("CYS", "SG"): "thiol",

#     # Lysine (Lys)
#     ("LYS", "NZ"): "amine",

#     # Arginine (Arg)
#     ("ARG", "NE"): "guanidinium",
#     ("ARG", "NH1"): "guanidinium",
#     ("ARG", "NH2"): "guanidinium",

#     # Asparagine (Asn)
#     ("ASN", "OD1"): "amide",
#     ("ASN", "ND2"): "amide",

#     # Glutamine (Gln)
#     ("GLN", "OE1"): "amide",
#     ("GLN", "NE2"): "amide",

#     # Methionine (Met)
#     ("MET", "SD"): "thioether",

#     # Backbone atoms (common to all residues)
#     ("ANY", "O"): "backbone",   # Double-bonded oxygen in peptide bond
#     ("ANY", "N"): "backbone",
    
#     # Termini (common to all)
#     ("ANY", "OXT"): "c_terminus",  # C-terminal carboxyl oxygen)

#     # Water molecule
#     ("HOH", "O"): "water",
    
#     ("WAT", "H1"): "water",
#     #("WAT", "H2"): "water",
#     ("WAT", "H2"): "water",





#     # O2 Oxy ligand
#     ("OY1", "O1"): "oxy",  # C-terminal carboxyl oxygen)
#     ("OY1", "O2"): "oxy",
    
    
            
#     # --- New entries for Arginine (ARG) ---
#     # Backbone atoms
#     ("ARG", "H"):  "backbone",
#     ("ARG", "CA"): "backbone",
#     ("ARG", "HA"): "backbone",
#     ("ARG", "C"):  "backbone",
#     # Side-chain atoms
#     ("ARG", "CB"):    "arginine_sidechain",
#     ("ARG", "HB2"):   "arginine_sidechain",
#     ("ARG", "HB3"):   "arginine_sidechain",
#     ("ARG", "CG"):    "arginine_sidechain",
#     ("ARG", "HG2"):   "arginine_sidechain",
#     ("ARG", "HG3"):   "arginine_sidechain",
#     ("ARG", "CD"):    "arginine_sidechain",
#     ("ARG", "HD2"):   "arginine_sidechain",
#     ("ARG", "HD3"):   "arginine_sidechain",
#     ("ARG", "HE"):    "arginine_sidechain",
#     ("ARG", "CZ"):    "arginine_sidechain",
#     ("ARG", "HH11"):  "arginine_sidechain",
#     ("ARG", "HH12"):  "arginine_sidechain",
#     ("ARG", "HH21"):  "arginine_sidechain",
#     ("ARG", "HH22"):  "arginine_sidechain",
    
#     # --- New entries for Aspartate (ASP) ---
#     # Backbone atoms
#     ("ASP", "H"):  "backbone",
#     ("ASP", "CA"): "backbone",
#     ("ASP", "HA"): "backbone",
#     ("ASP", "C"):  "backbone",
#     # Side-chain atoms (aside from the carboxylate already defined)
#     ("ASP", "CB"):   "aspartate_sidechain",
#     ("ASP", "HB2"):  "aspartate_sidechain",
#     ("ASP", "HB3"):  "aspartate_sidechain",



    
#     # For HM1: All other atoms in the heme porphyrin macrocycle
#     ("HM1", "CHA"): "heme_por",
#     ("HM1", "CHB"): "heme_por",
#     ("HM1", "CHC"): "heme_por",
#     ("HM1", "CHD"): "heme_por",
#     ("HM1", "NA"): "heme_por",
#     ("HM1", "C1A"): "heme_por",
#     ("HM1", "C2A"): "heme_por",
#     ("HM1", "C3A"): "heme_por",
#     ("HM1", "C4A"): "heme_por",
#     ("HM1", "CMA"): "heme_por",
#     ("HM1", "CAA"): "heme_por",
#     ("HM1", "CBA"): "heme_por",
#     ("HM1", "CGA"): "heme_por",
#     ("HM1", "O1A"): "heme_por",
#     ("HM1", "O2A"): "heme_por",
#     ("HM1", "NB"): "heme_por",
#     ("HM1", "C1B"): "heme_por",
#     ("HM1", "C2B"): "heme_por",
#     ("HM1", "C3B"): "heme_por",
#     ("HM1", "C4B"): "heme_por",
#     ("HM1", "CMB"): "heme_por",
#     ("HM1", "CAB"): "heme_por",
#     ("HM1", "CBB"): "heme_por",
#     ("HM1", "NC"): "heme_por",
#     ("HM1", "C1C"): "heme_por",
#     ("HM1", "C2C"): "heme_por",
#     ("HM1", "C3C"): "heme_por",
#     ("HM1", "C4C"): "heme_por",
#     ("HM1", "CMC"): "heme_por",
#     ("HM1", "CAC"): "heme_por",
#     ("HM1", "CBC"): "heme_por",
#     ("HM1", "ND"): "heme_por",
#     ("HM1", "C1D"): "heme_por",
#     ("HM1", "C2D"): "heme_por",
#     ("HM1", "C3D"): "heme_por",
#     ("HM1", "C4D"): "heme_por",
#     ("HM1", "CMD"): "heme_por",
#     ("HM1", "CAD"): "heme_por",
#     ("HM1", "CBD"): "heme_por",
#     ("HM1", "CGD"): "heme_por",
#     ("HM1", "O1D"): "heme_por",
#     ("HM1", "O2D"): "heme_por",
#     ("HM1", "HMA1"): "heme_por",
#     ("HM1", "HMA2"): "heme_por",
#     ("HM1", "HMA3"): "heme_por",
#     ("HM1", "HMB1"): "heme_por",
#     ("HM1", "HMB2"): "heme_por",
#     ("HM1", "HMB3"): "heme_por",
#     ("HM1", "HMC1"): "heme_por",
#     ("HM1", "HMC2"): "heme_por",
#     ("HM1", "HMC3"): "heme_por",
#     ("HM1", "HMD1"): "heme_por",
#     ("HM1", "HMD2"): "heme_por",
#     ("HM1", "HMD3"): "heme_por",
#     ("HM1", "HBB1"): "heme_por",
#     ("HM1", "HBB2"): "heme_por",
#     ("HM1", "HBC1"): "heme_por",
#     ("HM1", "HBC2"): "heme_por",
#     ("HM1", "HBA1"): "heme_por",
#     ("HM1", "HBA2"): "heme_por",
#     ("HM1", "HAA1"): "heme_por",
#     ("HM1", "HAA2"): "heme_por",
#     ("HM1", "HBD1"): "heme_por",
#     ("HM1", "HBD2"): "heme_por",
#     ("HM1", "HAD1"): "heme_por",
#     ("HM1", "HAD2"): "heme_por",
#     ("HM1", "HHA"): "heme_por",
#     ("HM1", "HHB"): "heme_por",
#     ("HM1", "HHC"): "heme_por",
#     ("HM1", "HHD"): "heme_por",
#     ("HM1", "HAB"): "heme_por",
#     ("HM1", "HAC"): "heme_por",

#     # For HM2: FE atoms (the central iron)
#     ("HM2", "FE"): "heme_fe",


    
#     # For HEM: All other atoms in the heme porphyrin macrocycle
#     ("HEM", "CHA"): "heme_por",
#     ("HEM", "CHB"): "heme_por",
#     ("HEM", "CHC"): "heme_por",
#     ("HEM", "CHD"): "heme_por",
#     ("HEM", "NA"): "heme_por",
#     ("HEM", "C1A"): "heme_por",
#     ("HEM", "C2A"): "heme_por",
#     ("HEM", "C3A"): "heme_por",
#     ("HEM", "C4A"): "heme_por",
#     ("HEM", "CMA"): "heme_por",
#     ("HEM", "CAA"): "heme_por",
#     ("HEM", "CBA"): "heme_por",
#     ("HEM", "CGA"): "heme_por",
#     ("HEM", "O1A"): "heme_por",
#     ("HEM", "O2A"): "heme_por",
#     ("HEM", "NB"): "heme_por",
#     ("HEM", "C1B"): "heme_por",
#     ("HEM", "C2B"): "heme_por",
#     ("HEM", "C3B"): "heme_por",
#     ("HEM", "C4B"): "heme_por",
#     ("HEM", "CMB"): "heme_por",
#     ("HEM", "CAB"): "heme_por",
#     ("HEM", "CBB"): "heme_por",
#     ("HEM", "NC"): "heme_por",
#     ("HEM", "C1C"): "heme_por",
#     ("HEM", "C2C"): "heme_por",
#     ("HEM", "C3C"): "heme_por",
#     ("HEM", "C4C"): "heme_por",
#     ("HEM", "CMC"): "heme_por",
#     ("HEM", "CAC"): "heme_por",
#     ("HEM", "CBC"): "heme_por",
#     ("HEM", "ND"): "heme_por",
#     ("HEM", "C1D"): "heme_por",
#     ("HEM", "C2D"): "heme_por",
#     ("HEM", "C3D"): "heme_por",
#     ("HEM", "C4D"): "heme_por",
#     ("HEM", "CMD"): "heme_por",
#     ("HEM", "CAD"): "heme_por",
#     ("HEM", "CBD"): "heme_por",
#     ("HEM", "CGD"): "heme_por",
#     ("HEM", "O1D"): "heme_por",
#     ("HEM", "O2D"): "heme_por",
#     ("HEM", "HMA1"): "heme_por",
#     ("HEM", "HMA2"): "heme_por",
#     ("HEM", "HMA3"): "heme_por",
#     ("HEM", "HMB1"): "heme_por",
#     ("HEM", "HMB2"): "heme_por",
#     ("HEM", "HMB3"): "heme_por",
#     ("HEM", "HMC1"): "heme_por",
#     ("HEM", "HMC2"): "heme_por",
#     ("HEM", "HMC3"): "heme_por",
#     ("HEM", "HMD1"): "heme_por",
#     ("HEM", "HMD2"): "heme_por",
#     ("HEM", "HMD3"): "heme_por",
#     ("HEM", "HBB1"): "heme_por",
#     ("HEM", "HBB2"): "heme_por",
#     ("HEM", "HBC1"): "heme_por",
#     ("HEM", "HBC2"): "heme_por",
#     ("HEM", "HBA1"): "heme_por",
#     ("HEM", "HBA2"): "heme_por",
#     ("HEM", "HAA1"): "heme_por",
#     ("HEM", "HAA2"): "heme_por",
#     ("HEM", "HBD1"): "heme_por",
#     ("HEM", "HBD2"): "heme_por",
#     ("HEM", "HAD1"): "heme_por",
#     ("HEM", "HAD2"): "heme_por",
#     ("HEM", "HHA"): "heme_por",
#     ("HEM", "HHB"): "heme_por",
#     ("HEM", "HHC"): "heme_por",
#     ("HEM", "HHD"): "heme_por",
#     ("HEM", "HAB"): "heme_por",
#     ("HEM", "HAC"): "heme_por",

#     ("HEM", "FE"): "heme_fe",
    
#     ("SO2", "S"): "sulfur_dioxide",
#     ("SO2", "O2"): "sulfur_dioxide",
#     ("SO2", "O1"): "sulfur_dioxide",


    
    
#     # FeMo-cofactor (ICS)
#     ("ICS", "MO"): "femoco_mo",        # Molybdenum
#     ("ICS", "FE1"): "femoco_fe",
#     ("ICS", "FE2"): "femoco_fe",
#     ("ICS", "FE3"): "femoco_fe",
#     ("ICS", "FE4"): "femoco_fe",
#     ("ICS", "FE5"): "femoco_fe",
#     ("ICS", "FE6"): "femoco_fe",
#     ("ICS", "FE7"): "femoco_fe",
#     ("ICS", "S1A"): "femoco_s",
#     ("ICS", "S2B"): "femoco_s",
#     ("ICS", "S3A"): "femoco_s",
#     ("ICS", "S3B"): "femoco_s",
#     ("ICS", "S4A"): "femoco_s",
#     ("ICS", "S4B"): "femoco_s",
#     ("ICS", "S5A"): "femoco_s",
#     ("ICS", "S5B"): "femoco_s",
#     ("ICS", "CX"): "femoco_carbide",   # central carbide
#     ("ICS", "O1"): "femoco_o",         # homocitrate coordination oxygens (sometimes labeled here)
#     ("ICS", "O2"): "femoco_o",
#     ("ICS", "O3"): "femoco_o",
#     ("ICS", "N1"): "femoco_n",         # belt nitrogens (occasionally appear in ICS block)
#     ("ICS", "N2"): "femoco_n",
    
    
        
        
#     # Homocitrate (HCA)
#     ("HCA", "C1"): "hca_c",
#     ("HCA", "C2"): "hca_c",
#     ("HCA", "C3"): "hca_c",
#     ("HCA", "C4"): "hca_c",
#     ("HCA", "C5"): "hca_c",
#     ("HCA", "C6"): "hca_c",
#     ("HCA", "O1"): "hca_o",        # α-carboxyl oxygen
#     ("HCA", "O2"): "hca_o",
#     ("HCA", "O3"): "hca_o",        # hydroxyl / terminal carboxyl oxygens
#     ("HCA", "O4"): "hca_o",
#     ("HCA", "O5"): "hca_o",
#     ("HCA", "O6"): "hca_o",
#     ("HCA", "H1"): "hca_h",
#     ("HCA", "H2"): "hca_h",
#     ("HCA", "H3"): "hca_h",
#     ("HCA", "H4"): "hca_h",
#     ("HCA", "H5"): "hca_h",
#     ("HCA", "H6"): "hca_h",
    
    
#     # --- New entries for Oxygen-Evolving Cluster (OEX) ---
#     ("OEX", "O1"):   "oxygen_evolving_cluster",
#     ("OEX", "CA1"):  "oxygen_evolving_cluster",
#     ("OEX", "MN1"):  "oxygen_evolving_cluster",
#     ("OEX", "O2"):   "oxygen_evolving_cluster",
#     ("OEX", "MN2"):  "oxygen_evolving_cluster",
#     ("OEX", "O3"):   "oxygen_evolving_cluster",
#     ("OEX", "MN3"):  "oxygen_evolving_cluster",
#     ("OEX", "O4"):   "oxygen_evolving_cluster",
#     ("OEX", "MN4"):  "oxygen_evolving_cluster",
#     ("OEX", "O5"):   "oxygen_evolving_cluster",


    

# }
    





chemical_moieties: Dict[tuple, str] = {




    # EXTRAS to add-in:
    ("CYS", "CB"): "cysteine_sidechain",
    ("PRO", "CB"): "proline_sidechain",


    # =========================
    # MISSING BACKBONES (filled)
    # =========================
    # Universal backbone set for standard residues (and common HIS variants)
    ("ALA","N"): "backbone", ("ALA","H"): "backbone", ("ALA","CA"): "backbone", ("ALA","HA"): "backbone", ("ALA","C"): "backbone", ("ALA","O"): "backbone", ("ALA","OXT"): "backbone",
    ("ARG","N"): "backbone", ("ARG","H"): "backbone", ("ARG","CA"): "backbone", ("ARG","HA"): "backbone", ("ARG","C"): "backbone", ("ARG","O"): "backbone", ("ARG","OXT"): "backbone",
    ("ASN","N"): "backbone", ("ASN","H"): "backbone", ("ASN","CA"): "backbone", ("ASN","HA"): "backbone", ("ASN","C"): "backbone", ("ASN","O"): "backbone", ("ASN","OXT"): "backbone",
    ("ASP","N"): "backbone", ("ASP","H"): "backbone", ("ASP","CA"): "backbone", ("ASP","HA"): "backbone", ("ASP","C"): "backbone", ("ASP","O"): "backbone", ("ASP","OXT"): "backbone",
    ("CYS","N"): "backbone", ("CYS","H"): "backbone", ("CYS","CA"): "backbone", ("CYS","HA"): "backbone", ("CYS","C"): "backbone", ("CYS","O"): "backbone", ("CYS","OXT"): "backbone",
    ("GLN","N"): "backbone", ("GLN","H"): "backbone", ("GLN","CA"): "backbone", ("GLN","HA"): "backbone", ("GLN","C"): "backbone", ("GLN","O"): "backbone", ("GLN","OXT"): "backbone",
    ("GLU","N"): "backbone", ("GLU","H"): "backbone", ("GLU","CA"): "backbone", ("GLU","HA"): "backbone", ("GLU","C"): "backbone", ("GLU","O"): "backbone", ("GLU","OXT"): "backbone",
    ("GLY","N"): "backbone", ("GLY","H"): "backbone", ("GLY","CA"): "backbone", ("GLY","HA"): "backbone", ("GLY","HA2"): "backbone", ("GLY","HA3"): "backbone", ("GLY","C"): "backbone", ("GLY","O"): "backbone", ("GLY","OXT"): "backbone",
    ("HIS","N"): "backbone", ("HIS","H"): "backbone", ("HIS","CA"): "backbone", ("HIS","HA"): "backbone", ("HIS","C"): "backbone", ("HIS","O"): "backbone", ("HIS","OXT"): "backbone",
    ("HID","N"): "backbone", ("HID","H"): "backbone", ("HID","CA"): "backbone", ("HID","HA"): "backbone", ("HID","C"): "backbone", ("HID","O"): "backbone", ("HID","OXT"): "backbone",
    ("HIE","N"): "backbone", ("HIE","H"): "backbone", ("HIE","CA"): "backbone", ("HIE","HA"): "backbone", ("HIE","C"): "backbone", ("HIE","O"): "backbone", ("HIE","OXT"): "backbone",
    ("HIP","N"): "backbone", ("HIP","H"): "backbone", ("HIP","CA"): "backbone", ("HIP","HA"): "backbone", ("HIP","C"): "backbone", ("HIP","O"): "backbone", ("HIP","OXT"): "backbone",
    ("ILE","N"): "backbone", ("ILE","H"): "backbone", ("ILE","CA"): "backbone", ("ILE","HA"): "backbone", ("ILE","C"): "backbone", ("ILE","O"): "backbone", ("ILE","OXT"): "backbone",
    ("LEU","N"): "backbone", ("LEU","H"): "backbone", ("LEU","CA"): "backbone", ("LEU","HA"): "backbone", ("LEU","C"): "backbone", ("LEU","O"): "backbone", ("LEU","OXT"): "backbone",
    ("LYS","N"): "backbone", ("LYS","H"): "backbone", ("LYS","CA"): "backbone", ("LYS","HA"): "backbone", ("LYS","C"): "backbone", ("LYS","O"): "backbone", ("LYS","OXT"): "backbone",
    ("MET","N"): "backbone", ("MET","H"): "backbone", ("MET","CA"): "backbone", ("MET","HA"): "backbone", ("MET","C"): "backbone", ("MET","O"): "backbone", ("MET","OXT"): "backbone",
    ("PHE","N"): "backbone", ("PHE","H"): "backbone", ("PHE","CA"): "backbone", ("PHE","HA"): "backbone", ("PHE","C"): "backbone", ("PHE","O"): "backbone", ("PHE","OXT"): "backbone",
    ("PRO","N"): "backbone", ("PRO","CA"): "backbone", ("PRO","HA"): "backbone", ("PRO","C"): "backbone", ("PRO","O"): "backbone", ("PRO","OXT"): "backbone",
    ("SER","N"): "backbone", ("SER","H"): "backbone", ("SER","CA"): "backbone", ("SER","HA"): "backbone", ("SER","C"): "backbone", ("SER","O"): "backbone", ("SER","OXT"): "backbone",
    ("THR","N"): "backbone", ("THR","H"): "backbone", ("THR","CA"): "backbone", ("THR","HA"): "backbone", ("THR","C"): "backbone", ("THR","O"): "backbone", ("THR","OXT"): "backbone",
    ("TRP","N"): "backbone", ("TRP","H"): "backbone", ("TRP","CA"): "backbone", ("TRP","HA"): "backbone", ("TRP","C"): "backbone", ("TRP","O"): "backbone", ("TRP","OXT"): "backbone",
    ("TYR","N"): "backbone", ("TYR","H"): "backbone", ("TYR","CA"): "backbone", ("TYR","HA"): "backbone", ("TYR","C"): "backbone", ("TYR","O"): "backbone", ("TYR","OXT"): "backbone",
    ("VAL","N"): "backbone", ("VAL","H"): "backbone", ("VAL","CA"): "backbone", ("VAL","HA"): "backbone", ("VAL","C"): "backbone", ("VAL","O"): "backbone", ("VAL","OXT"): "backbone",

    # =====================
    # OXYGEN (diatomic O2)
    # =====================
    ("OXY", "O1"): "oxy_ligand",
    ("OXY", "O2"): "oxy_ligand",

    # =========
    # HEME: HEA
    # =========
    ("HEA", "C16"): "heme",
    ("HEA", "C17"): "heme",
    ("HEA,","C18"): "heme",  # (kept as provided; if this comma is a typo in source, change to ("HEA","C18"))
    ("HEA", "C18"): "heme",
    ("HEA", "C19"): "heme",
    ("HEA", "C20"): "heme",
    ("HEA", "C21"): "heme",
    ("HEA", "C22"): "heme",
    ("HEA", "C23"): "heme",
    ("HEA", "C24"): "heme",
    ("HEA", "C25"): "heme",
    ("HEA", "C27"): "heme",

    # ("HIS", "CA"): "backbone",
    ("HIS", "CB"): "histidine_sidechain",

    ("HEA", "C13"): "heme",
    ("HEA", "C14"): "heme",
    ("HEA", "C15"): "heme",
    ("HEA", "C26"): "heme",
    ("HEA", "O1"): "heme",
    ("HEA", "O2"): "heme",
    ("HEA", "CAA"): "heme",
    ("HEA", "CAC"): "heme",
    ("HEA", "CAD"): "heme",
    ("HEA", "NA"): "heme",
    ("HEA", "CBA"): "heme",
    ("HEA", "CBC"): "heme",
    ("HEA", "CBD"): "heme",
    ("HEA", "NB"): "heme",
    ("HEA", "CGA"): "heme",
    ("HEA", "CGD"): "heme",
    ("HEA", "ND"): "heme",
    ("HEA", "CHA"): "heme",
    ("HEA", "CHB"): "heme",
    ("HEA", "CHC"): "heme",
    ("HEA", "CHD"): "heme",
    ("HEA", "CMA"): "heme",
    ("HEA", "CMB"): "heme",
    ("HEA", "CMC"): "heme",
    ("HEA", "CMD"): "heme",
    ("HEA", "OMA"): "heme",
    ("HEA", "C1A"): "heme",
    ("HEA", "C1B"): "heme",
    ("HEA", "C1C"): "heme",
    ("HEA", "C1D"): "heme",
    ("HEA", "O1A"): "heme",
    ("HEA", "O1D"): "heme",
    ("HEA", "C2A"): "heme",
    ("HEA", "C2B"): "heme",
    ("HEA", "C2C"): "heme",
    ("HEA", "C2D"): "heme",
    ("HEA", "O2A"): "heme",
    ("HEA", "O2D"): "heme",
    ("HEA", "C3A"): "heme",
    ("HEA", "C3B"): "heme",
    ("HEA", "C3C"): "heme",
    ("HEA", "C3D"): "heme",
    ("HEA", "C4A"): "heme",
    ("HEA", "C4B"): "heme",
    ("HEA", "C4C"): "heme",
    ("HEA", "C4D"): "heme",
    ("HEA", "C11"): "heme",
    ("HEA", "O11"): "heme",
    ("HEA", "C12"): "heme",
    ("HEA", "NC"): "heme",
    ("HEA", "FE"): "heme",

    # ==========
    # CU centers
    # ==========
    ("CU",  "CU"): "copper_center",
    ("CU1", "CU"): "copper_center",

    # ======================
    # TRP (sidechain extras)
    # ======================
    ("TRP", "CZ3"): "tryptophan_sidechain",

    # ======================
    # GLN (newer entries)
    # ======================
    ("GLN", "CB"):   "glutamine_sidechain",
    ("GLN", "HB2"):  "glutamine_sidechain",
    ("GLN", "HB3"):  "glutamine_sidechain",
    ("GLN", "CG"):   "glutamine_sidechain",
    ("GLN", "HG2"):  "glutamine_sidechain",
    ("GLN", "HG3"):  "glutamine_sidechain",
    ("GLN", "CD"):   "glutamine_sidechain",
    ("GLN", "HE21"): "glutamine_sidechain",
    ("GLN", "HE22"): "glutamine_sidechain",

    # =======================
    # PRO (missing sidechain)
    # =======================
    ("PRO", "HG3"): "proline_sidechain",

    # ======================
    # TRP (full sidechain)
    # ======================
    ("TRP", "CB"):   "tryptophan_sidechain",
    ("TRP", "HB2"):  "tryptophan_sidechain",
    ("TRP", "HB3"):  "tryptophan_sidechain",
    ("TRP", "CG"):   "tryptophan_sidechain",
    ("TRP", "CD1"):  "tryptophan_sidechain",
    ("TRP", "HD1"):  "tryptophan_sidechain",
    ("TRP", "NE1"):  "tryptophan_sidechain",
    ("TRP", "HE1"):  "tryptophan_sidechain",
    ("TRP", "CE2"):  "tryptophan_sidechain",
    ("TRP", "CZ2"):  "tryptophan_sidechain",
    ("TRP", "HZ2"):  "tryptophan_sidechain",
    ("TRP", "CE3"):  "tryptophan_sidechain",
    ("TRP", "HE3"):  "tryptophan_sidechain",
    ("TRP", "CD2"):  "tryptophan_sidechain",
    ("TRP", "HZ3"):  "tryptophan_sidechain",
    ("TRP", "HH2"):  "tryptophan_sidechain",
    ("TRP", "CH2"):  "tryptophan_sidechain",
    ("TRP", "CZ3"):  "tryptophan_sidechain",

    # ======================
    # ASN (additions)
    # ======================
    ("ASN", "CB"):    "asparagine_sidechain",
    ("ASN", "HB2"):   "asparagine_sidechain",
    ("ASN", "HB3"):   "asparagine_sidechain",
    ("ASN", "CG"):    "asparagine_sidechain",
    ("ASN", "HD21"):  "asparagine_sidechain",
    ("ASN", "HD22"):  "asparagine_sidechain",

    # =========
    # ALA (new)
    # =========
    ("ALA", "CB"):   "alanine_sidechain",
    ("ALA", "HB1"):  "alanine_sidechain",
    ("ALA", "HB2"):  "alanine_sidechain",
    ("ALA", "HB3"):  "alanine_sidechain",

    # =========
    # TYR (new)
    # =========
    ("TYR", "CB"):   "tyrosine_sidechain",
    ("TYR", "HB2"):  "tyrosine_sidechain",
    ("TYR", "HB3"):  "tyrosine_sidechain",
    ("TYR", "CG"):   "tyrosine_sidechain",
    ("TYR", "CD1"):  "tyrosine_sidechain",
    ("TYR", "HD1"):  "tyrosine_sidechain",
    ("TYR", "CE1"):  "tyrosine_sidechain",
    ("TYR", "HE1"):  "tyrosine_sidechain",
    ("TYR", "CZ"):   "tyrosine_sidechain",
    ("TYR", "HH"):   "tyrosine_sidechain",
    ("TYR", "CE2"):  "tyrosine_sidechain",
    ("TYR", "HE2"):  "tyrosine_sidechain",
    ("TYR", "CD2"):  "tyrosine_sidechain",
    ("TYR", "HD2"):  "tyrosine_sidechain",

    # =========
    # MET (new)
    # =========
    ("MET", "CB"):   "methionine_sidechain",
    ("MET", "HB2"):  "methionine_sidechain",
    ("MET", "HB3"):  "methionine_sidechain",
    ("MET", "CG"):   "methionine_sidechain",
    ("MET", "HG2"):  "methionine_sidechain",
    ("MET", "HG3"):  "methionine_sidechain",
    ("MET", "SD"):   "thioether",
    ("MET", "CE"):   "methionine_sidechain",
    ("MET", "HE1"):  "methionine_sidechain",
    ("MET", "HE2"):  "methionine_sidechain",
    ("MET", "HE3"):  "methionine_sidechain",

    # =========
    # ILE (new)
    # =========
    ("ILE", "CB"):   "isoleucine_sidechain",
    ("ILE", "HB"):   "isoleucine_sidechain",
    ("ILE", "CG2"):  "isoleucine_sidechain",
    ("ILE", "HG21"): "isoleucine_sidechain",
    ("ILE", "HG22"): "isoleucine_sidechain",
    ("ILE", "HG23"): "isoleucine_sidechain",
    ("ILE", "HD12"): "isoleucine_sidechain",
    ("ILE", "CD1"):  "isoleucine_sidechain",
    ("ILE", "HD13"): "isoleucine_sidechain",
    ("ILE", "HD11"): "isoleucine_sidechain",
    ("ILE", "HG12"): "isoleucine_sidechain",
    ("ILE", "CG1"):  "isoleucine_sidechain",
    ("ILE", "HG13"): "isoleucine_sidechain",

    # =========
    # LEU (new)
    # =========
    ("LEU", "CB"):   "leucine_sidechain",
    ("LEU", "HB2"):  "leucine_sidechain",
    ("LEU", "HB3"):  "leucine_sidechain",
    ("LEU", "CG"):   "leucine_sidechain",
    ("LEU", "HG"):   "leucine_sidechain",
    ("LEU", "CD1"):  "leucine_sidechain",
    ("LEU", "HD11"): "leucine_sidechain",
    ("LEU", "HD12"): "leucine_sidechain",
    ("LEU", "HD13"): "leucine_sidechain",
    ("LEU", "CD2"):  "leucine_sidechain",
    ("LEU", "HD21"): "leucine_sidechain",
    ("LEU", "HD22"): "leucine_sidechain",
    ("LEU", "HD23"): "leucine_sidechain",

    # =========
    # LYS (new)
    # =========
    ("LYS", "CB"):   "lysine_sidechain",
    ("LYS", "HB2"):  "lysine_sidechain",
    ("LYS", "HB3"):  "lysine_sidechain",
    ("LYS", "CG"):   "lysine_sidechain",
    ("LYS", "HG2"):  "lysine_sidechain",
    ("LYS", "HG3"):  "lysine_sidechain",
    ("LYS", "CD"):   "lysine_sidechain",
    ("LYS", "HD2"):  "lysine_sidechain",
    ("LYS", "HD3"):  "lysine_sidechain",
    ("LYS", "CE"):   "lysine_sidechain",
    ("LYS", "HE2"):  "lysine_sidechain",
    ("LYS", "HE3"):  "lysine_sidechain",
    ("LYS", "NZ"):   "amine",
    ("LYS", "HZ1"):  "lysine_sidechain",
    ("LYS", "HZ2"):  "lysine_sidechain",
    ("LYS", "HZ3"):  "lysine_sidechain",

    # =========
    # PHE (new)
    # =========
    ("PHE", "CB"):  "phenyl_sidechain",
    ("PHE", "HB2"): "phenyl_sidechain",
    ("PHE", "HB3"): "phenyl_sidechain",
    ("PHE", "CG"):  "phenyl_sidechain",
    ("PHE", "CD1"): "phenyl_sidechain",
    ("PHE", "HD1"): "phenyl_sidechain",
    ("PHE", "CD2"): "phenyl_sidechain",
    ("PHE", "HD2"): "phenyl_sidechain",
    ("PHE", "CE1"): "phenyl_sidechain",
    ("PHE", "HE1"): "phenyl_sidechain",
    ("PHE", "CE2"): "phenyl_sidechain",
    ("PHE", "HE2"): "phenyl_sidechain",
    ("PHE", "CZ"):  "phenyl_sidechain",
    ("PHE", "HZ"):  "phenyl_sidechain",

    # =========
    # HIE (extra)
    # =========
    ("HIE", "ND1"): "imidazole",
    ("HIE", "NE2"): "imidazole",
    ("HIE", "CG"):  "imidazole",
    ("HIE", "CD2"): "imidazole",
    ("HIE", "CE1"): "imidazole",
    ("HIE", "HE1"): "imidazole",
    ("HIE", "CB"):  "histidine_sidechain",
    ("HIE", "HB2"): "histidine_HIE_sidechain",
    ("HIE", "HB3"): "histidine_HIE_sidechain",
    ("HIE", "HE2"): "histidine_HIE_sidechain",
    ("HIE", "HD2"): "histidine_HIE_sidechain",

    # =========
    # VAL (new)
    # =========
    ("VAL", "CB"):   "valine_sidechain",
    ("VAL", "HB"):   "valine_sidechain",
    ("VAL", "CG1"):  "valine_sidechain",
    ("VAL", "HG11"): "valine_sidechain",
    ("VAL", "HG12"): "valine_sidechain",
    ("VAL", "HG13"): "valine_sidechain",
    ("VAL", "CG2"):  "valine_sidechain",
    ("VAL", "HG21"): "valine_sidechain",
    ("VAL", "HG22"): "valine_sidechain",
    ("VAL", "HG23"): "valine_sidechain",

    # =========
    # THR (new)
    # =========
    ("THR", "CB"):   "threonine_sidechain",
    ("THR", "HB"):   "threonine_sidechain",
    ("THR", "CG2"):  "threonine_sidechain",
    ("THR", "HG21"): "threonine_sidechain",
    ("THR", "HG22"): "threonine_sidechain",
    ("THR", "HG23"): "threonine_sidechain",
    ("THR", "OG1"):  "hydroxyl",
    ("THR", "HG1"):  "threonine_sidechain",

    # =========
    # SER (new)
    # =========
    ("SER", "CB"):  "serine_sidechain",
    ("SER", "HB2"): "serine_sidechain",
    ("SER", "HB3"): "serine_sidechain",
    ("SER", "OG"):  "hydroxyl",
    ("SER", "HG"):  "serine_sidechain",

    # ==========================
    # COO moieties (Asp/Glu)
    # ==========================
    ("ASP", "OD1"): "COO",
    ("ASP", "OD2"): "COO",
    ("ASP", "CG"):  "COO",
    ("GLU", "OE1"): "COO",
    ("GLU", "OE2"): "COO",
    ("GLU", "CD"):  "COO",

    # ==========================
    # Imidazole moiety (His)
    # ==========================
    ("HIS", "ND1"): "imidazole",
    ("HIS", "NE2"): "imidazole",
    ("HIS", "CG"):  "imidazole",
    ("HIS", "CD2"): "imidazole",
    ("HIS", "CE1"): "imidazole",

    # ==========================
    # Termini / backbone commons
    # ==========================
    ("ANY", "OXT"): "c_terminus",  # C-terminal oxygen
    ("ANY", "O"):   "backbone",
    ("ANY", "N"):   "backbone",

    # ==========================
    # Specific functional groups
    # ==========================
    ("TYR", "OH"): "phenol",
    ("CYS", "SG"): "thiol",
    ("LYS", "NZ"): "amine",
    ("ARG", "NE"): "guanidinium",
    ("ARG", "NH1"): "guanidinium",
    ("ARG", "NH2"): "guanidinium",
    ("ASN", "OD1"): "amide",
    ("ASN", "ND2"): "amide",
    ("GLN", "OE1"): "amide",
    ("GLN", "NE2"): "amide",
    ("MET", "SD"):  "thioether",

    # ======
    # Water
    # ======
    ("HOH", "O"): "water",
    ("WAT", "H1"): "water",
    ("WAT", "H2"): "water",

    # ==========
    # O2 ligand
    # ==========
    ("OY1", "O1"): "oxy",
    ("OY1", "O2"): "oxy",

    # ==========
    # ARG (full)
    # ==========
    ("ARG", "CB"):   "arginine_sidechain",
    ("ARG", "HB2"):  "arginine_sidechain",
    ("ARG", "HB3"):  "arginine_sidechain",
    ("ARG", "CG"):   "arginine_sidechain",
    ("ARG", "HG2"):  "arginine_sidechain",
    ("ARG", "HG3"):  "arginine_sidechain",
    ("ARG", "CD"):   "arginine_sidechain",
    ("ARG", "HD2"):  "arginine_sidechain",
    ("ARG", "HD3"):  "arginine_sidechain",
    ("ARG", "CZ"):   "arginine_sidechain",
    ("ARG", "HE"):   "arginine_sidechain",
    ("ARG", "HH11"): "arginine_sidechain",
    ("ARG", "HH12"): "arginine_sidechain",
    ("ARG", "HH21"): "arginine_sidechain",
    ("ARG", "HH22"): "arginine_sidechain",

    # ==========
    # ASP (full)
    # ==========
    ("ASP", "CB"):  "aspartate_sidechain",
    ("ASP", "HB2"): "aspartate_sidechain",
    ("ASP", "HB3"): "aspartate_sidechain",

    # ===================
    # HM1/HM2 (heme sets)
    # ===================
    ("HM1", "CHA"): "heme_por",
    ("HM1", "CHB"): "heme_por",
    ("HM1", "CHC"): "heme_por",
    ("HM1", "CHD"): "heme_por",
    ("HM1", "NA"):  "heme_por",
    ("HM1", "C1A"): "heme_por",
    ("HM1", "C2A"): "heme_por",
    ("HM1", "C3A"): "heme_por",
    ("HM1", "C4A"): "heme_por",
    ("HM1", "CMA"): "heme_por",
    ("HM1", "CAA"): "heme_por",
    ("HM1", "CBA"): "heme_por",
    ("HM1", "CGA"): "heme_por",
    ("HM1", "O1A"): "heme_por",
    ("HM1", "O2A"): "heme_por",
    ("HM1", "NB"):  "heme_por",
    ("HM1", "C1B"): "heme_por",
    ("HM1", "C2B"): "heme_por",
    ("HM1", "C3B"): "heme_por",
    ("HM1", "C4B"): "heme_por",
    ("HM1", "CMB"): "heme_por",
    ("HM1", "CAB"): "heme_por",
    ("HM1", "CBB"): "heme_por",
    ("HM1", "NC"):  "heme_por",
    ("HM1", "C1C"): "heme_por",
    ("HM1", "C2C"): "heme_por",
    ("HM1", "C3C"): "heme_por",
    ("HM1", "C4C"): "heme_por",
    ("HM1", "CMC"): "heme_por",
    ("HM1", "CAC"): "heme_por",
    ("HM1", "CBC"): "heme_por",
    ("HM1", "ND"):  "heme_por",
    ("HM1", "C1D"): "heme_por",
    ("HM1", "C2D"): "heme_por",
    ("HM1", "C3D"): "heme_por",
    ("HM1", "C4D"): "heme_por",
    ("HM1", "CMD"): "heme_por",
    ("HM1", "CAD"): "heme_por",
    ("HM1", "CBD"): "heme_por",
    ("HM1", "CGD"): "heme_por",
    ("HM1", "O1D"): "heme_por",
    ("HM1", "O2D"): "heme_por",
    ("HM1", "HMA1"): "heme_por",
    ("HM1", "HMA2"): "heme_por",
    ("HM1", "HMA3"): "heme_por",
    ("HM1", "HMB1"): "heme_por",
    ("HM1", "HMB2"): "heme_por",
    ("HM1", "HMB3"): "heme_por",
    ("HM1", "HMC1"): "heme_por",
    ("HM1", "HMC2"): "heme_por",
    ("HM1", "HMC3"): "heme_por",
    ("HM1", "HMD1"): "heme_por",
    ("HM1", "HMD2"): "heme_por",
    ("HM1", "HMD3"): "heme_por",
    ("HM1", "HBB1"): "heme_por",
    ("HM1", "HBB2"): "heme_por",
    ("HM1", "HBC1"): "heme_por",
    ("HM1", "HBC2"): "heme_por",
    ("HM1", "HBA1"): "heme_por",
    ("HM1", "HBA2"): "heme_por",
    ("HM1", "HAA1"): "heme_por",
    ("HM1", "HAA2"): "heme_por",
    ("HM1", "HBD1"): "heme_por",
    ("HM1", "HBD2"): "heme_por",
    ("HM1", "HAD1"): "heme_por",
    ("HM1", "HAD2"): "heme_por",
    ("HM1", "HHA"):  "heme_por",
    ("HM1", "HHB"):  "heme_por",
    ("HM1", "HHC"):  "heme_por",
    ("HM1", "HHD"):  "heme_por",
    ("HM1", "HAB"):  "heme_por",
    ("HM1", "HAC"):  "heme_por",

    ("HM2", "FE"): "heme_fe",

    # ===================
    # HEM (heme 3-letter)
    # ===================
    ("HEM", "CHA"): "heme_por",
    ("HEM", "CHB"): "heme_por",
    ("HEM", "CHC"): "heme_por",
    ("HEM", "CHD"): "heme_por",
    ("HEM", "NA"):  "heme_por",
    ("HEM", "C1A"): "heme_por",
    ("HEM", "C2A"): "heme_por",
    ("HEM", "C3A"): "heme_por",
    ("HEM", "C4A"): "heme_por",
    ("HEM", "CMA"): "heme_por",
    ("HEM", "CAA"): "heme_por",
    ("HEM", "CBA"): "heme_por",
    ("HEM", "CGA"): "heme_por",
    ("HEM", "O1A"): "heme_por",
    ("HEM", "O2A"): "heme_por",
    ("HEM", "NB"):  "heme_por",
    ("HEM", "C1B"): "heme_por",
    ("HEM", "C2B"): "heme_por",
    ("HEM", "C3B"): "heme_por",
    ("HEM", "C4B"): "heme_por",
    ("HEM", "CMB"): "heme_por",
    ("HEM", "CAB"): "heme_por",
    ("HEM", "CBB"): "heme_por",
    ("HEM", "NC"):  "heme_por",
    ("HEM", "C1C"): "heme_por",
    ("HEM", "C2C"): "heme_por",
    ("HEM", "C3C"): "heme_por",
    ("HEM", "C4C"): "heme_por",
    ("HEM", "CMC"): "heme_por",
    ("HEM", "CAC"): "heme_por",
    ("HEM", "CBC"): "heme_por",
    ("HEM", "ND"):  "heme_por",
    ("HEM", "C1D"): "heme_por",
    ("HEM", "C2D"): "heme_por",
    ("HEM", "C3D"): "heme_por",
    ("HEM", "C4D"): "heme_por",
    ("HEM", "CMD"): "heme_por",
    ("HEM", "CAD"): "heme_por",
    ("HEM", "CBD"): "heme_por",
    ("HEM", "CGD"): "heme_por",
    ("HEM", "O1D"): "heme_por",
    ("HEM", "O2D"): "heme_por",
    ("HEM", "HMA1"): "heme_por",
    ("HEM", "HMA2"): "heme_por",
    ("HEM", "HMA3"): "heme_por",
    ("HEM", "HMB1"): "heme_por",
    ("HEM", "HMB2"): "heme_por",
    ("HEM", "HMB3"): "heme_por",
    ("HEM", "HMC1"): "heme_por",
    ("HEM", "HMC2"): "heme_por",
    ("HEM", "HMC3"): "heme_por",
    ("HEM", "HMD1"): "heme_por",
    ("HEM", "HMD2"): "heme_por",
    ("HEM", "HMD3"): "heme_por",
    ("HEM", "HBB1"): "heme_por",
    ("HEM", "HBB2"): "heme_por",
    ("HEM", "HBC1"): "heme_por",
    ("HEM", "HBC2"): "heme_por",
    ("HEM", "HBA1"): "heme_por",
    ("HEM", "HBA2"): "heme_por",
    ("HEM", "HAA1"): "heme_por",
    ("HEM", "HAA2"): "heme_por",
    ("HEM", "HBD1"): "heme_por",
    ("HEM", "HBD2"): "heme_por",
    ("HEM", "HAD1"): "heme_por",
    ("HEM", "HAD2"): "heme_por",
    ("HEM", "HHA"):  "heme_por",
    ("HEM", "HHB"):  "heme_por",
    ("HEM", "HHC"):  "heme_por",
    ("HEM", "HHD"):  "heme_por",
    ("HEM", "HAB"):  "heme_por",
    ("HEM", "HAC"):  "heme_por",
    ("HEM", "FE"):   "heme_fe",

    # ===============
    # Sulfur dioxide
    # ===============
    ("SO2", "S"):  "sulfur_dioxide",
    ("SO2", "O1"): "sulfur_dioxide",
    ("SO2", "O2"): "sulfur_dioxide",

    # ==========================
    # FeMo-cofactor (ICS) atoms
    # ==========================
    ("ICS", "MO"):  "femoco_mo",
    ("ICS", "FE1"): "femoco_fe",
    ("ICS", "FE2"): "femoco_fe",
    ("ICS", "FE3"): "femoco_fe",
    ("ICS", "FE4"): "femoco_fe",
    ("ICS", "FE5"): "femoco_fe",
    ("ICS", "FE6"): "femoco_fe",
    ("ICS", "FE7"): "femoco_fe",
    ("ICS", "S1A"): "femoco_s",
    ("ICS", "S2B"): "femoco_s",
    ("ICS", "S3A"): "femoco_s",
    ("ICS", "S3B"): "femoco_s",
    ("ICS", "S4A"): "femoco_s",
    ("ICS", "S4B"): "femoco_s",
    ("ICS", "S5A"): "femoco_s",
    ("ICS", "S5B"): "femoco_s",
    ("ICS", "CX"):  "femoco_carbide",
    ("ICS", "O1"):  "femoco_o",
    ("ICS", "O2"):  "femoco_o",
    ("ICS", "O3"):  "femoco_o",
    ("ICS", "N1"):  "femoco_n",
    ("ICS", "N2"):  "femoco_n",

    # ==================
    # Homocitrate (HCA)
    # ==================
    ("HCA", "C1"): "hca_c",
    ("HCA", "C2"): "hca_c",
    ("HCA", "C3"): "hca_c",
    ("HCA", "C4"): "hca_c",
    ("HCA", "C5"): "hca_c",
    ("HCA", "C6"): "hca_c",
    ("HCA", "O1"): "hca_o",
    ("HCA", "O2"): "hca_o",
    ("HCA", "O3"): "hca_o",
    ("HCA", "O4"): "hca_o",
    ("HCA", "O5"): "hca_o",
    ("HCA", "O6"): "hca_o",
    ("HCA", "H1"): "hca_h",
    ("HCA", "H2"): "hca_h",
    ("HCA", "H3"): "hca_h",
    ("HCA", "H4"): "hca_h",
    ("HCA", "H5"): "hca_h",
    ("HCA", "H6"): "hca_h",

    # =====================================
    # Oxygen-Evolving Cluster (PSII) – OEX
    # =====================================
    ("OEX", "O1"):  "oxygen_evolving_cluster",
    ("OEX", "CA1"): "oxygen_evolving_cluster",
    ("OEX", "MN1"): "oxygen_evolving_cluster",
    ("OEX", "O2"):  "oxygen_evolving_cluster",
    ("OEX", "MN2"): "oxygen_evolving_cluster",
    ("OEX", "O3"):  "oxygen_evolving_cluster",
    ("OEX", "MN3"): "oxygen_evolving_cluster",
    ("OEX", "O4"):  "oxygen_evolving_cluster",
    ("OEX", "MN4"): "oxygen_evolving_cluster",
    ("OEX", "O5"):  "oxygen_evolving_cluster",
}


    
    
    









    
    









# ---------------------------------------------------------
# Moiety utilities
# ---------------------------------------------------------
def get_moiety_atoms(
    atom: Atom.Atom,
    residue: Residue.Residue,
    exclude_residue: Optional[Residue.Residue] = None
) -> List[Atom.Atom]:
    """
    Find all atoms in the same residue that belong to the same chemical moiety
    as the given atom. Handles special cases for O/OXT and water.

    Parameters
    ----------
    atom : Bio.PDB.Atom.Atom
        Atom object.
    residue : Bio.PDB.Residue.Residue
        Residue containing the atom.
    exclude_residue : Residue, optional
        If provided, suppress warnings for this residue.

    Returns
    -------
    List[Atom.Atom]
        List of atoms in the same moiety as `atom`.
    """
    residue_name = residue.get_resname()
    atom_name = atom.get_name()

    # --- special handling ---
    if atom_name == "OXT":
        if any(a.get_name() == "O" for a in residue):
            moiety = "c_terminus"
        else:
            moiety = None
    elif atom_name == "O":
        if any(a.get_name() == "OXT" for a in residue):
            moiety = "c_terminus"
        elif residue_name == "HOH":
            moiety = "water"
        else:
            moiety = "backbone"
    else:
        moiety = chemical_moieties.get((residue_name, atom_name))
        if moiety is None:
            moiety = chemical_moieties.get(("ANY", atom_name))

    # --- return atoms in same moiety ---
    if moiety:
        return [
            a for a in residue
            if (residue.get_resname(), a.get_name()) in chemical_moieties and
               chemical_moieties[(residue.get_resname(), a.get_name())] == moiety
        ]

    # warnings / skip rules
    if exclude_residue and residue == exclude_residue:
        return []
    if residue_name == "HOH":
        return []

    print(f"[WARNING] Atom {atom_name} in residue {residue_name} not found in MOIETY lookup table.")
    return []












# BOND LOOKUP TABLE!






# bond_lookup = {
    
    
    
    
#     # New entries for non-standard residues
#     "CU": [],           # Residue CU: no bonds defined (or define as needed)
#     "OXY": [("O1", "O2")],  # Residue OXY: bond between O1 and O2



#     # Standard amino acids
#     "ALA": [("N", "CA"), ("CA", "C"), ("C", "O")],
#     "GLY": [("N", "CA"), ("CA", "C"), ("C", "O")],
#     "ASP": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2"), ("CA", "C"), ("C", "O")],
#     "GLU": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2"), ("CA", "C"), ("C", "O")],
#     "HIS": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"), ("CD2", "NE2"), ("CA", "C"), ("C", "O"), ("CE1", "NE2")],
#     "HIE": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"), ("CD2", "NE2"), ("CA", "C"), ("C", "O"), ("CE1", "NE2")],
#     "HD1": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"), ("ND1", "CE1"), ("CD2", "NE2"), ("CA", "C"), ("C", "O"), ("CE1", "NE2")],
#     "CYS": [("N", "CA"), ("CA", "CB"), ("CB", "SG"), ("CA", "C"), ("C", "O")],
#     "LYS": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ"), ("CA", "C"), ("C", "O")],
#     "ARG": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"), ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2"), ("CA", "C"), ("C", "O")],
#     "SER": [("N", "CA"), ("CA", "CB"), ("CB", "OG"), ("CA", "C"), ("C", "O")],
#     "THR": [("N", "CA"), ("CA", "CB"), ("CB", "OG1"), ("CB", "CG2"), ("CA", "C"), ("C", "O")],
#     "TYR": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CZ", "OH"), ("CA", "C"), ("C", "O"), ("CZ", "CE2")],
#     "PHE": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CA", "C"), ("C", "O")],
#     "TRP": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CD1", "NE1"), ("CG", "CD2"), ("CD2", "CE2"), ("CD2", "CE3"), ("CE2", "CZ2"), ("CE3", "CZ3"), ("CZ2", "CH2"), ("CA", "C"), ("C", "O")],
#     "VAL": [("N", "CA"), ("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CA", "C"), ("C", "O")],
#     "LEU": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CA", "C"), ("C", "O")],
#     "ILE": [("N", "CA"), ("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1"), ("CA", "C"), ("C", "O")],
#     "PRO": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "N")],
#     "MET": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "SD"), ("SD", "CE")],
#     "ASN": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2"), ("CA", "C"), ("C", "O")],
#     "GLN": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2"), ("CA", "C"), ("C", "O")],

#     #SO3 (sulf)
#     "SO3": [("S", "O1"), ("S", "O2"), ("S", "O3")],
    
    
    
#     # HEM (HemeB) S
#     "HEM": [
#             ("CHA", "C1A"), ("C1A", "C2A"), ("CMA", "CDA"), ("CDA", "CEA"),
#             ("CEA", "O3A"), ("CEA", "O4A"), ("C2A", "C3A"), ("CAA", "CBA"), ("C3A", "C4A"),
#             ("CBA", "CCA"), ("CCA", "O1A"), ("CCA", "O2A"), ("C4A", "NA"), ("C1A", "NA"),
#             ("CHB", "C1B"), ("C1B", "C2B"), ("C2B", "CMB"), ("CMB", "CDB"), ("CDB", "CEB"),
#             ("CEB", "O3B"), ("CEB", "O4B"), ("C2B", "C3B"), ("C3B", "CAB"), ("CAB", "CBB"),
#             ("CBB", "CCB"), ("CCB", "O1B"), ("CCB", "O2B"), ("C4B", "NB"), ("C1B", "NB"), ("C3B", "C4B"),
#             ("CHC", "C1C"), ("C1C", "C2C"), ("C2C", "CDC"), ("CDC", "CEC"), ("CEC", "O3C"),
#             ("CEC", "O4C"), ("C2C", "C3C"), ("C3C", "CAC"), ("CAC", "CBC"), ("CBC", "CCC"),
#             ("CCC", "O1C"), ("CCC", "O2C"), ("C4C", "NC"), ("C1C", "NC"), ("C3C", "C4C"),
#             ("CHD", "C1D"), ("C1D", "C2D"), ("C2D", "CAD"), ("CAD", "CBD"), ("CBD", "CCD"),
#             ("CCD", "O1D"), ("CCD", "O2D"), ("C2D", "C3D"), ("C3D", "CDD"), ("CDD", "CED"),
#             ("CED", "O3D"), ("CED", "O4D"), ("C4D", "ND"), ("C1D", "ND"), ("C3D", "C4D"),
#             ("FE", "NA"), ("FE", "NB"), ("FE", "NC"), ("FE", "ND"),  # FE-NA bonds / Central iron and surrounding rings
#             ("CHA", "C4D"),  # CHA to C4D
#             ("CHD", "C4C"),  # CHD to C4C
#             ("CHC", "C4B"),  # CHC to C4B
#             ("CHB", "C4A"),   # CHB to C4A
#             ("CMD", "C2D"),
#             ("C2A", "CAA"),
#             ("CMC", "C2C"),
#             ("CBD", "CGD"),
#             ("CGD", "O2D"),
#             ("CGD", "O1D"),
#             ("CBA", "CGA"),
#             ("C3A", "CMA"),
#             ("O2A", "CGA"),
#             ("CGA", "O1A")],
          
#     # HEM (HemeB) S
#     "HEA": [
#             ("CHA", "C1A"), ("C1A", "C2A"), ("CMA", "CDA"), ("CDA", "CEA"),
#             ("CEA", "O3A"), ("CEA", "O4A"), ("C2A", "C3A"), ("CAA", "CBA"), ("C3A", "C4A"),
#             ("CBA", "CCA"), ("CCA", "O1A"), ("CCA", "O2A"), ("C4A", "NA"), ("C1A", "NA"),
#             ("CHB", "C1B"), ("C1B", "C2B"), ("C2B", "CMB"), ("CMB", "CDB"), ("CDB", "CEB"),
#             ("CEB", "O3B"), ("CEB", "O4B"), ("C2B", "C3B"), ("C3B", "CAB"), ("CAB", "CBB"),
#             ("CBB", "CCB"), ("CCB", "O1B"), ("CCB", "O2B"), ("C4B", "NB"), ("C1B", "NB"), ("C3B", "C4B"),
#             ("CHC", "C1C"), ("C1C", "C2C"), ("C2C", "CDC"), ("CDC", "CEC"), ("CEC", "O3C"),
#             ("CEC", "O4C"), ("C2C", "C3C"), ("C3C", "CAC"), ("CAC", "CBC"), ("CBC", "CCC"),
#             ("CCC", "O1C"), ("CCC", "O2C"), ("C4C", "NC"), ("C1C", "NC"), ("C3C", "C4C"),
#             ("CHD", "C1D"), ("C1D", "C2D"), ("C2D", "CAD"), ("CAD", "CBD"), ("CBD", "CCD"),
#             ("CCD", "O1D"), ("CCD", "O2D"), ("C2D", "C3D"), ("C3D", "CDD"), ("CDD", "CED"),
#             ("CED", "O3D"), ("CED", "O4D"), ("C4D", "ND"), ("C1D", "ND"), ("C3D", "C4D"),
#             ("FE", "NA"), ("FE", "NB"), ("FE", "NC"), ("FE", "ND"),  # FE-NA bonds / Central iron and surrounding rings
#             ("CHA", "C4D"),  # CHA to C4D
#             ("CHD", "C4C"),  # CHD to C4C
#             ("CHC", "C4B"),  # CHC to C4B
#             ("CHB", "C4A"),   # CHB to C4A
#             ("CMD", "C2D"),
#             ("C2A", "CAA"),
#             ("CMC", "C2C"),
#             ("CBD", "CGD"),
#             ("CGD", "O2D"),
#             ("CGD", "O1D"),
#             ("CBA", "CGA"),
#             ("C3A", "CMA"),
#             ("O2A", "CGA"),
#             ("CGA", "O1A"), ("C3B", "C11"), ("C11", "C12"), ("C11","O11"), ("CMA","OMA")],
    
    
    
    
    
    
#     # HEM (HemeB) S
#     "HM1": [
#             ("CHA", "C1A"), ("C1A", "C2A"), ("CMA", "CDA"), ("CDA", "CEA"),
#             ("CEA", "O3A"), ("CEA", "O4A"), ("C2A", "C3A"), ("CAA", "CBA"), ("C3A", "C4A"),
#             ("CBA", "CCA"), ("CCA", "O1A"), ("CCA", "O2A"), ("C4A", "NA"), ("C1A", "NA"),
#             ("CHB", "C1B"), ("C1B", "C2B"), ("C2B", "CMB"), ("CMB", "CDB"), ("CDB", "CEB"),
#             ("CEB", "O3B"), ("CEB", "O4B"), ("C2B", "C3B"), ("C3B", "CAB"), ("CAB", "CBB"),
#             ("CBB", "CCB"), ("CCB", "O1B"), ("CCB", "O2B"), ("C4B", "NB"), ("C1B", "NB"), ("C3B", "C4B"),
#             ("CHC", "C1C"), ("C1C", "C2C"), ("C2C", "CDC"), ("CDC", "CEC"), ("CEC", "O3C"),
#             ("CEC", "O4C"), ("C2C", "C3C"), ("C3C", "CAC"), ("CAC", "CBC"), ("CBC", "CCC"),
#             ("CCC", "O1C"), ("CCC", "O2C"), ("C4C", "NC"), ("C1C", "NC"), ("C3C", "C4C"),
#             ("CHD", "C1D"), ("C1D", "C2D"), ("C2D", "CAD"), ("CAD", "CBD"), ("CBD", "CCD"),
#             ("CCD", "O1D"), ("CCD", "O2D"), ("C2D", "C3D"), ("C3D", "CDD"), ("CDD", "CED"),
#             ("CED", "O3D"), ("CED", "O4D"), ("C4D", "ND"), ("C1D", "ND"), ("C3D", "C4D"),
#             ("FE", "NA"), ("FE", "NB"), ("FE", "NC"), ("FE", "ND"),  # FE-NA bonds / Central iron and surrounding rings
#             ("CHA", "C4D"),  # CHA to C4D
#             ("CHD", "C4C"),  # CHD to C4C
#             ("CHC", "C4B"),  # CHC to C4B
#             ("CHB", "C4A"),   # CHB to C4A
#             ("CMD", "C2D"),
#             ("C2A", "CAA"),
#             ("CMC", "C2C"),
#             ("CBD", "CGD"),
#             ("CGD", "O2D"),
#             ("CGD", "O1D"),
#             ("CBA", "CGA"),
#             ("C3A", "CMA"),
#             ("O2A", "CGA"),
#             ("CGA", "O1A"), ("C3B", "C11"), ("C11", "C12"), ("C11","O11"), ("CMA","OMA")],
          
    
        
    
#     # HEM (HemeB) S
#     "HM2": [
#             ("CHA", "C1A"), ("C1A", "C2A"), ("CMA", "CDA"), ("CDA", "CEA"),
#             ("CEA", "O3A"), ("CEA", "O4A"), ("C2A", "C3A"), ("CAA", "CBA"), ("C3A", "C4A"),
#             ("CBA", "CCA"), ("CCA", "O1A"), ("CCA", "O2A"), ("C4A", "NA"), ("C1A", "NA"),
#             ("CHB", "C1B"), ("C1B", "C2B"), ("C2B", "CMB"), ("CMB", "CDB"), ("CDB", "CEB"),
#             ("CEB", "O3B"), ("CEB", "O4B"), ("C2B", "C3B"), ("C3B", "CAB"), ("CAB", "CBB"),
#             ("CBB", "CCB"), ("CCB", "O1B"), ("CCB", "O2B"), ("C4B", "NB"), ("C1B", "NB"), ("C3B", "C4B"),
#             ("CHC", "C1C"), ("C1C", "C2C"), ("C2C", "CDC"), ("CDC", "CEC"), ("CEC", "O3C"),
#             ("CEC", "O4C"), ("C2C", "C3C"), ("C3C", "CAC"), ("CAC", "CBC"), ("CBC", "CCC"),
#             ("CCC", "O1C"), ("CCC", "O2C"), ("C4C", "NC"), ("C1C", "NC"), ("C3C", "C4C"),
#             ("CHD", "C1D"), ("C1D", "C2D"), ("C2D", "CAD"), ("CAD", "CBD"), ("CBD", "CCD"),
#             ("CCD", "O1D"), ("CCD", "O2D"), ("C2D", "C3D"), ("C3D", "CDD"), ("CDD", "CED"),
#             ("CED", "O3D"), ("CED", "O4D"), ("C4D", "ND"), ("C1D", "ND"), ("C3D", "C4D"),
#             ("FE", "NA"), ("FE", "NB"), ("FE", "NC"), ("FE", "ND"),  # FE-NA bonds / Central iron and surrounding rings
#             ("CHA", "C4D"),  # CHA to C4D
#             ("CHD", "C4C"),  # CHD to C4C
#             ("CHC", "C4B"),  # CHC to C4B
#             ("CHB", "C4A"),   # CHB to C4A
#             ("CMD", "C2D"),
#             ("C2A", "CAA"),
#             ("CMC", "C2C"),
#             ("CBD", "CGD"),
#             ("CGD", "O2D"),
#             ("CGD", "O1D"),
#             ("CBA", "CGA"),
#             ("C3A", "CMA"),
#             ("O2A", "CGA"),
#             ("CGA", "O1A"), ("C3B", "C11"), ("C11", "C12"), ("C11","O11"), ("CMA","OMA")],
          
    
    
    
    
#     # SRM (siroheme) S
#     "SRM": [
#             ("CHA", "C1A"), ("C1A", "C2A"), ("C2A", "CMA"), ("CMA", "CDA"), ("CDA", "CEA"),
#             ("CEA", "O3A"), ("CEA", "O4A"), ("C2A", "C3A"), ("C3A", "CAA"), ("CAA", "CBA"), ("C3A", "C4A"),
#             ("CBA", "CCA"), ("CCA", "O1A"), ("CCA", "O2A"), ("C4A", "NA"), ("C1A", "NA"),
#             ("CHB", "C1B"), ("C1B", "C2B"), ("C2B", "CMB"), ("CMB", "CDB"), ("CDB", "CEB"),
#             ("CEB", "O3B"), ("CEB", "O4B"), ("C2B", "C3B"), ("C3B", "CAB"), ("CAB", "CBB"),
#             ("CBB", "CCB"), ("CCB", "O1B"), ("CCB", "O2B"), ("C4B", "NB"), ("C1B", "NB"), ("C3B", "C4B"),
#             ("CHC", "C1C"), ("C1C", "C2C"), ("C2C", "CDC"), ("CDC", "CEC"), ("CEC", "O3C"),
#             ("CEC", "O4C"), ("C2C", "C3C"), ("C3C", "CAC"), ("CAC", "CBC"), ("CBC", "CCC"),
#             ("CCC", "O1C"), ("CCC", "O2C"), ("C4C", "NC"), ("C1C", "NC"), ("C3C", "C4C"),
#             ("CHD", "C1D"), ("C1D", "C2D"), ("C2D", "CAD"), ("CAD", "CBD"), ("CBD", "CCD"),
#             ("CCD", "O1D"), ("CCD", "O2D"), ("C2D", "C3D"), ("C3D", "CDD"), ("CDD", "CED"),
#             ("CED", "O3D"), ("CED", "O4D"), ("C4D", "ND"), ("C1D", "ND"), ("C3D", "C4D"),
#             ("FE", "NA"), ("FE", "NB"), ("FE", "NC"), ("FE", "ND"),  # FE-NA bonds / Central iron and surrounding rings
#             ("CHA", "C4D"),  # CHA to C4D
#             ("CHD", "C4C"),  # CHD to C4C
#             ("CHC", "C4B"),  # CHC to C4B
#             ("CHB", "C4A")],   # CHB to C4A
            
    
#     # SF4 (4Fe4S) and SO3 (sulf)
#     "SF4": [("FE1", "S1"), ("FE2", "S1"), ("FE1", "S2"), ("FE3", "S2"), ("FE2", "S3"), ("FE4", "S3"), ("FE3", "S4"),
#             ("FE4", "S1"), 
#             ("FE4", "S2"), 
#             ("FE1", "S3"), 
#             ("FE2", "S4"), 
#             ("FE3", "S1")],
            
    
    
    
#     # Oxygen evolving cluster:
#     "OEX": [("MN4", "O4"),
#             ("MN4", "O5"),
#             ("MN3", "O2"),
#             ("MN3", "O5"),
#             ("MN3", "O3"),
#             ("MN3", "O4"),
#             ("MN2", "O3"),
#             ("MN2", "O2"),
#             ("MN2", "O1"),
#             ("MN1", "O5"),
#             ("MN1", "O3"),
#             ("MN1", "O1"),
#             ("CA1", "O1"),
#             ("CA1", "O2"),
#             ("CA1", "O5")],
    
    
#     "OY1": [("O1", "O2")],



#     "SO2": [("S", "O1"),
#             ("S", "O2")],




#     "CU1": [("C1", "CU")],

    


#     "ICS": [
#         # Central carbide bonds
#         ("CX", "FE1"), ("CX", "FE2"), ("CX", "FE3"), ("CX", "FE4"), ("CX", "FE5"), ("CX", "FE6"),
#         # Fe–S belts (each Fe typically linked to 3 sulfurs)
#         ("FE1", "S1A"), ("FE1", "S2A"), ("FE1", "S4A"),
#         ("FE2", "S2A"), ("FE2", "S2B"), ("FE2", "S1A"),
#         ("FE3", "S5A"), ("FE3", "S4B"), ("FE3", "S2A"),
#         ("FE4", "S1A"), ("FE4", "S3A"), ("FE4", "S4A"),
#         ("FE5", "S4B"), ("FE5", "S3A"), ("FE5", "S1B"),
#         ("FE6", "S3B"), ("FE6", "S2B"), ("FE6", "S1B"),
#         ("FE7", "S3B"), ("FE7", "S4B"), ("FE7", "S5A"),
#         # Mo connections
#         ("MO1", "S4B"), ("MO1", "S1B"), ("MO1", "S3B"),
#     ],


#     "HCA": [
#         # Backbone chain
#         ("C1", "C2"), ("C2", "C3"), ("C3", "C4"), ("C3", "C7"), ("C4", "C5"), ("C5", "C6"),

#         # Carboxyl groups
#         ("C1", "O1"), ("C1", "O2"),       # α-carboxyl
#         ("C3", "O7"),
#         ("C7", "O6"), ("C7", "O5"),
#         ("C6", "O4"), ("C6", "O3"),        # terminal carboxyl 
#     ],




    
# }



# BOND LOOKUP TABLE!



bond_lookup = {

    # --- Non-standard / small molecules ---
    "CU": [],                                # single-atom cofactor
    "OXY": [("O1", "O2")],                   # O2
    "OY1": [("O1", "O2")],                   # alt O2 label
    "SO2": [("S", "O1"), ("S", "O2")],
    "SO3": [("S", "O1"), ("S", "O2"), ("S", "O3")],
    "CU1": [("C1", "CU")],                   # site-specific label

    # --- Standard amino acids with termini support (N-H*, C-OXT) ---

    "ALA": [
        ("N","CA"), ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "GLY": [
        ("N","CA"), ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "ASP": [
        ("N","CA"), ("CA","CB"), ("CB","CG"),
        ("CG","OD1"), ("CG","OD2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "GLU": [
        ("N","CA"), ("CA","CB"), ("CB","CG"), ("CG","CD"),
        ("CD","OE1"), ("CD","OE2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "HIS": [
        ("N","CA"), ("CA","CB"), ("CB","CG"),
        ("CG","ND1"), ("CG","CD2"), ("ND1","CE1"), ("CD2","NE2"), ("CE1","NE2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "HIE": [
        ("N","CA"), ("CA","CB"), ("CB","CG"),
        ("CG","ND1"), ("CG","CD2"), ("ND1","CE1"), ("CD2","NE2"), ("CE1","NE2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "HD1": [
        ("N","CA"), ("CA","CB"), ("CB","CG"),
        ("CG","ND1"), ("CG","CD2"), ("ND1","CE1"), ("CD2","NE2"), ("CE1","NE2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "CYS": [
        ("N","CA"), ("CA","CB"), ("CB","SG"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "LYS": [
        ("N","CA"), ("CA","CB"), ("CB","CG"), ("CG","CD"), ("CD","CE"), ("CE","NZ"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "ARG": [
        ("N","CA"), ("CA","CB"), ("CB","CG"), ("CG","CD"),
        ("CD","NE"), ("NE","CZ"), ("CZ","NH1"), ("CZ","NH2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "SER": [
        ("N","CA"), ("CA","CB"), ("CB","OG"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "THR": [
        ("N","CA"), ("CA","CB"), ("CB","OG1"), ("CB","CG2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "TYR": [
        ("N","CA"), ("CA","CB"), ("CB","CG"),
        ("CG","CD1"), ("CG","CD2"), ("CD1","CE1"), ("CD2","CE2"),
        ("CE1","CZ"), ("CE2","CZ"), ("CZ","OH"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "PHE": [
        ("N","CA"), ("CA","CB"), ("CB","CG"),
        ("CG","CD1"), ("CG","CD2"), ("CD1","CE1"), ("CD2","CE2"),
        ("CE1","CZ"), ("CE2","CZ"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "TRP": [
        ("N","CA"), ("CA","CB"), ("CB","CG"),
        ("CG","CD1"), ("CD1","NE1"), ("CG","CD2"),
        ("CD2","CE2"), ("CD2","CE3"),
        ("CE2","CZ2"), ("CE3","CZ3"), ("CZ2","CH2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "VAL": [
        ("N","CA"), ("CA","CB"), ("CB","CG1"), ("CB","CG2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "LEU": [
        ("N","CA"), ("CA","CB"), ("CB","CG"),
        ("CG","CD1"), ("CG","CD2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "ILE": [
        ("N","CA"), ("CA","CB"), ("CB","CG1"), ("CB","CG2"),
        ("CG1","CD1"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    # PRO: ring closes to backbone N; usually no backbone H, but N–H* bonds are harmless if absent
    "PRO": [
        ("N","CA"), ("CA","CB"), ("CB","CG"), ("CG","CD"), ("CD","N"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),  # will be ignored if not present
        ("C","OXT")
    ],

    "MET": [
        ("N","CA"), ("CA","CB"), ("CB","CG"), ("CG","SD"), ("SD","CE"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "ASN": [
        ("N","CA"), ("CA","CB"), ("CB","CG"),
        ("CG","OD1"), ("CG","ND2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],

    "GLN": [
        ("N","CA"), ("CA","CB"), ("CB","CG"), ("CG","CD"),
        ("CD","OE1"), ("CD","NE2"),
        ("CA","C"), ("C","O"),
        ("N","H"), ("N","H1"), ("N","H2"), ("N","H3"),
        ("C","OXT")
    ],




    # --- Heme macrocycles & related cofactors (kept as provided) ---
    "HEM": [
        ("CHA","C1A"), ("C1A","C2A"), ("CMA","CDA"), ("CDA","CEA"),
        ("CEA","O3A"), ("CEA","O4A"), ("C2A","C3A"), ("CAA","CBA"), ("C3A","C4A"),
        ("CBA","CCA"), ("CCA","O1A"), ("CCA","O2A"), ("C4A","NA"), ("C1A","NA"),
        ("CHB","C1B"), ("C1B","C2B"), ("C2B","CMB"), ("CMB","CDB"), ("CDB","CEB"),
        ("CEB","O3B"), ("CEB","O4B"), ("C2B","C3B"), ("C3B","CAB"), ("CAB","CBB"),
        ("CBB","CCB"), ("CCB","O1B"), ("CCB","O2B"), ("C4B","NB"), ("C1B","NB"), ("C3B","C4B"),
        ("CHC","C1C"), ("C1C","C2C"), ("C2C","CDC"), ("CDC","CEC"), ("CEC","O3C"),
        ("CEC","O4C"), ("C2C","C3C"), ("C3C","CAC"), ("CAC","CBC"), ("CBC","CCC"),
        ("CCC","O1C"), ("CCC","O2C"), ("C4C","NC"), ("C1C","NC"), ("C3C","C4C"),
        ("CHD","C1D"), ("C1D","C2D"), ("C2D","CAD"), ("CAD","CBD"), ("CBD","CCD"),
        ("CCD","O1D"), ("CCD","O2D"), ("C2D","C3D"), ("C3D","CDD"), ("CDD","CED"),
        ("CED","O3D"), ("CED","O4D"), ("C4D","ND"), ("C1D","ND"), ("C3D","C4D"),
        ("FE","NA"), ("FE","NB"), ("FE","NC"), ("FE","ND"),
        ("CHA","C4D"), ("CHD","C4C"), ("CHC","C4B"), ("CHB","C4A"),
        ("CMD","C2D"), ("C2A","CAA"), ("CMC","C2C"),
        ("CBD","CGD"), ("CGD","O2D"), ("CGD","O1D"),
        ("CBA","CGA"), ("C3A","CMA"), ("O2A","CGA"), ("CGA","O1A")
    ],

    "HEA": [
        ("CHA","C1A"), ("C1A","C2A"), ("CMA","CDA"), ("CDA","CEA"),
        ("CEA","O3A"), ("CEA","O4A"), ("C2A","C3A"), ("CAA","CBA"), ("C3A","C4A"),
        ("CBA","CCA"), ("CCA","O1A"), ("CCA","O2A"), ("C4A","NA"), ("C1A","NA"),
        ("CHB","C1B"), ("C1B","C2B"), ("C2B","CMB"), ("CMB","CDB"), ("CDB","CEB"),
        ("CEB","O3B"), ("CEB","O4B"), ("C2B","C3B"), ("C3B","CAB"), ("CAB","CBB"),
        ("CBB","CCB"), ("CCB","O1B"), ("CCB","O2B"), ("C4B","NB"), ("C1B","NB"), ("C3B","C4B"),
        ("CHC","C1C"), ("C1C","C2C"), ("C2C","CDC"), ("CDC","CEC"), ("CEC","O3C"),
        ("CEC","O4C"), ("C2C","C3C"), ("C3C","CAC"), ("CAC","CBC"), ("CBC","CCC"),
        ("CCC","O1C"), ("CCC","O2C"), ("C4C","NC"), ("C1C","NC"), ("C3C","C4C"),
        ("CHD","C1D"), ("C1D","C2D"), ("C2D","CAD"), ("CAD","CBD"), ("CBD","CCD"),
        ("CCD","O1D"), ("CCD","O2D"), ("C2D","C3D"), ("C3D","CDD"), ("CDD","CED"),
        ("CED","O3D"), ("CED","O4D"), ("C4D","ND"), ("C1D","ND"), ("C3D","C4D"),
        ("FE","NA"), ("FE","NB"), ("FE","NC"), ("FE","ND"),
        ("CHA","C4D"), ("CHD","C4C"), ("CHC","C4B"), ("CHB","C4A"),
        ("CMD","C2D"), ("C2A","CAA"), ("CMC","C2C"),
        ("CBD","CGD"), ("CGD","O2D"), ("CGD","O1D"),
        ("CBA","CGA"), ("C3A","CMA"), ("O2A","CGA"), ("CGA","O1A"),
        ("C3B","C11"), ("C11","C12"), ("C11","O11"), ("CMA","OMA")
    ],

    "HM1": [
        ("CHA","C1A"), ("C1A","C2A"), ("C2A","CMA"), ("CMA","CDA"), ("CDA","CEA"),
        ("CEA","O3A"), ("CEA","O4A"), ("C2A","C3A"), ("C3A","CAA"), ("CAA","CBA"), ("C3A","C4A"),
        ("CBA","CCA"), ("CCA","O1A"), ("CCA","O2A"), ("C4A","NA"), ("C1A","NA"),
        ("CHB","C1B"), ("C1B","C2B"), ("C2B","CMB"), ("CMB","CDB"), ("CDB","CEB"),
        ("CEB","O3B"), ("CEB","O4B"), ("C2B","C3B"), ("C3B","CAB"), ("CAB","CBB"),
        ("CBB","CCB"), ("CCB","O1B"), ("CCB","O2B"), ("C4B","NB"), ("C1B","NB"), ("C3B","C4B"),
        ("CHC","C1C"), ("C1C","C2C"), ("C2C","CDC"), ("CDC","CEC"), ("CEC","O3C"),
        ("CEC","O4C"), ("C2C","C3C"), ("C3C","CAC"), ("CAC","CBC"), ("CBC","CCC"),
        ("CCC","O1C"), ("CCC","O2C"), ("C4C","NC"), ("C1C","NC"), ("C3C","C4C"),
        ("CHD","C1D"), ("C1D","C2D"), ("C2D","CAD"), ("CAD","CBD"), ("CBD","CCD"),
        ("CCD","O1D"), ("CCD","O2D"), ("C2D","C3D"), ("C3D","CDD"), ("CDD","CED"),
        ("CED","O3D"), ("CED","O4D"), ("C4D","ND"), ("C1D","ND"), ("C3D","C4D"),
        ("FE","NA"), ("FE","NB"), ("FE","NC"), ("FE","ND"),
        ("CHA","C4D"), ("CHD","C4C"), ("CHC","C4B"), ("CHB","C4A"),
        ("CMD","C2D"), ("C2A","CAA"), ("CMC","C2C"),
        ("CBD","CGD"), ("CGD","O2D"), ("CGD","O1D"),
        ("CBA","CGA"), ("C3A","CMA"), ("O2A","CGA"), ("CGA","O1A"),
        ("C3B","C11"), ("C11","C12"), ("C11","O11"), ("CMA","OMA")
    ],

    "HM2": [
        ("CHA","C1A"), ("C1A","C2A"), ("CMA","CDA"), ("CDA","CEA"),
        ("CEA","O3A"), ("CEA","O4A"), ("C2A","C3A"), ("CAA","CBA"), ("C3A","C4A"),
        ("CBA","CCA"), ("CCA","O1A"), ("CCA","O2A"), ("C4A","NA"), ("C1A","NA"),
        ("CHB","C1B"), ("C1B","C2B"), ("C2B","CMB"), ("CMB","CDB"), ("CDB","CEB"),
        ("CEB","O3B"), ("CEB","O4B"), ("C2B","C3B"), ("C3B","CAB"), ("CAB","CBB"),
        ("CBB","CCB"), ("CCB","O1B"), ("CCB","O2B"), ("C4B","NB"), ("C1B","NB"), ("C3B","C4B"),
        ("CHC","C1C"), ("C1C","C2C"), ("C2C","CDC"), ("CDC","CEC"), ("CEC","O3C"),
        ("CEC","O4C"), ("C2C","C3C"), ("C3C","CAC"), ("CAC","CBC"), ("CBC","CCC"),
        ("CCC","O1C"), ("CCC","O2C"), ("C4C","NC"), ("C1C","NC"), ("C3C","C4C"),
        ("CHD","C1D"), ("C1D","C2D"), ("C2D","CAD"), ("CAD","CBD"), ("CBD","CCD"),
        ("CCD","O1D"), ("CCD","O2D"), ("C2D","C3D"), ("C3D","CDD"), ("CDD","CED"),
        ("CED","O3D"), ("CED","O4D"), ("C4D","ND"), ("C1D","ND"), ("C3D","C4D"),
        ("FE","NA"), ("FE","NB"), ("FE","NC"), ("FE","ND"),
        ("CHA","C4D"), ("CHD","C4C"), ("CHC","C4B"), ("CHB","C4A"),
        ("CMD","C2D"), ("C2A","CAA"), ("CMC","C2C"),
        ("CBD","CGD"), ("CGD","O2D"), ("CGD","O1D"),
        ("CBA","CGA"), ("C3A","CMA"), ("O2A","CGA"), ("CGA","O1A"),
        ("C3B","C11"), ("C11","C12"), ("C11","O11"), ("CMA","OMA")
    ],

    "SRM": [
        ("CHA","C1A"), ("C1A","C2A"), ("C2A","CMA"), ("CMA","CDA"), ("CDA","CEA"),
        ("CEA","O3A"), ("CEA","O4A"), ("C2A","C3A"), ("C3A","CAA"), ("CAA","CBA"), ("C3A","C4A"),
        ("CBA","CCA"), ("CCA","O1A"), ("CCA","O2A"), ("C4A","NA"), ("C1A","NA"),
        ("CHB","C1B"), ("C1B","C2B"), ("C2B","CMB"), ("CMB","CDB"), ("CDB","CEB"),
        ("CEB","O3B"), ("CEB","O4B"), ("C2B","C3B"), ("C3B","CAB"), ("CAB","CBB"),
        ("CBB","CCB"), ("CCB","O1B"), ("CCB","O2B"), ("C4B","NB"), ("C1B","NB"), ("C3B","C4B"),
        ("CHC","C1C"), ("C1C","C2C"), ("C2C","CDC"), ("CDC","CEC"), ("CEC","O3C"),
        ("CEC","O4C"), ("C2C","C3C"), ("C3C","CAC"), ("CAC","CBC"), ("CBC","CCC"),
        ("CCC","O1C"), ("CCC","O2C"), ("C4C","NC"), ("C1C","NC"), ("C3C","C4C"),
        ("CHD","C1D"), ("C1D","C2D"), ("C2D","CAD"), ("CAD","CBD"), ("CBD","CCD"),
        ("CCD","O1D"), ("CCD","O2D"), ("C2D","C3D"), ("C3D","CDD"), ("CDD","CED"),
        ("CED","O3D"), ("CED","O4D"), ("C4D","ND"), ("C1D","ND"), ("C3D","C4D"),
        ("FE","NA"), ("FE","NB"), ("FE","NC"), ("FE","ND"),
        ("CHA","C4D"), ("CHD","C4C"), ("CHC","C4B"), ("CHB","C4A")
    ],

    # --- Fe-S cluster ---
    "SF4": [
        ("FE1","S1"), ("FE2","S1"),
        ("FE1","S2"), ("FE3","S2"),
        ("FE2","S3"), ("FE4","S3"),
        ("FE3","S4"), ("FE4","S1"),
        ("FE4","S2"), ("FE1","S3"), ("FE2","S4"), ("FE3","S1")
    ],

    # --- Oxygen evolving cluster ---
    "OEX": [
        ("MN4","O4"), ("MN4","O5"),
        ("MN3","O2"), ("MN3","O5"), ("MN3","O3"), ("MN3","O4"),
        ("MN2","O3"), ("MN2","O2"), ("MN2","O1"),
        ("MN1","O5"), ("MN1","O3"), ("MN1","O1"),
        ("CA1","O1"), ("CA1","O2"), ("CA1","O5")
    ],

    # --- FeMo-cofactor (ICS) ---
    "ICS": [
        # central carbide
        ("CX","FE1"), ("CX","FE2"), ("CX","FE3"), ("CX","FE4"), ("CX","FE5"), ("CX","FE6"),
        # Fe–S belt (typical connectivity)
        ("FE1","S1A"), ("FE1","S2A"), ("FE1","S4A"),
        ("FE2","S2A"), ("FE2","S2B"), ("FE2","S1A"),
        ("FE3","S5A"), ("FE3","S4B"), ("FE3","S2A"),
        ("FE4","S1A"), ("FE4","S3A"), ("FE4","S4A"),
        ("FE5","S4B"), ("FE5","S3A"), ("FE5","S1B"),
        ("FE6","S3B"), ("FE6","S2B"), ("FE6","S1B"),
        ("FE7","S3B"), ("FE7","S4B"), ("FE7","S5A"),
        # Mo connections (support MO and MO1 label variants)
        ("MO","S4B"), ("MO","S1B"), ("MO","S3B"),
        ("MO1","S4B"), ("MO1","S1B"), ("MO1","S3B"),
    ],




    # --- Homocitrate (HCA) ---
    # Uses atom names C1–C6 and O1–O6 as in your moieties table.
    # If your structure uses a C7/O7 label, replace the C5/C6 terminal pairs with (C7,"O5"/"O6") and add ("C3","C7").
    "HCA": [
        # carbon backbone
        ("C1","C2"), ("C2","C3"), ("C3","C4"), ("C4","C5"), ("C5","C6"),
        # carboxylate at C1 (alpha)
        ("C1","O1"), ("C1","O2"),
        # second carboxylate at C5 (variant: sometimes annotated at C6)
        ("C5","O3"), ("C5","O4"),
        # terminal carboxylate at C6
        ("C6","O5"), ("C6","O6"),
    ],
}


