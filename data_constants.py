#
# Global data for structure manager
#

#=========================================================================
# Defaults. Don't touch unless knowing how
#=========================================================================

#Atom groups
metal_ats = [
    "MG", "MN", "MO", "Mg2", "Mn2", "Ca2", "ZN", "NI", "FE",
    "Zn2", "Ni2", "Fe2", "CO", "CU", "HG", "Co2", "Cu2",
    "CD", "AG", "AU", "Cd2", "CA", "Ca"
]
apolar_elements = ["C", "S"]

polar_acceptor = {
    '*': ['O', "O4'", "05'"],
    'MET': ['SD'],
    'ASP': ['OD1', 'OD2'],
    'ASN': ['OD1'],
    'GLU': ['OE1', 'OE2'],
    'GLN': ['OE1'],
    'DG': ['N7', 'O6', 'N3'],
    'DC': ['N3'],
    'DA': ['N7', 'N3', 'N1'],
    'DT': ['O4'],
    'G': ['N7', 'O6', 'N3'],
    'C': ['N3'],
    'A': ['N7', 'N3', 'N1'],
    'U': ['O4'],
}

polar_donor = {
    '*': ['N'],
    'ARG': ['NE', 'NH1', 'NH2'],
    'LYS': ['NZ'],
    'HIS': ['ND1', 'NE2'],
    'TRP': ['NE1'],
    'ASN': ['ND2'],
    'GLN': ['NE2'],
    'DG': ['N2'],
    'DC': ['N4'],
    'DA': ['N6'],
    'DT': ['N3'],
    'G': ['N2'],
    'C': ['N4'],
    'A': ['N6'],
    'U': ['N3']
}

#ionics
pos_ats = {
    '*' :[], 
    'LYS': ['NZ'],
    'ARG': ['NE', 'NH1', 'NH2']
}

neg_ats = {
    'ASP': ["OD1", "OD2"],
    'GLU': ["OE1", "OE2"],
    '*': ['OXT', 'OP1', "OP2"]
}

amide_res = {
    'ASN':['OD1', 'ND2'], 
    'GLN':['OE1', 'NE2']
}

amide_atoms = ['OD1', 'ND2', 'OE1', 'NE2']

chiral_res = {
    'THR':['OG1', 'CG2'], 
    'ILE':['CG1', 'CG2']
}


# Relevant Distances
SS_DIST = 2.5

CA_DIST = 3.8
CA_DIST_ERR = 1.

CLASH_DIST = {
    'severe': 1.,
    'apolar': 2.9,
    'donor': 3.1,
    'acceptor': 3.1,
    'positive': 3.5,
    'negative': 3.5
}

MAX_DIST = 4.5;

R_R_CUTOFF = 15

