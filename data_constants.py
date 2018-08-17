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
apolar_elements = ["C","S"]

polar_acceptor = {
    1:["O", "SD","OD1","OD2","OE1","OE2","O4'","O5'"],
    2:["DG.N7", "DG.O6","DG.N3","DC.N3","DA.N7","DA.N3","DA.N1","DT.O4","O4'","O5'"],
    3:["G.N7", "G.O6","G.N3","C.N3","A.N7","A.N3","A.N1","U.O4","O4'","O5'"]
}

polar_donor = {
    1:["N","NE","NZ","NH1","NH2","ND1","ND2","NE1","NE2"],
    2:["DG.N2", "DC.N4","DA.N6","DT.N3"],
    3:["G.N2", "C.N4","A.N6","U.N3"]
}

#ionic
pos_ats = ["NZ", "NE", "NH1", "NH2"]
neg_ats = ["OD1", "OD2", "OE1", "OE2","OP2"]

amide_res = {'ASN':['OD1','ND2'],'GLN':['OE1','NE2']}
amide_atoms = ['OD1','ND2','OE1','NE2']

chiral_res = {'THR':['OG1','CG2'], 'ILE':['CG1','CG2']}


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

