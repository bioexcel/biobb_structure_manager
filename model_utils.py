"""
  Utility functions to manipulate structures
"""

import sys
import re
import math
import numpy as np
from numpy import sin, cos, pi, dot, clip, arccos
from numpy.linalg import norm

#Model Types
NMR=1
TRAJ=2
UNKNOWN=0

# TODO: replace by Bio.PDB equivalent
one_letter_residue_code = {
    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G',
    'HIS':'H', 'HID':'H', 'HIE':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
    'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S',
    'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
}

three_letter_residue_code = {
    'A':'ALA', 'C': 'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY',
    'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN',
    'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
    'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'
}

# Residue ids

def residue_id(r, models=False):
    return '{:>3} {}'.format(r.get_resname(), residuenum(r,models))

def residue_num (r, models=False):
    rn = str(r.get_parent().id) + str(r.id[1])
    if models:
        rn += "/" + str(r.get_parent().get_parent().id)
    return rn

def atom_id(at, models=False):
    return '{}.{}'.format(residueid(at.get_parent(), models),at.id)

# Id Checks
def residue_check(r):
    r = r.upper()
    rid = ''
    if r in three_letter_residue_code.keys():
        rid = three_letter_residue_code[r]
    elif r in one_letter_residue_code.keys():
        rid = r
    else:
        print ('#ERROR: unknown residue id ' + r)
        sys.exit(1)

    return rid

def same_residue (at1, at2):
    return at1.get_parent() == at2.get_parent()

def same_model(r1, r2):
    return r1.get_parent().get_parent() == r2.get_parent().get_parent()

def same_chain(r1, r2):
    return r1.get_parent() == r2.get_parent() and same_model(r1, r2)

def seq_consecutive(r1, r2):
    resnum1 = r1.id[1]
    resnum2 = r2.id[1]
    return same_chain(r1, r2) and abs(resnum1-resnum2) == 1

def is_wat(r):
    return r.id[0] == 'W'

def is_hetatm(r):
    return re.match('H_', r.id[0]) or re.match('W', r.id[0])

def check_chiral_residue(r,chiral_data):
    """
    Checks proper chirality of side chain atoms as defined in 
    chiral_data
    """
    ok = True
    if r.get_resname() in chiral_data:
        atids = chiral_data[r.get_resname()]
        ok = check_chiral(r, 'CB', 'CA', atids[0], atids[1], -1.)
    return ok

def check_chiral_ca(r):
    """
    Checks proper (L) quirality of CB atom with respect to backbone
    """
    return r.get_resname() == 'GLY' or check_chiral(r, 'N', 'CA','C','CB')

def check_chiral(r,at1,at2,at3,at4, sign=1.):
    """
    Checks proper chirality. 
    at1-at3 define reference plane. 
    Position of at4 with respect to the plane is checked. 
    Sign (+1,-1) allows to check for a specific enantiomer.
    """
    # check all atoms are present
    at_ok=True
    chi_ok=True
    for at in [at1,at2,at3,at4]:
        at_ok = at_ok and at in r
        if not at_ok:
            print ('Warning: atom {} not found in {}'.format(at,residueid(r)))
    if at_ok:
        v1=r[at1].coord-r[at2].coord
        v2=r[at3].coord-r[at2].coord
        vp = np.cross(v1,v2)
        v3=r[at4].coord-r[at2].coord
        chi_ok = sign * (_calc_v_angle(vp,v3) - 90.) < 0.
    return chi_ok

def get_altloc_residues(st):
    res_list = {}
    for at in st.get_atoms():
        r = at.get_parent()
        rid = residue_id(r)
        if at.get_altloc() != ' ':
            if rid not in res_list:
                res_list[rid] = []
            res_list[rid].append(at)
    return res_list

def get_metal_atoms(self, metal_ats):
    met_list = []
    for at in self.st.get_atoms():
        if not re.match('H_', at.get_parent().id[0]):
            continue
        if at.id in metal_ats:
            met_list.append(at)
    return met_list

def get_ligands(st, incl_water=False):
    lig_list = []
    for r in self.st.get_residues():
        if re.match('H_', r.id[0]) or (incl_water and re.match('W', r.id[0])):
            lig_list.append(r)
    return lig_list

def get_residues_with_H(st):
    resh_list = []
    for r in st.get_residues():
        has_h=0
        for a in r.get_atoms():
            if a.element == 'H':
                has_h +=1
        if has_h:
            resh_list.append({'r':r, 'n_h':has_h})
    return resh_list

# Residue manipulation =======================================================
def remove_H_from_r (r, verbose=False):
    H_list = []
    for at in r.get_atoms():
        if at.element == 'H':
            H_list.append(at.id)
    for at_id in H_list:
        if verbose:
            print ("  Deleting atom " + at_id)
        r.detach_child(at_id)

def remove_residue(r):
    r.get_parent().detach_child(r.id)
        
def swap_atom_names(at1,at2):
    at1_id = at1.id
    at1_full_id = at1.full_id
    at1_element = at1.element
    at1_name = at1.name
    at1_fullname = at1.fullname

    at1.id = at2.id
    at1.full_id = at2.full_id
    at1.element = at2.element
    at1.name = at2.name
    at1.fullname = at2.fullname

    at2.id = at1_id
    at2.full_id = at1_full_id
    at2.element = at1_element
    at2.name = at1_name
    at2.fullname = at1_fullname

def invert_side_atoms(r, res_data):
    res_type=r.get_resname()
    if not res_type in res_data:
        sys.stderr.write('Error: {} is not a valid residue'.format(res_type))
    util.swap_atom_names(r[res_data[res_type][0]],r[res_data[res_type][1]])

# Atom building ===============================================================
def buildCoordsOther(r, res_lib, new_res, at_id):

    resid_def = res_lib.residues[new_res]
    i = 1
    while resid_def.ats[i].id != at_id and i < len(resid_def.ats):
        i = i + 1
    if resid_def.ats[i].id == at_id:
        return buildCoords(
                           r[resid_def.ats[resid_def.ats[i].link[0]].id].get_coord(),
                           r[resid_def.ats[resid_def.ats[i].link[1]].id].get_coord(),
                           r[resid_def.ats[resid_def.ats[i].link[2]].id].get_coord(),
                           resid_def.ats[i].geom
                           )
    else:
        print ("#ERROR: Unknown target atom")
        sys.exit(1)

def buildCoordsCB(r): # Get CB from Backbone

    return buildCoords(
                       r['CA'].get_coord(),
                       r['N'].get_coord(),
                       r['C'].get_coord(),
                       [1.5, 115.5, -123.]
                       )

def buildCoords(avec, bvec, cvec, geom):

    dst = geom[0]
    ang = geom[1] * pi / 180.
    tor = geom[2] * pi / 180.0

    v1 = avec-bvec
    v2 = avec-cvec

    n = np.cross(v1, v2)
    nn = np.cross(v1, n)

    n /= norm(n)
    nn /= norm(nn)

    n *= -sin(tor)
    nn *= cos(tor)

    v3 = n + nn
    v3 /= norm(v3)
    v3 *= dst * sin(ang)

    v1 /= norm(v1)
    v1 *= dst * cos(ang)

    return avec + v3 - v1

# Metrics =============================================================
def calc_at_dist(at1,at2):
    return norm(at1.coord-at2.coord) #TODO look for a faster alternative

def calc_bond_angle(at1,at2, at3):
    v1 = at1.coord - at2.coord
    v2 = at3.coord - at2.coord
    return _calc_v_angle(v1,v2)

def calc_bond_dihedral(at1,at2,at3,at4):
    ab = at1.coord-at2.coord
    cb = at3.coord-at2.coord
    db = at4.coord-at3.coord
    u = np.cross(ab,cb)
    v = np.cross(db,cb)
    w = np.cross(u,v)
    angle_uv = _calc_v_angle(u,v)
    angle_cbw = _calc_v_angle(cb,w)
    try:
        if angle_cbw > 0.001:
            angle_uv = -angle_uv
    except ZeroDivisionError:
        pass
    return angle_uv

def get_all_at2at_distances(st, at_ids = ['all'], d_cutoff=0.):
    dist_mat = []
    at_list = []
    for at in st.get_atoms():
        if at.id in at_ids or at_ids == 'all':
            at_list.append(at)
    for i in range(0, len(at_list)-1):
        for j in range(i + 1, len(at_list)):
            d = calc_at_dist(at_list[i],at_list[j])
            if d_cutoff > 0. and d < d_cutoff:
                dist_mat.append ([at_list[i], at_list[j], d])
    return dist_mat

    
def get_all_r2r_distances(r_ids = ['all'], d_cutoff=0.):
    # Uses distances between the first atom of each residue as r-r distance
    dist_mat = []
    r_list = []
    for r in st.get_residues():
        if r.resname in r_ids or r_ids == 'all':
            r_list.append(r)
    for i in range(0, len(r_list)-1):
        ati = r_list[i].child_list[0]
        for j in range(i + 1, len(r_list)):
            atj = r_list[j].child_list[0]
            d = calc_at_dist(ati,atj)
            if d_cutoff > 0. and d < d_cutoff:
                dist_mat.append ([r_list[i], r_list[j], d])
    return dist_mat
        
def calc_RMSd_ats (ats1, ats2):
    if len(ats1) != len(ats2):
        print ("Warning: atom lists of different length when calculating RMSd", file=sys.stderr())
    
    rmsd = 0
    i = 0
    while i < len(ats1) and i < len(ats2):
        d = ats1[i] - ats2[i]
        rmsd = rmsd + d * d / len(ats1)
        i = i + 1

    return (math.sqrt(rmsd))

def calc_RMSd_all_ats (st1, st2):
    ats1 = []
    ats2 = []

    for at in st1.get_atoms():
        ats1.append(at)
    for at in st2.get_atoms():
        ats2.append(at)
        
    return calc_RMSd.ats(ats1,ats2)    

def get_all_rr_distances(r1, r2, with_h=False):
    dist_mat = []
    for at1 in r1.get_atoms():
        if at1.element == 'H' and not with_h:
            continue
        for at2 in r2.get_atoms():
            if at2.element == 'H' and not with_h:
                continue
            if at1.serial_number < at2.serial_number:
                dist_mat.append ([at1, at2, at1-at2])
    return dist_mat

def guess_models_type(st, threshold):
    if len(self.st)> 1:
        if util.calc_RMSd_all_ats(self.st[0], self.st[1]) < threshold:
            return NMR
        else:
            return TRAJ
    else:
        return UNKNOWN

#===============================================================================
def _calc_v_angle(v1,v2,deg=True): 
    angle = arccos(clip(dot(v1,v2)/norm(v1)/norm(v2),-1.,1.))
    if deg:
        angle *= 180./pi
    return angle
            