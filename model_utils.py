"""
  Utility functions to manipulate structures. based on Bio.PDB data model
"""

from Bio.PDB.Atom import Atom
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.vectors import Vector,rotaxis
import math
import numpy as np
from numpy import arccos
from numpy import clip
from numpy import cos
from numpy import dot
from numpy import pi
from numpy import sin
from numpy import sqrt
from numpy.linalg import norm
import re
import sys

UNKNOWN = 0

#chain types
PROTEIN = 1
DNA = 2
RNA = 3
NA = 4
SEQ_THRESHOLD = 0.8
chain_type_labels = {PROTEIN:'Protein', DNA:'DNA', RNA:'RNA', UNKNOWN:'Unknown'}

#Model Types
ENSM = 1
BUNIT = 2
MODELS_MAXRMS = 15.0    # Threshold value to detect NMR models (angs)
model_type_labels = {ENSM:'Ensembl/NMR', BUNIT:'BioUnit', UNKNOWN:'Unknown'}

#HetAtm Types
MODRES = 1
METAL = 2
ORGANIC = 3
COVORGANIC = 4
WAT = 5

#Geometry
CBINTERNALS = [1.5, 115.5, -123.]
OINTERNALS = [1.229, 120.500,0.000]
SP3ANGLE = 109.470
SP3DIHS = [60.0,180.0,300.0]


# TODO: consider replace by Bio.PDB equivalent
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

dna_residue_code = ['DA', 'DC', 'DG', 'DT']
rna_residue_code = ['A', 'C', 'G', 'U']
na_residue_code = dna_residue_code + rna_residue_code


# Residue & Atom friendly ids

def residue_id(r, models='auto'):
    """
    Friendly replacement for residue ids like ASN A324/0
    """
    return '{:>3} {}'.format(r.get_resname(), residue_num(r, models))

def residue_num (r, models='auto'):
    """
    Shortcut for getting residue num including chain id, includes model number if any
    """
    if models == 'auto':
        models = len(r.get_parent().get_parent().get_parent()) > 1

    if has_ins_code(r):
        rn = str(r.get_parent().id) + str(r.id[1]) + r.id[2]
    else:
        rn = str(r.get_parent().id) + str(r.id[1])
    if models:
        rn += "/" + str(r.get_parent().get_parent().id + 1)
    return rn

def atom_id(at, models='auto'):
    """
    Friendly replacement for atom ids like ASN A324/0.CA
    """
    return '{}.{}'.format(residue_id(at.get_parent(), models), at.id)

# Id Checks
def protein_residue_check(r):
    """
    Checks whether is a valid protein residue id, either one or three letter code
    return upper case three-letter code
    """
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
    """
    Checks whether atoms belong to the same residue
    """
    return at1.get_parent() == at2.get_parent()

def same_model(r1, r2):
    """
    Checks whether residues belong to the same model
    """
    return r1.get_parent().get_parent() == r2.get_parent().get_parent()

def same_chain(r1, r2):
    """
    Checks whether residues belong to the same chain
    """
    return (r1.get_parent() == r2.get_parent()) and same_model(r1, r2)

def seq_consecutive(r1, r2):
    """
    Checks whether residues belong to the same chain and are consecutive in sequences,
    taken from residue number
    """
    resnum1 = r1.id[1]
    resnum2 = r2.id[1]
    return same_chain(r1, r2) and abs(resnum1-resnum2) == 1

def is_wat(r):
    """
    Shortcut to check for water residues
    """
    return r.id[0] == 'W'

def is_hetatm(r):
    """
    Shortcut to check for HETATM residues
    """
    return r.id[0][0:2] == 'H_' or 'W' == r.id[0]

def is_at_in_list(at, at_list, rname=None):
    """Checks if **at is in **at_list in the context of **rname residue
    """
    if rname is None:
        rname = at.get_parent().get_resname().replace(' ', '')
    if rname in at_list:
        return at.id in at_list[rname] or at.id in at_list['*']
    else:
        return at.id in at_list['*']

def has_ins_code(r):
    """Checks whether residue **r has insertion code
    """
    return r.id[2] != ' '

#===============================================================================
# Guesses

def guess_models_type(st, threshold=MODELS_MAXRMS):
    """Guesses the type of models according to all atom RMSd of first 2 models
       Considers Ensembl/NMR if rmsd value is less than **threshold
    """
    if len(st) == 1:
        return 0
    rmsd = calc_RMSd_all_ats(st[0], st[1])
    rmsd = round(rmsd,4)
    if len(st) > 1:
        if rmsd < threshold:
            return {'type':ENSM, 'rmsd':rmsd}
        else:
            return {'type':BUNIT, 'rmsd':rmsd}
    else:
        return {'type':UNKNOWN, 'rmsd':rmsd}

def guess_chain_type(ch, thres=SEQ_THRESHOLD):
    """
    Guesses chain type (protein, dna, or rna) from residue composition
    Allow for non-std residues.
    """
    #TODO improve guessing for hybrid chains
    prot = 0.
    dna = 0.
    rna = 0.
    total = 0.
    for r in ch.get_residues():
        total += 1
        if is_wat(r):
            continue
        rname = r.get_resname().replace(' ', '')
        if rname in three_letter_residue_code.values():
            prot += 1
        elif rname in dna_residue_code:
            dna += 1
        elif rname in rna_residue_code:
            rna += 1
    prot = prot / total
    dna = dna / total
    rna = rna / total
    other = 1 - prot - dna - rna
    if prot > thres or prot > dna + rna:
        return PROTEIN
    elif dna > thres or dna > prot + rna:
        return DNA
    elif rna > thres or rna > prot + dna:
        return RNA
    else:
        return [prot, dna, rna, other]

#===============================================================================

def check_chiral_residue(r, chiral_data):
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
    return r.get_resname() == 'GLY' or check_chiral(r, 'N', 'CA', 'C', 'CB')

def check_chiral(r, at1, at2, at3, at4, sign=1.):
    """
    Checks proper chirality of at2 atom.
    at1-at3 define reference plane.
    Position of at4 with respect to the plane is checked.
    Sign (+1,-1) allows to check for a specific enantiomer.
    """
    # check all atoms are present
    at_ok = True
    chi_ok = True
    for at in [at1, at2, at3, at4]:
        at_ok = at_ok and at in r
        if not at_ok:
            print ('Warning: atom {:3} not found in {}'.format(at, residue_id(r)),file=sys.stderr)
    if at_ok:
        v1 = r[at1].coord-r[at2].coord
        v2 = r[at3].coord-r[at2].coord
        vp = np.cross(v1, v2)
        v3 = r[at4].coord-r[at2].coord
        chi_ok = sign * (_calc_v_angle(vp, v3) - 90.) < 0.
    return chi_ok

def check_all_at_in_r(r, at_list):
    miss_at = {}
    for group in ['backbone', 'side']:
        miss_at[group] = []
        for at_id in at_list[group]:
            if not at_id in r:
                miss_at[group].append(at_id)
    if len(miss_at['backbone'] + miss_at['side']) > 0:
        return miss_at
    else:
        return {}

def check_r_list_clashes(r_list, rr_list, CLASH_DIST, atom_lists, join_models=True):
    clash_list = {'severe':{}}
    for cls in atom_lists:
        clash_list[cls] = {}
    for r_pair in rr_list:
        [r1, r2, d] = r_pair

        if (r1 in r_list or r2 in r_list) and not is_wat(r1) and not is_wat(r2):
            c_list = check_rr_clashes(r1, r2, CLASH_DIST, atom_lists, join_models)
            rkey = residue_id(r1) + '-' + residue_id(r2)
            for cls in c_list:
                if len(c_list[cls]):
                    clash_list[cls][rkey] = c_list[cls]
    return clash_list

def check_rr_clashes(r1, r2, CLASH_DIST, atom_lists, join_models=True):
    clash_list = {}
    min_dist2 = {}
    CLASH_DIST2 = {}
    ats_list1 = {}
    ats_list2 = {}
    for cls in atom_lists:
        clash_list[cls] = []
        min_dist2[cls] = 99999.
        CLASH_DIST2[cls] = CLASH_DIST[cls]**2
        ats_list1[cls]=set()
        ats_list2[cls]=set()
    for at in r1.get_atoms():
        for cls in atom_lists:
            if is_at_in_list(at,atom_lists[cls], r1.get_resname()):
                ats_list1[cls].add(at.id)
    for at in r2.get_atoms():
        for cls in atom_lists:
            if is_at_in_list(at,atom_lists[cls], r2.get_resname()):
                ats_list2[cls].add(at.id)

    if r1 != r2 and not seq_consecutive(r1, r2) \
                and (join_models or same_model(r1, r2)):
        for at_pair in get_all_rr_distances(r1, r2):
            [at1, at2, dist2] = at_pair
            if 'severe' in atom_lists and dist2 < CLASH_DIST2['severe']:
                if dist2 < min_dist2:
                    clash_list['severe'] = at_pair
                    min_dist2['severe'] = dist2
            else:
                for cls in atom_lists:
                    if cls == 'apolar':
                        #Only one of the atoms should be apolar
                        if not at1.id in ats_list1[cls] and \
                           not at2.id in ats_list2[cls]:
                            continue
                        #Remove n->n+2 backbone clashes. TODO Improve
                        if abs(r1.index - r2.index) <= 2:
                            continue
                        #Remove Ca2+ looking like backbone CA's
                        if at1.id == 'CA' and is_hetatm(r1) or \
                           at2.id == 'CA' and is_hetatm(r2):
                           continue
                    else:
                        # Both atoms should be of the same kind
                        if not at1.id in ats_list1[cls] or \
                           not at2.id in ats_list2[cls]:
                            continue
                    if dist2 < CLASH_DIST2[cls]:
                        if dist2 < min_dist2[cls]:
                            clash_list[cls] = at_pair
                            min_dist2[cls] = dist2
    return clash_list

#===============================================================================
def get_altloc_residues(st):
    """ Gets list of residue  with atoms with alternative location labels
    """
    res_list = {}
    for r in st.get_residues():
        for at in r.get_atoms():
#        r = at.get_parent()
           if at.get_altloc() != ' ':
                if r not in res_list:
                    res_list[r] = []
                res_list[r].append(at)
    return res_list

def get_metal_atoms(st, metal_ats):
    """
    Gets list of metal atoms

    """
    met_list = []
    for at in st.get_atoms():
        if not re.match('H_', at.get_parent().id[0]):
            continue
        if at.id in metal_ats:
            met_list.append(at)
    return met_list

def get_ligands(st, incl_water=False):
    """
    Gets lists of ligands, water molecules can be excluded
    """
    lig_list = []
    for r in st.get_residues():
        if re.match('H_', r.id[0]) or (incl_water and re.match('W', r.id[0])):
            lig_list.append(r)
    return lig_list

def get_residues_with_H(st):
    """
    Get residues containing Hydrogen atoms

    """
    resh_list = []
    for r in st.get_residues():
        has_h = 0
        for a in r.get_atoms():
            if a.element == 'H':
                has_h += 1
        if has_h:
            resh_list.append({'r':r, 'n_h':has_h})
    return resh_list


def get_backbone_links(st, backbone_atoms, COVLNK, join_models=True): #TODO differenciate Protein and NA
    cov_links = []
    for m in st:
        bckats = []
        for at in st[m.id].get_atoms():
            if at.id in backbone_atoms:
                if at.disordered_flag:
                    bckats.append(at.selected_child)
                else:
                    bckats.append(at)
        if bckats:
            nbsearch = NeighborSearch(bckats)
            for at1, at2 in nbsearch.search_all(COVLNK):
                if not same_residue(at1, at2) \
                   and (join_models or same_model(at1.get_parent(),at2.get_parent())):
                    cov_links.append(sorted([at1, at2], key=lambda x: x.serial_number))
        else:
            print ("Warning: No backbone atoms defined", file=sys.stderr)

    return cov_links

# Residue manipulation =======================================================
def remove_H_from_r (r, verbose=False):
    """
    Removes Hydrogen atoms from given residue
    """
    H_list = []
    for at in r.get_atoms():
        if at.element == 'H':
            H_list.append(at.id)
    for at_id in H_list:
        if verbose:
            print ("  Deleting atom " + at_id)
        r.detach_child(at_id)

def remove_residue(r):
    """
    Removes residue completely
    """
    r.get_parent().detach_child(r.id)

def swap_atoms(at1, at2):
    """
    Swaps names for two given atoms. Useful to fix labelling issues
    """
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

def invert_side_atoms(r, res_data): #TODO check merging with swat_atom_names
    """
     Swaps atoms according to res_data. Useful for fixing amides and chirality issues
    """
    #TODO reconstruct CD1 atom in Ile
    res_type = r.get_resname()
    if not res_type in res_data:
        sys.stderr.write('Error: {} is not a valid residue'.format(res_type))
    swap_atoms(r[res_data[res_type][0]], r[res_data[res_type][1]])

def invert_chirality (r, at1, at2, at3, at4):
    """Inverts chirality of at2 by rotating at4, and the associated end chain atoms
    """
    #TODO
    pass

def invert_chiral_ca(r):
    """
     Inverts CA Chirality.
    """
    #TODO
    return None

# Atom management ==============================================================
def rename_atom(r, old_at, new_at):
    at = r[old_at]
    r.detach_child(at.id)
    at.id = new_at
    at.full_id = new_at
    at.element = new_at[0:1]
    at.fullname = ' ' + new_at
    r.add(at)

def delete_atom(r, at_id):
    r.detach_child(at_id)


def build_atom(r, at_id, res_lib, new_res_id):
    if at_id == 'CB':
        coords = build_coords_CB(r)
    else:
        coords = build_coords_from_lib(r, res_lib, new_res_id, at_id)
    add_new_atom_to_residue(r,at_id,coords)

def add_new_atom_to_residue(r, at_id, coords):
    r.add(Atom(at_id,coords,99.0,1.0,' ',' ' + at_id + ' ',0,at_id[0:1]))

def build_coords_from_lib(r, res_lib, new_res, at_id):
    """
     Calculates cartesian coordinates for a new atom from internal coordinates definition.

    """
    atom_def = res_lib.get_atom_def(new_res, at_id)
    if atom_def == None:
        print ("#ERROR: Unknown target atom")
        sys.exit(1)

    return build_coords_from_ats_internal(
            r[atom_def.link_ats[0]],
            r[atom_def.link_ats[1]],
            r[atom_def.link_ats[2]],
             atom_def.geom
        )


def build_coords_CB(r): # Get CB from Backbone
    """
     Calculates cartesian coordinates for a new CB atom from backbone.

    """
    return build_coords_from_ats_internal(r['CA'],r['N'],r['C'],CBINTERNALS)

def build_coords_O(r): # Get O from Backbone
    """
     Calculates cartesian coordinates for a new O atom from backbone.

    """
    return build_coords_from_ats_internal(r['C'],r['CA'],r['N'],OINTERNALS
                       )

def build_coords_3xSP3(dst,at,at1,at2):
    """
        Generates coordinates for 3 SP3 atoms
        **dst** bond distance
        **at**  central atom
        **at1** atom to define bond angles
        **at2** atom to define dihedrals
    """
    #TODO try a pure geometrical generation to avoid at2
    crs=[]
    for i in range(0,3):
        crs.append(build_coords_from_ats_internal(at,at1,at2,[dst,SP3ANGLE,SP3DIHS[i]]))
    return crs


def build_coords_2xSP3(dst,at,at1,at2):
    """
        Generates coordinates for two SP3 bonds given the other two
        **dst** Bond distance
        **at** Central atom
        **at1** atom with existing bond
        **at2** atom with existing bond
    """
    cr0 = Vector(at.get_coord())
    cr1 = Vector(at1.get_coord())
    cr2 = Vector(at2.get_coord())
    axe = cr0-cr1
    m = rotaxis(120.*pi/180., axe)
    bond = cr2-cr0
    bond.normalize()
    bond._ar= bond._ar * dst
    cr3=cr0+bond.left_multiply(m)
    cr4=cr0+bond.left_multiply(m).left_multiply(m)
    crs=[]
    crs.append(cr3._ar)
    crs.append(cr4._ar)
    return crs

def build_coords_1xSP3(dst, at, at1, at2, at3):
    """
      Calculated cartesian coordinates to complete a SP3 group
    """
    cr0=at.get_coord()
    cr1=at1.get_coord()
    cr2=at2.get_coord()
    cr3=at3.get_coord()
    avg = cr1 + cr2
    avg = avg + cr3
    avg /= 3.
    avec = cr0 - avg
    avec /= norm(avec)
    avec *= dst
    return cr0 + avec

def build_coords_SP2(dst, at, at1, at2):
    """
      Calculates cartesian coordinaties to complete a SP2 group
    """
    cr0 = at.get_coord()
    cr1 = at1.get_coord()
    cr2 = at2.get_coord()

    avg = cr1 + cr2
    avg /= 2.
    avec = cr0 - avg
    avec /= norm(avec)
    avec *= dst
    return cr0 + avec


def build_coords_from_ats_internal(at1, at2, at3, geom):
    return build_coords_from_internal(
        at1.get_coord(),
        at2.get_coord(),
        at3.get_coord(),
        geom)

def build_coords_from_internal(at1c, at2c, at3c, geom):
    """
     Calculates cartesian coordinates for a new atom from internal coordinates.

    """
    dst = geom[0]
    ang = geom[1] * pi / 180.
    tor = geom[2] * pi / 180.0

    v1 = at1c-at2c
    v2 = at1c-at3c

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

    return at1c + v3 - v1

# Metrics =============================================================
def calc_at_dist(at1, at2):
    """
    Calculates distance between two atoms
    """
    return np.sqrt(calc_at_sq_dist(at1, at2))

def calc_at_sq_dist(at1, at2):
    """
    Calculates distance between two atoms
    """
    v = at1.coord-at2.coord
    return np.dot(v, v)

def calc_bond_angle(at1, at2, at3):
    """
    Calculates angle among three atoms at1-at2-at3
    """
    v1 = at1.coord - at2.coord
    v2 = at3.coord - at2.coord
    return _calc_v_angle(v1, v2)

def calc_bond_dihedral(at1, at2, at3, at4):
    """
    Calculates dihedral angles at1-at2-at3-at4
    """
    ab = at1.coord-at2.coord
    cb = at3.coord-at2.coord
    db = at4.coord-at3.coord
    u = np.cross(ab, cb)
    v = np.cross(db, cb)
    w = np.cross(u, v)
    angle_uv = _calc_v_angle(u, v)
    if norm(w)== 0.:
        angle_cbw = 0.
    else:
        angle_cbw = _calc_v_angle(cb, w)
    try:
        if angle_cbw > 0.001:
            angle_uv = -angle_uv
    except ZeroDivisionError:
        pass
    return angle_uv

def get_all_at2at_distances(st, at_ids='all', d_cutoff=0., join_models=False):
    """
    Gets a list of all at-at distances below a cutoff, at ids can be limited
    """
    if not isinstance(at_ids, list):
        at_ids = at_ids.split(',')

    dist_mat = []
    at_list = []
    d_cut2 = d_cutoff**2
    for at in st.get_atoms():
        if at.id in at_ids or at_ids == ['all']:
            at_list.append(at)
    for i in range(len(at_list)-1):
        for j in range(i + 1, len(at_list)):
            if join_models or same_model(at_list[i].get_parent(), at_list[j].get_parent()):
                d = calc_at_sq_dist(at_list[i], at_list[j])
                if d_cutoff > 0. and d < d_cut2:
                    dist_mat.append ([at_list[i], at_list[j], d])
    return dist_mat


def get_all_r2r_distances(st, r_ids='all', d_cutoff=0., join_models=False):
    # Uses distances between the first atom of each residue as r-r distance
    if not isinstance(r_ids, list):
        r_ids = r_ids.split(',')
    dist_mat = []
    check_ats={}
    for md_id in range(len(st)):
        check_ats[md_id]=[]
        for r in st[md_id].get_residues():
            if r.resname in r_ids or r_ids == ['all']:
                if join_models:
                    check_ats[0].append(r.child_list[0])
                else:
                    check_ats[md_id].append(r.child_list[0])
    for md_id in range(len(st)):
        if len(check_ats[md_id]):
            dist_mat += _get_contacts(check_ats[md_id],d_cutoff)
    return dist_mat


def _get_contacts(ats_list,d_cutoff):
    contact_list=[]
    nbsearch = NeighborSearch(ats_list)
    for at1, at2 in nbsearch.search_all(d_cutoff):
        contact_list.append([at1.get_parent(), at2.get_parent(), at1-at2])
    return contact_list


def calc_RMSd_ats (ats1, ats2):
    if len(ats1) != len(ats2):
        print ("Warning: atom lists of different length when calculating RMSd ({}, {})".format(len(ats1),len(ats2)), file=sys.stderr)
    rmsd = 0
    i = 0
    while i < len(ats1) and i < len(ats2):
        d2 = calc_at_sq_dist(ats1[i], ats2[i])
        rmsd = rmsd + d2
        i = i + 1

    return (sqrt(rmsd / i))

def calc_RMSd_all_ats (st1, st2):
    ats1 = []
    ats2 = []

    for at in st1.get_atoms():
        ats1.append(at)
    for at in st2.get_atoms():
        ats2.append(at)

    return calc_RMSd_ats(ats1, ats2)

def get_all_rr_distances(r1, r2, with_h=False):
    dist_mat = []
    for at1 in r1.get_atoms():
        if at1.element == 'H' and not with_h:
            continue
        for at2 in r2.get_atoms():
            if at2.element == 'H' and not with_h:
                continue
            d2 = calc_at_sq_dist(at1,at2)
#            if at1.serial_number < at2.serial_number:
            dist_mat.append ([at1, at2, d2])
#            else:
#                dist_mat.append ([at2, at1, d2])
    return dist_mat


#===============================================================================
def _calc_v_angle(v1, v2, deg=True):
    angle = arccos(clip(dot(v1, v2) / norm(v1)/norm(v2),-1.,1.))
    if deg:
        angle *= 180./pi
    return angle

