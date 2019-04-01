"""
  Utility functions to manipulate structures. based on Bio.PDB data model
"""

import re
import sys
import numpy as np
from numpy import arccos, clip, cos, dot, pi, sin, sqrt
from numpy.linalg import norm

from Bio.PDB.Atom import Atom
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.vectors import Vector, rotaxis

UNKNOWN = 0

#chain types
PROTEIN = 1
DNA = 2
RNA = 3
NA = 4
SEQ_THRESHOLD = 0.8
CHAIN_TYPE_LABELS = {PROTEIN:'Protein', DNA:'DNA', RNA:'RNA', UNKNOWN:'Unknown'}

#Model Types
ENSM = 1
BUNIT = 2
MODELS_MAXRMS = 15.0    # Threshold value to detect NMR models (angs)
MODEL_TYPE_LABELS = {ENSM:'Ensembl/NMR', BUNIT:'BioUnit', UNKNOWN:'Unknown'}

#HetAtm Types
MODRES = 1
METAL = 2
ORGANIC = 3
COVORGANIC = 4
WAT = 5

#Geometry
CBINTERNALS = [1.5, 115.5, -123.]
OINTERNALS = [1.229, 120.500, 0.000]
SP3ANGLE = 109.470
SP3DIHS = [60.0, 180.0, 300.0]
HDIS = 1.08


# TODO: consider replace by Bio.PDB equivalent
ONE_LETTER_RESIDUE_CODE = {
    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G',
    'HIS':'H', 'HID':'H', 'HIE':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
    'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S',
    'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
}

THREE_LETTER_RESIDUE_CODE = {
    'A':'ALA', 'C': 'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY',
    'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN',
    'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
    'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'
}

DNA_RESIDUE_CODE = ['DA', 'DC', 'DG', 'DT']
RNA_RESIDUE_CODE = ['A', 'C', 'G', 'U']
NA_RESIDUE_CODE = DNA_RESIDUE_CODE + RNA_RESIDUE_CODE

# Residue & Atom friendly ids
def residue_id(res, models='auto'):
    """
    Friendly replacement for residue ids like ASN A324/0
    """
    return '{:>3} {}'.format(res.get_resname(), residue_num(res, models))

def residue_num(res, models='auto'):
    """
    Shortcut for getting residue num including chain id, includes model number if any
    """
    if models == 'auto':
        models = len(res.get_parent().get_parent().get_parent()) > 1

    if has_ins_code(res):
        rnum = str(res.get_parent().id) + str(res.id[1]) + res.id[2]
    else:
        rnum = str(res.get_parent().id) + str(res.id[1])
    if models:
        rnum += "/" + str(res.get_parent().get_parent().id + 1)
    return rnum

def atom_id(atm, models='auto'):
    """
    Friendly replacement for atom ids like ASN A324/0.CA
    """
    return '{}.{}'.format(residue_id(atm.get_parent(), models), atm.id)

# Id Checks
def protein_residue_check(res):
    """
    Checks whether is a valid protein residue id, either one or three letter code
    return upper case three-letter code
    """
    res = res.upper()
    rid = ''
    if res in THREE_LETTER_RESIDUE_CODE.keys():
        rid = THREE_LETTER_RESIDUE_CODE[res]
    elif res in ONE_LETTER_RESIDUE_CODE.keys():
        rid = res
    else:
        return False

    return rid

def same_residue(at1, at2):
    """
    Checks whether atoms belong to the same residue
    """
    return at1.get_parent() == at2.get_parent()

def same_model(res1, res2):
    """
    Checks whether residues belong to the same model
    """
    return res1.get_parent().get_parent() == res2.get_parent().get_parent()

def same_chain(res1, res2):
    """
    Checks whether residues belong to the same chain
    """
    return (res1.get_parent() == res2.get_parent()) and same_model(res1, res2)

def seq_consecutive(res1, res2):
    """
    Checks whether residues belong to the same chain and are consecutive in sequences,
    taken from residue number
    """
    rnum1 = res1.id[1]
    rnum2 = res2.id[1]
    return same_chain(res1, res2) and abs(rnum1 - rnum2) == 1

def is_wat(res):
    """
    Shortcut to check for water residues
    """
    return res.id[0].startswith('W')

def is_hetatm(res):
    """
    Shortcut to check for HETATM residues
    """
    return res.id[0].startswith('H_') or res.id[0].startswith('W')

def is_at_in_list(atm, at_list, rname=None):
    """Checks if **at is in **at_list in the context of **rname residue
    """
    if rname is None:
        rname = atm.get_parent().get_resname().replace(' ', '')
    if rname in at_list:
        return atm.id in at_list[rname] or atm.id in at_list['*']
    return atm.id in at_list['*']

def has_ins_code(res):
    """Checks whether residue **r has insertion code
    """
    return res.id[2] != ' '

#===============================================================================
# Guesses

def guess_models_type(struc, threshold=MODELS_MAXRMS):
    """Guesses the type of models according to all atom RMSd of first 2 models
       Considers Ensembl/NMR if rmsd value is less than **threshold
    """
    if len(struc) == 1:
        return 0
    rmsd = calc_RMSd_all_ats(struc[0], struc[1])
    rmsd = round(rmsd, 4)
    if rmsd < threshold:
        return {'type':ENSM, 'rmsd':rmsd}
    return {'type':BUNIT, 'rmsd':rmsd}

def guess_chain_type(chn, thres=SEQ_THRESHOLD):
    """
    Guesses chain type (protein, dna, or rna) from residue composition
    Allow for non-std residues.
    """
    #TODO improve guessing for hybrid chains
    prot = 0.
    dna = 0.
    rna = 0.
    total = 0.
    for res in chn.get_residues():
        total += 1
        if is_wat(res):
            continue
        rname = res.get_resname().replace(' ', '')
        if rname in THREE_LETTER_RESIDUE_CODE.values():
            prot += 1
        elif rname in DNA_RESIDUE_CODE:
            dna += 1
        elif rname in RNA_RESIDUE_CODE:
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
    return [prot, dna, rna, other]

#===============================================================================
def check_chiral_residue(res, chiral_data):
    """
        Checks proper chirality of side chain atoms as defined in
        chiral_data
    """
    chi_ok = True
    if res.get_resname() in chiral_data:
        at_ids = chiral_data[res.get_resname()]
        chi_ok = check_chiral(res, 'CB', 'CA', at_ids[0], at_ids[1], -1.)
    return chi_ok

def check_chiral_ca(res):
    """
    Checks proper (L) quirality of CB atom with respect to backbone
    """
    return res.get_resname() == 'GLY' or check_chiral(res, 'N', 'CA', 'C', 'CB')

def check_chiral(res, at1, at2, at3, at4, sign=1.):
    """
    Checks proper chirality of at2 atom.
    at1-at3 define reference plane.
    Position of at4 with respect to the plane is checked.
    Sign (+1,-1) allows to check for a specific enantiomer.
    """
    # check all atoms are present
    at_ok = True
    chi_ok = True
    for atm in [at1, at2, at3, at4]:
        at_ok = at_ok and atm in res
        if not at_ok:
            print('Warning: atom {:3} not found in {}'.format(atm, residue_id(res)))
    if at_ok:
        vec1 = res[at1].coord - res[at2].coord
        vec2 = res[at3].coord - res[at2].coord
        vecp = np.cross(vec1, vec2)
        vec3 = res[at4].coord - res[at2].coord
        chi_ok = sign * (_calc_v_angle(vecp, vec3) - 90.) < 0.
    return chi_ok

def check_all_at_in_r(res, at_list):
    """ Check whether all residue atoms are present. """
    miss_at = {}
    for group in ['backbone', 'side']:
        miss_at[group] = []
        for at_id in at_list[group]:
            if not at_id in res:
                miss_at[group].append(at_id)
    if miss_at['backbone'] + miss_at['side']:
        return miss_at
    return {}

def check_r_list_clashes(r_list, rr_list, clash_dist, atom_lists, join_models=True):
    """ Check clashes generated by a list of residues """
    clash_list = {'severe':{}}
    for cls in atom_lists:
        clash_list[cls] = {}
    for r_pair in rr_list:
        [res1, res2] = r_pair[0:2]

        if (res1 in r_list or res2 in r_list) and not is_wat(res1) and not is_wat(res2):
            c_list = check_rr_clashes(res1, res2, clash_dist, atom_lists, join_models)
            rkey = residue_id(res1) + '-' + residue_id(res2)
            for cls in c_list:
                if c_list[cls]:
                    clash_list[cls][rkey] = c_list[cls]
    return clash_list

def check_rr_clashes(res1, res2, clash_dist, atom_lists, join_models=True):
    """ Check all clashes between two residues """
    clash_list = {}
    min_dist2 = {}
    clash_dist2 = {}
    ats_list1 = {}
    ats_list2 = {}
    for cls in atom_lists:
        clash_list[cls] = []
        min_dist2[cls] = 99999.
        clash_dist2[cls] = clash_dist[cls]**2
        ats_list1[cls] = set()
        ats_list2[cls] = set()
    for atm in res1.get_atoms():
        for cls in atom_lists:
            if is_at_in_list(atm, atom_lists[cls], res1.get_resname()):
                ats_list1[cls].add(atm.id)
    for atm in res2.get_atoms():
        for cls in atom_lists:
            if is_at_in_list(atm, atom_lists[cls], res2.get_resname()):
                ats_list2[cls].add(atm.id)

    if res1 != res2 and not seq_consecutive(res1, res2) \
                and (join_models or same_model(res1, res2)):
        for at_pair in get_all_rr_distances(res1, res2):
            [at1, at2, dist2] = at_pair
            if 'severe' in atom_lists and dist2 < clash_dist2['severe']:
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
                        if abs(res1.index - res2.index) <= 2:
                            continue
                        #Remove Ca2+ looking like backbone CA's
                        if at1.id == 'CA' and is_hetatm(res1) or \
                                at2.id == 'CA' and is_hetatm(res2):
                            continue
                    else:
                        # Both atoms should be of the same kind
                        if not at1.id in ats_list1[cls] or \
                                not at2.id in ats_list2[cls]:
                            continue
                    if dist2 < clash_dist2[cls]:
                        if dist2 < min_dist2[cls]:
                            clash_list[cls] = at_pair
                            min_dist2[cls] = dist2
    return clash_list

#===============================================================================
def get_altloc_residues(struc):
    """ Gets list of residue  with atoms with alternative location labels
    """
    res_list = {}
    for res in struc.get_residues():
        for atm in res.get_atoms():
            if atm.get_altloc() != ' ':
                if res not in res_list:
                    res_list[res] = []
                res_list[res].append(atm)
    return res_list

def get_metal_atoms(struc, metal_ats):
    """
    Gets list of metal atoms

    """
    met_list = []
    for atm in struc.get_atoms():
        if not re.match('H_', atm.get_parent().id[0]):
            continue
        if atm.id in metal_ats:
            met_list.append(atm)
    return met_list

def get_ligands(struc, incl_water=False):
    """
    Gets lists of ligands, water molecules can be excluded
    """
    lig_list = []
    for res in struc.get_residues():
        if is_hetatm(res) and (incl_water or not is_wat(res)):
            lig_list.append(res)
    return lig_list

def get_residues_with_H(struc):
    """ Get residues containing Hydrogen atoms
    """
    resh_list = []
    for res in struc.get_residues():
        has_h = 0
        for atm in res.get_atoms():
            if atm.element == 'H':
                has_h += 1
        if has_h:
            resh_list.append({'r':res, 'n_h':has_h})
    return resh_list


def get_backbone_links(struc, backbone_atoms, covlnk, join_models=True):
    """ Get links making the main chain """
    #TODO differenciate Protein and NA
    cov_links = []
    for mod in struc:
        bckats = []
        for atm in struc[mod.id].get_atoms():
            if atm.id in backbone_atoms:
                if atm.disordered_flag:
                    bckats.append(atm.selected_child)
                else:
                    bckats.append(atm)
        if bckats:
            nbsearch = NeighborSearch(bckats)
            for at1, at2 in nbsearch.search_all(covlnk):
                if not same_residue(at1, at2) \
                        and (join_models or same_model(at1.get_parent(), at2.get_parent())):
                    cov_links.append(sorted([at1, at2], key=lambda x: x.serial_number))
        else:
            print("Warning: No backbone atoms defined")

    return cov_links

# Residue manipulation =======================================================
def remove_H_from_r(res, verbose=False):
    """
    Removes Hydrogen atoms from given residue
    """
    h_list = []
    for atm in res.get_atoms():
        if atm.element == 'H':
            h_list.append(atm.id)
    for at_id in h_list:
        if verbose:
            print("  Deleting atom " + at_id)
        res.detach_child(at_id)

def remove_residue(res):
    """
    Removes residue completely
    """
    res.get_parent().detach_child(res.id)

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

#def invert_chirality(res, at1, at2, at3, at4):
#    """Inverts chirality of at2 by rotating at4, and the associated end chain atoms
    #"""
#    #TODO
#
#def invert_chiral_ca(res):
#    """
#    Inverts CA Chirality.
    #"""
#    #TODO

# Atom management ==============================================================
def rename_atom(res, old_at, new_at):
    atm = res[old_at]
    res.detach_child(atm.id)
    atm.id = new_at
    atm.full_id = new_at
    atm.element = new_at[0:1]
    atm.fullname = ' ' + new_at
    res.add(atm)

def delete_atom(res, at_id):
    res.detach_child(at_id)

def add_hydrogens_backbone(res, res_1):
    """ Add hydrogen atoms to the backbone"""

    # only proteins

    if not protein_residue_check(res.get_resname()):
        return "Warning: Residue not valid in this context "

    error_msg = "Warning: not enough atoms to build backbone hydrogen atoms on"

    if 'N' not in res or 'CA' not in res:
        return error_msg

    if res_1 is None:
        # Nterminal TODO  Neutral NTerm
        if res.get_resname() == 'PRO':
            if 'CD' not in res:
                return error_msg

            add_new_atom_to_residue(
                res,
                'H',
                build_coords_SP2(HDIS, res['N'], res['CA'], res['CD'])
            )
        else:
            if 'C' not in res:
                return error_msg

            crs = build_coords_3xSP3(HDIS, res['N'], res['CA'], res['C'])
            add_new_atom_to_residue(res, 'H1', crs[0])
            add_new_atom_to_residue(res, 'H2', crs[1])
            add_new_atom_to_residue(res, 'H3', crs[2])

    elif res.get_resname() != 'PRO':
        if 'C' not in res_1:
            return error_msg

        add_new_atom_to_residue(
            res,
            'H',
            build_coords_SP2(HDIS, res['N'], res['CA'], res_1['C'])
        )

    if res.get_resname() == 'GLY':
        if 'C' not in res:
            return error_msg

        crs = build_coords_2xSP3(HDIS, res['CA'], res['N'], res['C'])
        add_new_atom_to_residue(res, 'HA2', crs[0])
        add_new_atom_to_residue(res, 'HA3', crs[1])

    else:
        if 'C' not in res or 'CB' not in res:
            return error_msg

        add_new_atom_to_residue(
            res,
            'HA',
            build_coords_1xSP3(HDIS, res['CA'], res['N'], res['C'], res['CB'])
        )

    return False



def add_hydrogens_side(res, res_library, opt, rules):
    """ Add hydrogens to side chains"""
    if 'N' not in res or 'CA' not in res or 'C' not in res:
        return "Warning: not enough atoms to build side chain hydrogen atoms on"

    for key_rule in rules.keys():
        rule = rules[key_rule]
        if rule['mode'] == 'B2':
            crs = build_coords_2xSP3(
                rule['dist'],
                res[key_rule],
                res[rule['ref_ats'][0]],
                res[rule['ref_ats'][1]]
            )
            add_new_atom_to_residue(res, rule['ats'][0], crs[0])
            add_new_atom_to_residue(res, rule['ats'][1], crs[1])

        elif rule['mode'] == "B1":
            crs = build_coords_1xSP3(
                rule['dist'],
                res[key_rule],
                res[rule['ref_ats'][0]],
                res[rule['ref_ats'][1]],
                res[rule['ref_ats'][2]]
            )
            add_new_atom_to_residue(res, rule['ats'][0], crs)

        elif rule['mode'] == 'S2':
            crs = build_coords_SP2(
                rule['dist'],
                res[key_rule],
                res[rule['ref_ats'][0]],
                res[rule['ref_ats'][1]],
            )
            add_new_atom_to_residue(res, rule['ats'][0], crs)

        elif rule['mode'] == 'L':
            for at_id in rule['ats']:
                crs = build_coords_from_lib(res, res_library, opt, at_id)
                add_new_atom_to_residue(res, at_id, crs)
    return False

def build_atom(res, at_id, res_lib, new_res_id):
    if at_id == 'CB':
        coords = build_coords_CB(res)
    else:
        coords = build_coords_from_lib(res, res_lib, new_res_id, at_id)
    add_new_atom_to_residue(res, at_id, coords)


def add_new_atom_to_residue(res, at_id, coords):
    res.add(Atom(at_id, coords, 99.0, 1.0, ' ', ' ' + at_id + ' ', 0, at_id[0:1]))


def build_coords_from_lib(res, res_lib, new_res, at_id):
    """
     Calculates cartesian coordinates for a new atom from internal coordinates definition.

    """
    atom_def = res_lib.get_atom_def(new_res, at_id)
    if atom_def is None:
        print("#ERROR: Unknown target atom")
        sys.exit(1)

    return build_coords_from_ats_internal(
        res[atom_def.link_ats[0]],
        res[atom_def.link_ats[1]],
        res[atom_def.link_ats[2]],
        atom_def.geom
    )

def build_coords_CB(res): # Get CB from Backbone
    """
     Calculates cartesian coordinates for a new CB atom from backbone.

    """
    return build_coords_from_ats_internal(res['CA'], res['N'], res['C'], CBINTERNALS)

def build_coords_O(res): # Get O from Backbone
    """
     Calculates cartesian coordinates for a new O atom from backbone.

    """
    return build_coords_from_ats_internal(res['C'], res['CA'], res['N'], OINTERNALS)

def build_coords_3xSP3(dst, at0, at1, at2):
    """
        Generates coordinates for 3 SP3 atoms
        **dst** bond distance
        **at0**  central atom
        **at1** atom to define bond angles
        **at2** atom to define dihedrals
    """
    #TODO try a pure geometrical generation to avoid at2
    crs = []
    for i in range(0, 3):
        crs.append(build_coords_from_ats_internal(at0, at1, at2, [dst, SP3ANGLE, SP3DIHS[i]]))
    return crs

def build_coords_2xSP3(dst, at0, at1, at2):
    """
        Generates coordinates for two SP3 bonds given the other two
        **dst** Bond distance
        **at0** Central atom
        **at1** atom with existing bond
        **at2** atom with existing bond
    """
    cr0 = Vector(at0.get_coord())
    cr1 = Vector(at1.get_coord())
    cr2 = Vector(at2.get_coord())
    axe = cr0 - cr1
    mat = rotaxis(120.*pi/180., axe)
    bond = cr2 - cr0
    bond.normalize()
    bond._ar = bond._ar * dst
    cr3 = cr0 + bond.left_multiply(mat)
    cr4 = cr0 + bond.left_multiply(mat).left_multiply(mat)
    crs = []
    crs.append(cr3._ar)
    crs.append(cr4._ar)
    return crs

def build_coords_1xSP3(dst, at0, at1, at2, at3):
    """
      Calculated cartesian coordinates to complete a SP3 group
    """
    cr0 = at0.get_coord()
    cr1 = at1.get_coord()
    cr2 = at2.get_coord()
    cr3 = at3.get_coord()
    avg = cr1 + cr2
    avg = avg + cr3
    avg /= 3.
    avec = cr0 - avg
    avec /= norm(avec)
    avec *= dst
    return cr0 + avec

def build_coords_SP2(dst, at0, at1, at2):
    """
      Calculates cartesian coordinaties to complete a SP2 group
    """
    cr0 = at0.get_coord()
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

    vec1 = at1c - at2c
    vec2 = at1c - at3c

    vcr12 = np.cross(vec1, vec2)
    vcr112 = np.cross(vec1, vcr12)

    vcr12 /= norm(vcr12)
    vcr112 /= norm(vcr112)

    vcr12 *= -sin(tor)
    vcr112 *= cos(tor)

    vec3 = vcr12 + vcr112
    vec3 /= norm(vec3)
    vec3 *= dst * sin(ang)

    vec1 /= norm(vec1)
    vec1 *= dst * cos(ang)

    return at1c + vec3 - vec1

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
    vec = at1.coord - at2.coord
    return np.dot(vec, vec)

def calc_bond_angle(at1, at2, at3):
    """
    Calculates angle among three atoms at1-at2-at3
    """
    vec1 = at1.coord - at2.coord
    vec2 = at3.coord - at2.coord
    return _calc_v_angle(vec1, vec2)

def calc_bond_dihedral(at1, at2, at3, at4):
    """
    Calculates dihedral angles at1-at2-at3-at4
    """
    abv = at1.coord - at2.coord
    cbv = at3.coord - at2.coord
    dbv = at4.coord - at3.coord
    uvec = np.cross(abv, cbv)
    vvec = np.cross(dbv, cbv)
    wvec = np.cross(uvec, vvec)
    angle_uv = _calc_v_angle(uvec, vvec)
    if norm(wvec) == 0.:
        angle_cbw = 0.
    else:
        angle_cbw = _calc_v_angle(cbv, wvec)
    try:
        if angle_cbw > 0.001:
            angle_uv = -angle_uv
    except ZeroDivisionError:
        pass
    return angle_uv

def get_all_at2at_distances(struc, at_ids='all', d_cutoff=0., join_models=False):
    """
    Gets a list of all at-at distances below a cutoff, at ids can be limited
    """
    if not isinstance(at_ids, list):
        at_ids = at_ids.split(',')

    dist_mat = []
    at_list = []
    d_cut2 = d_cutoff**2
    for atm in struc.get_atoms():
        if atm.id in at_ids or at_ids == ['all']:
            at_list.append(atm)
    for i in range(len(at_list)-1):
        for j in range(i + 1, len(at_list)):
            if join_models or same_model(at_list[i].get_parent(), at_list[j].get_parent()):
                dist2 = calc_at_sq_dist(at_list[i], at_list[j])
                if d_cutoff > 0. and dist2 < d_cut2:
                    dist_mat.append([at_list[i], at_list[j], dist2])
    return dist_mat


def get_all_r2r_distances(struc, r_ids='all', d_cutoff=0., join_models=False):
    # Uses distances between the first atom of each residue as r-r distance
    if not isinstance(r_ids, list):
        r_ids = r_ids.split(',')
    dist_mat = []
    check_ats = {}
    for md_id in range(len(struc)):
        check_ats[md_id] = []
        for res in struc[md_id].get_residues():
            if res.resname in r_ids or r_ids == ['all']:
                if join_models:
                    check_ats[0].append(res.child_list[0])
                else:
                    check_ats[md_id].append(res.child_list[0])
    for md_id in range(len(struc)):
        if check_ats[md_id]:
            dist_mat += _get_contacts(check_ats[md_id], d_cutoff)
    return dist_mat

def _get_contacts(ats_list, d_cutoff):
    contact_list = []
    nbsearch = NeighborSearch(ats_list)
    for at1, at2 in nbsearch.search_all(d_cutoff):
        contact_list.append([at1.get_parent(), at2.get_parent(), at1 - at2])
    return contact_list

def calc_RMSd_ats(ats1, ats2):
    if len(ats1) != len(ats2):
        print(
            "Warning: atom lists of different length when calculating RMSd ({}, {})".format(
                len(ats1), len(ats2)
            )
        )
    rmsd = 0
    i = 0
    while i < len(ats1) and i < len(ats2):
        dist2 = calc_at_sq_dist(ats1[i], ats2[i])
        rmsd = rmsd + dist2
        i = i + 1

    return sqrt(rmsd / i)

def calc_RMSd_all_ats(st1, st2):
    ats1 = []
    ats2 = []

    for atm in st1.get_atoms():
        ats1.append(atm)
    for atm in st2.get_atoms():
        ats2.append(atm)

    return calc_RMSd_ats(ats1, ats2)

def get_all_rr_distances(res1, res2, with_h=False):
    dist_mat = []
    for at1 in res1.get_atoms():
        if at1.element == 'H' and not with_h:
            continue
        for at2 in res2.get_atoms():
            if at2.element == 'H' and not with_h:
                continue
            dist2 = calc_at_sq_dist(at1, at2)
            dist_mat.append([at1, at2, dist2])
    return dist_mat
#===============================================================================
def _calc_v_angle(vec1, vec2, deg=True):
    angle = arccos(clip(dot(vec1, vec2)/norm(vec1)/norm(vec2), -1., 1.))
    if deg:
        angle *= 180./pi
    return angle
