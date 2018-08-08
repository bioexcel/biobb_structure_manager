"""
    StructureManager: module to handle structure data.
"""
import re
import sys
import warnings

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.parse_pdb_header import parse_pdb_header
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from Bio import BiopythonWarning

import structure_manager.util as util

MODELS_MAXRMS = 5.0    # Threshold value to detect NMR models (angs)

class StructureManager():

    def __init__(self):
        self.model_type = ''
        self.num_ats = 0

    def loadStructure(self, input_pdb_path, use_models, remove_H, debug=False):
# TODO support biounits
        if "pdb:"in input_pdb_path:
            pdbl = PDBList(pdb='tmpPDB')

            try:
                input_pdb_path = input_pdb_path[4:].upper()
                real_pdb_path = pdbl.retrieve_pdb_file(input_pdb_path,file_format='mmCif')
                parser = MMCIFParser()
                input_format = 'cif'

            except IOError:
                sys.stderr.write('ERROR: fetching structure at {}'.format(input_pdb_path))
#                print ('ERROR: fetching structure at {}'.format(input_pdb_path), file=sys.stderr)
                sys.exit(2)
        else:
            real_pdb_path = input_pdb_path
            if '.pdb' in real_pdb_path:
                parser = PDBParser(PERMISSIVE=1)
                input_format = 'pdb'
            elif '.cif' in real_pdb_path:
                parser = MMCIFParser()
                input_format = 'cif'
            else:
#                print ('ERROR: unknown filetype', file=sys.stderr)
                sys.stderr.write('ERROR: unknown filetype')
                sys.exit(2)
        try:
            warnings.simplefilter('ignore', BiopythonWarning)
            self.st = parser.get_structure('st', real_pdb_path)
            if input_format =='pdb':
                self.headers = parse_pdb_header(real_pdb_path)
            else:
                self.headers = MMCIF2Dict(real_pdb_path)

        except OSError:
            #print ("#ERROR: parsing PDB", file=sys.stderr)
            sys.stderr.write ("#ERROR: parsing PDB")
            sys.exit(2)

        #====== Internal residue renumbering =========================================
        i = 1
        for r in self.st.get_residues():
            r.index = i
            i += 1

        #Atom renumbering for mmCIF,
        if input_format == 'cif':
            i = 1
            for at in self.st.get_atoms(): # Check numbering in models
                at.serial_number = i
                if hasattr(at, 'selected_child'):
                    at.selected_child.serial_number = i
                i += 1

        for at in self.st[0].get_atoms():
            self.num_ats += 1
        #checking models type
        if len(self.st) > 1:

            if use_models == 'no':
                print ("WARNING: Input Structure contains models, but using only first one due to use_models settings")
                remove_models = True

            elif use_models == 'auto':
                if debug:
                    print ("DEBUG: Found " + str(len(self.st)) + " models")
                    print ("DEBUG: RMS " + str(util.calcRMSdAll(self.st[0], self.st[1])))

                if util.calcRMSdAll(self.st[0], self.st[1]) < MODELS_MAXRMS:
                    if debug:
                        print ("DEBUG: Models look like alternative conformations, will consider only one")
                    self.model_type = 'alt'
                    print ("WARNING: Input Structure contains models, but they look as NMR models, using the first one (override with force)")
                    remove_models = True

                else:
                    self.model_type = 'traj'
                    remove_models = False

            elif use_models == 'force':
                if util.calcRMSdAll(self.st[0], self.st[1]) < MODELS_MAXRMS:
                    print ('#WARNING: Models found look like NMR models, but using all due to use_models = force')
                remove_models = False

            else:
#                print ("#ERROR: Unknown use_models option", file=sys.stderr)
                sys.stderr.write ("#ERROR: Unknown use_models option")
                sys.exit(1)

            if remove_models:
                print ("Removing models")
                ids = []

                for md in self.st.get_models():
                    ids.append(md.id)

                for i in range(1, len(ids)):
                    self.st.detach_child(ids[i])
        # Hydrogens

        if remove_H == 'all':
            print ("Removing H atoms")

            for r in self.st.get_residues():
                util.removeHFromRes(r)

    def get_structure(self):
        return self.st

    def saveStructure(self, output_pdb_path):

        pdbio = PDBIO()
        pdbio.set_structure(self.st)

        try:
            pdbio.save(output_pdb_path)

        except OSError:
#            print ("#ERROR: unable to save PDB data on " + output_path, file=sys.stderr)
            sys.stderr.write ("#ERROR: unable to save PDB data on " + output_path)


    def get_all_at2at_distances(self, at_ids = ['all'], d_cutoff=0.):
        dist_mat = []
        at_list = []
        for at in self.st.get_atoms():
            if at.id in at_ids or at_ids == 'all':
                at_list.append(at)
        for i in range(0, len(at_list)-1):
            for j in range(i + 1, len(at_list)):
                d = at_list[i] - at_list[j]
                if d_cutoff > 0. and d < d_cutoff:
                    dist_mat.append ([at_list[i], at_list[j], d])
        return dist_mat
    
    def get_all_r2r_distances(self, r_ids = ['all'], d_cutoff=0.):
        # Uses distances between the first atom of each residue as r-r distance
        dist_mat = []
        r_list = []
        for r in self.st.get_residues():
            if r.resname in r_ids or r_ids == 'all':
                r_list.append(r)
        for i in range(0, len(r_list)-1):
            ati = r_list[i].child_list[0]
            for j in range(i + 1, len(r_list)):
                atj = r_list[j].child_list[0]
                d = ati - atj
                if d_cutoff > 0. and d < d_cutoff:
                    dist_mat.append ([r_list[i], r_list[j], d])
        return dist_mat
        
    def get_nmodels(self):
        return len(self.st)

    def select_model(self, nm):
        self.st = self.st[nm-1]

    def get_chain_ids(self,models=False):
        chain_ids = []
        for ch in self.st.get_chains():
            id = ch.id
            if models:
                id += "/{}".format(ch.get_parent().id + 1)
            chain_ids.append(id)
        return chain_ids

    def select_chains(self, select_chains):
        chain_ids = self.get_chain_ids()
        ch_ok = select_chains.split(',')
        for ch in ch_ok:
            if not ch in chain_ids:
#                print ("Error: request chain not present", ch, file=sys.stderr)
                sys.stderr.write ('Error: requested chain {} not present'.format(ch))
                select_chains = ''
        for ch in chain_ids:
            if ch not in ch_ok:
                self.st[0].detach_child(ch)

    def get_altloc_residues(self):
        res_list = {}
        for at in self.st.get_atoms():
            r = at.get_parent()
            rid = util.residueid(r)
            if at.get_altloc() != ' ':
                if rid not in res_list:
                    res_list[rid] = []
                res_list[rid].append(at)
        return res_list

    def select_altloc_residues(self, r, select_altloc):
        alt_loc_res = self.get_altloc_residues()
        for at in alt_loc_res[r]:
            res = at.get_parent()
            if select_altloc == 'occupancy':
                newat = at.selected_child
            else:
                newat = at.child_dict[select_altloc]
            newat.disordered_flag = 0
            newat.altloc = ' '
            res.detach_child(at.id)
            res.add (newat)
        res.disordered = 0

    def get_metals(self, metal_ats):
        met_list = []
        for at in self.st.get_atoms():
            if not re.match('H_', at.get_parent().id[0]):
                continue
            if at.id in metal_ats:
                met_list.append(at)
        return met_list

    def remove_residue(self,r):
        r.get_parent().detach_child(r.id)
        
    def get_ligands(self, incl_water=False):
        lig_list = []
        for r in self.st.get_residues():
            if re.match('H_', r.id[0]) or (incl_water and re.match('W', r.id[0])):
               lig_list.append(r)
        return lig_list
        
    def get_residues_with_H(self):
        resh_list = []
        for r in self.st.get_residues():
            has_h=0
            for a in r.get_atoms():
                if a.element == 'H':
                    has_h +=1
            if has_h:
                resh_list.append({'r':r, 'n_h':has_h})
        return resh_list
        
    def invert_amide_residue(self, r):
        res_type=r.get_resname()
        if not res_type in ['ASN','GLN']:
            sys.stderr.write('Error: {} is not an amide residue',format(res_type))
        if res_type == 'ASN':
            util.swap_atom_names(r['OD1'],r['ND2'])
        else: 
            util.swap_atom_names(r['OE1'],r['NE2'])
            
