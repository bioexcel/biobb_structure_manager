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

import structure_manager.model_utils as mu

MODELS_MAXRMS = 5.0    # Threshold value to detect NMR models (angs)

class StructureManager():

    def __init__(self, input_pdb_path, debug=False):
        self.input_format=''
        self.model_type = ''
        self.num_ats = 0
        self.nmodels = 0    
        self.chain_ids=[]
        
# TODO support biounits
        if "pdb:"in input_pdb_path:
            pdbl = PDBList(pdb='tmpPDB')

            try:
                input_pdb_path = input_pdb_path[4:].upper()
                real_pdb_path = pdbl.retrieve_pdb_file(input_pdb_path,file_format='mmCif')
                parser = MMCIFParser()
                self.input_format = 'cif'

            except IOError:
                sys.stderr.write("ERROR: fetching structure at {}\n".format(input_pdb_path))
#                print ('ERROR: fetching structure at {}'.format(input_pdb_path), file=sys.stderr)
                sys.exit(2)
        else:
            real_pdb_path = input_pdb_path
            if '.pdb' in real_pdb_path:
                parser = PDBParser(PERMISSIVE=1)
                self.input_format = 'pdb'
            elif '.cif' in real_pdb_path:
                parser = MMCIFParser()
                self.input_format = 'cif'
            else:
#                print ('ERROR: unknown filetype', file=sys.stderr)
                sys.stderr.write("ERROR: unknown filetype\n")
                sys.exit(2)
        try:
            warnings.simplefilter('ignore', BiopythonWarning)
            self.st = parser.get_structure('st', real_pdb_path)
            if self.input_format =='pdb':
                self.headers = parse_pdb_header(real_pdb_path)
            else:
                self.headers = MMCIF2Dict(real_pdb_path)

        except OSError:
            #print ("#ERROR: parsing PDB", file=sys.stderr)
            sys.stderr.write ("#ERROR: parsing PDB\n")
            sys.exit(2)
        
        self.residue_renumbering()
        
        #Atom renumbering for mmCIF,
        if self.input_format == 'cif':
            self.atom_renumbering()

        #checking models type
        self.nmodels = len(self.st)
        self.models_type=mu.guess_models_type(self.st, MODELS_MAXRMS)
        
        self.set_chain_ids()
        self.set_num_ats()
        
    def residue_renumbering(self):
        i = 1
        for r in self.st.get_residues():
            r.index = i
            i += 1
    
    def atom_renumbering(self):
        i = 1
        for at in self.st.get_atoms(): # Check numbering in models
            at.serial_number = i
            if hasattr(at, 'selected_child'):
                at.selected_child.serial_number = i
            i += 1
        
    def set_num_ats(self):
        i = 0
        hi = 0
        wi = 0
        for at in self.st.get_atoms(): # Check numbering in models
            i += 1
            if mu.is_hetatm(at.get_parent()):
                hi += 1
            if mu.is_wat(at.get_parent()):
                wi += 1
        self.num_ats = i
        self.num_hetats = hi
        self.num_wats = wi
    
    def get_structure(self):
        return self.st

    def save_structure(self, output_pdb_path):
        if output_pdb_path:
            pdbio = PDBIO()
            pdbio.set_structure(self.st)
            try:
                pdbio.save(output_pdb_path)

            except OSError:
                sys.stderr.write ("#ERROR: unable to save PDB data on " + output_path)
        else:
            sys.stderr.write ("Error: output_pdb_path not provided \n")
            sys.exit(1)

    def select_model(self, nm):
        ids=[]
        for md in self.st.get_models():
            ids.append(md.id)
        for i in range(1, len(ids)):
            if i != nm-1:
                self.st.detach_child(ids[i])
        self.nmodels=1
        self.residue_renumbering()
        self.atom_renumbering()
        self.set_num_ats()
        self.set_chain_ids()

    def set_chain_ids(self):
        self.chain_ids = []
        for ch in self.st.get_chains():
            id = ch.id
            if self.nmodels > 1:
                id += "/{}".format(ch.get_parent().id + 1)
            self.chain_ids.append(id)
        self.chains_ids = sorted(self.chain_ids)
        
    def select_chains(self, select_chains):
        if not self.chain_ids:
            self.set_chain_ids()
        ch_ok = select_chains.split(',')
        for ch in ch_ok:
            if not ch in self.chain_ids:
                sys.stderr.write ('Error: requested chain {} not present'.format(ch))
                select_chains = ''
        for ch in self.chain_ids:
            if ch not in ch_ok:
                self.st[0].detach_child(ch)
        self.set_chain_ids()
        self.residue_renumbering()
        self.atom_renumbering()
        self.set_num_ats()

    def select_altloc_residues(self, r, select_altloc):
        alt_loc_res = mu.get_altloc_residues(self.st)
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
        self.atom_renumbering()
        self.set_num_ats()            
        
    def remove_residue(self,r):
        mu.remove_residue(r)
        self.atom_renumbering()
        self.residue_renumbering()
        self.set_num_ats()            
    