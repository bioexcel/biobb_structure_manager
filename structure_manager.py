"""
    StructureManager: module to handle structure data.
"""
from Bio import BiopythonWarning
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.parse_pdb_header import parse_pdb_header
import re
import structure_manager.model_utils as mu
import sys
import warnings


class StructureManager():

    def __init__(self, input_pdb_path, debug=False):
        self.input_format = ''
        self.model_type = ''
        self.num_ats = 0
        self.nmodels = 0
        self.chain_ids = []
        self.modified = False

# TODO support biounits
        if "pdb:"in input_pdb_path:
            pdbl = PDBList(pdb='tmpPDB')

            try:
                input_pdb_path = input_pdb_path[4:].upper()
                real_pdb_path = pdbl.retrieve_pdb_file(input_pdb_path, file_format='mmCif')
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
            if self.input_format == 'pdb':
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
        if self.nmodels > 1:
            self.models_type = mu.guess_models_type(self.st)
        else:
            self.models_type = 0

        self.set_chain_ids()
        self.calc_stats()

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

    def calc_stats(self):
        nr = 0
        na = 0
        hi = 0
        wi = 0
        insi = 0
        for r in self.st.get_residues():
            nr += 1
            if mu.is_wat(r):
                wi += 1
            if mu.is_hetatm(r):
                hi += 1
            if mu.has_ins_code(r):
                insi += 1
            for at in r.get_atoms(): # Check numbering in models
                na += 1
        self.num_res = nr
        self.num_ats = na
        self.res_insc = insi
        self.res_hetats = hi
        self.res_ligands = hi-wi
        self.num_wat = wi

    def get_stats(self):
        return {
            'nmodels': self.nmodels,
            'nchains': len(self.chain_ids),
            'chain_ids': self.chain_ids,
            'num_res': self.num_res,
            'num_ats': self.num_ats,
            'res_insc': self.res_insc,
            'res_hetats': self.res_hetats,
            'res_ligands': self.res_ligands,
            'num_wat': self.num_wat
        }
    def print_stats(self, prefix=''):
        stats = self.get_stats()
        print ('{} Num. models: {}'.format(prefix, stats['nmodels']))
        print ('{} Num. chains: {} ({})'.format(prefix, stats['nchains'], ','.join(sorted(stats['chain_ids']))))
        print ('{} Num. residues:  {}'.format(prefix, stats['num_res']))
        print ('{} Num. residues with ins. codes:  {}'.format(prefix, stats['res_insc']))
        print ('{} Num. HETATM residues:  {}'.format(prefix, stats['res_hetats']))
        print ('{} Num. ligand or modified residues:  {}'.format(prefix, stats['res_ligands']))
        print ('{} Num. water mol.:  {}'.format(prefix, stats['num_wat']))
        print ('{} Num. atoms:  {}'.format(prefix, stats['num_ats']))
    

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
        ids = []
        for md in self.st.get_models():
            ids.append(md.id)
        for i in range(0, len(ids)):
            if i != nm-1:
                self.st.detach_child(ids[i])
        if nm != 1:
            self.st[nm-1].id = 0
        self.nmodels = 1
        self.models_type = 0
        self.set_chain_ids()
        self.modified=True

    def set_chain_ids(self):
        self.chain_ids = {}
        for ch in self.st[0].get_chains():
            self.chain_ids[ch.id] = mu.guess_chain_type(ch)
        self.chains_ids = sorted(self.chain_ids)

    def select_chains(self, select_chains):
        if not self.chain_ids:
            self.set_chain_ids()
        ch_ok = select_chains.split(',')
        for ch in ch_ok:
            if not ch in self.chain_ids:
                sys.stderr.write ('Error: requested chain {} not present'.format(ch))
                select_chains = ''
        for md in self.st:
            for ch in self.chain_ids:
                if ch not in ch_ok:
                    self.st[md.id].detach_child(ch)
        self.set_chain_ids()
        self.modified=True

    def select_altloc_residues(self, r, select_altloc):
        alt_loc_res = mu.get_altloc_residues(self.st)
        for at in alt_loc_res[r]:
            res = at.get_parent()
            if select_altloc.lower() == 'occupancy':
                newat = at.selected_child
            else:
                newat = at.child_dict[select_altloc]
            newat.disordered_flag = 0
            newat.altloc = ' '
            res.detach_child(at.id)
            res.add (newat)
        res.disordered = 0
        self.modified=True

    def remove_residue(self, r):
        mu.remove_residue(r)
        self.modified=True
        
