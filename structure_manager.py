"""
    StructureManager: module to handle structure data.
"""
import sys
import warnings

from Bio import BiopythonWarning
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.parse_pdb_header import parse_pdb_header

import structure_manager.model_utils as mu


class StructureManager():

    def __init__(self, input_pdb_path, debug=False):
        self.input_format = ''
        self.model_type = ''
        self.num_ats = 0
        self.nmodels = 0
        self.chain_ids = []
        self.modified = False
        self.all_residues = []

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
            self.all_residues.append(r)
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
            na += len(r.get_list()) # Check numbering in models
        self.num_res = nr
        self.num_ats = na
        self.res_insc = insi
        self.res_hetats = hi
        self.res_ligands = hi-wi
        self.num_wat = wi

    def check_missing_atoms(self, valid_codes, residue_data):
        miss_at_list = []
        for r in self.st.get_residues():
            if r.get_resname() in valid_codes and not mu.is_hetatm(r):
                miss_at = mu.check_all_at_in_r(r, residue_data[r.get_resname().replace(' ','')])
                if len(miss_at) > 0:
                    miss_at_list.append([r,miss_at])
        return miss_at_list

    def get_missing_side_chain_atoms(self, valid_codes, residue_data):
        miss_side = []
        for res in self.check_missing_atoms(valid_codes, residue_data):
            [r,at_list]=res
            if not len(at_list['backbone']):
                miss_side.append([r,at_list['side']])
        return miss_side

    def get_missing_backbone_atoms(self, valid_codes, residue_data):
        miss_bck = []
        for res in self.check_missing_atoms(valid_codes, residue_data):
            [r,at_list]=res
            if at_list['backbone']:
                miss_bck.append([r,at_list['backbone']])
        return miss_bck

    def get_bck_breaks(self):
        if not self.residue_link_to:
            self.check_backbone_connect(backbone_atoms, self.data_library.get_distances('COVLNK'))
        self.bck_breaks_list=[]
        self.wrong_link_list=[]
        self.not_link_seq_list=[]
        self.modified_residue_list=[]
        for i in range(0,len(self.all_residues)-1):
            r1 = self.all_residues[i]
            r2 = self.all_residues[i+1]
            if mu.is_hetatm(r1) or mu.is_hetatm(r2):
                if r1 in self.residue_link_to:
                    if self.residue_link_to[r1] == r2:
                        if mu.is_hetatm(r1):
                            self.modified_residue_list.append(r1)
                else:
                    continue
            #Skip NA #TODO include NA Backbone
            if self.chain_ids[r1.get_parent().id] != mu.PROTEIN:
                continue

            if r1 not in self.residue_link_to:
                self.bck_breaks_list.append([r1,r2])
                if mu.seq_consecutive(r1,r2):
                    d=0.
                    if 'N' in r1 and 'C' in r2:
                        d = r1['N']-r2['C']
                    self.not_link_seq_list.append([r1,r2,d])


            else:
                if r2 != self.residue_link_to[r1]:
                    self.wrong_link_list.append([r1,self.residue_link_to[r1],r2])
#            print (i, mu.residue_id(self.all_residues[i]))

    def check_backbone_connect(self, backbone_atoms, COVLNK):
        backbone_links = mu.get_backbone_links(self.get_structure(), backbone_atoms, COVLNK)
        self.residue_link_to = {}
        for lnk in backbone_links:
            [at1, at2] = lnk
            self.residue_link_to[at1.get_parent()]=at2.get_parent()
            # missing residues
            # diff chain
            # no n->n+1
            # missing n->n+1








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
        for i,v in enumerate(ids):
            if i != nm-1:
                self.st.detach_child(v)
        if nm != 1:
            self.st[nm-1].id = 0
        self.nmodels = 1
        self.models_type = 0
        self.set_chain_ids()
        self.modified=True

    def has_models(self):
        return self.nmodels > 1

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

    def select_altloc_residue(self, res, to_fix):
        for at in to_fix['ats']:
            if to_fix['select'].lower() == 'occupancy':
                newat = at.selected_child
            else:
                if to_fix['select'] in at.child_dict.keys():
                    newat = at.child_dict[to_fix['select']]
                else:
                    print ('Warning: unknown alternative {} in {}'.format(to_fix['select'],mu.atom_id(at)), file=sys.stderr)
                    continue
            newat.disordered_flag = 0
            newat.altloc = ' '
            res.detach_child(at.id)
            res.add (newat)
        res.disordered = 0
        self.modified=True

    def remove_residue(self, r):
        mu.remove_residue(r)
        self.modified=True

    def fix_side_chain(self,r_at, res_library):
        for at_id in r_at[1]:
            print ("  Adding new atom " + at_id)
            if at_id == 'CB':
                coords = mu.buildCoordsCB(r_at[0])
            else:
                coords = mu.buildCoordsOther(r_at[0], res_library, r_at[0].get_resname(), at_id)

            r_at[0].add(Atom(
                at_id,
                coords,
                99.0,
                1.0,
                ' ',
                ' ' + at_id + ' ',
                0,
                at_id[0:1]
                ))


        self.modified=True


