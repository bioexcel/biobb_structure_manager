"""Module to manage structure, based on BioPython Bio.PDB
"""
import sys
import warnings
import os
from Bio import BiopythonWarning
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.parse_pdb_header import parse_pdb_header
from mmb_server import MMBPDBList

import biobb_structure_manager.model_utils as mu

class StructureManager():
    """Main Class wrapping Bio.PDB structure object
    """

    def __init__(self, input_pdb_path, pdb_server='ftp://ftp.wwpdb.org', cache_dir='tmpPDB', file_format='mmCif'):
        """Class constructor. Sets an empty object and loads a structure
        according to parameters

        Args:
            **input_pdb_path** (str): path to input structure either in pdb or
            mmCIF format. Format is taken from file extension.
            Alternatively **pdb:pdbId** fetches the mmCIF file from RCSB

            **pdb_server** (str) : **default** for Bio.PDB defaults (RCSB), **mmb** for MMB PDB API

        Object structure:
            {
                "input_format" (str): Format of input file pdb|cif
                "model_type" (int): Guessed model type when num. models > 1 NMR|Bunit
                "num_ats" (int): Total Number of atoms
                "nmodels" (int): Number of Models
                "chain_ids" (List): Chain composition as
                    [chain_id:{}]
                "modified" (Boolean): Flag to indicated that structure has been modified
                "all_residues" (list): List of pointer to Bio.PDB.Residue objects, ordered
                    acording to input file
                "num_res" (int): Number of residues
                "res_insc" (int): Number of residues with insertion codes
                "res_hetats" (int): Number of residues flagged as HETATM
                "res_ligands" (int): Number of non water residues flagged as HETATM
                "num_wat" (int): Number of water residues
                "ca_only" (boolead): Flag to indicate a possible CA-only structure
                "bck_breaks_list" []: List of backbone breaks as [r1,r2] tuples
                "wrong_link_list" []: List of wrongly connected residues as [r1, r2, r3] tuples
                    where r1-r2 is the actual link whereas r1-r3 was expected
                "not_link_seq_list" []: List of sequence consecutive residues that are not linked as
                    [r1,r2] tuples. Most probable due to bond distance above the Covalent distance
                "modified_residue_list" []: List of residues being connected HETAM as PDB.Bio.Residue
                "backbone_links" []: List of found canonical backbone links as [at1, at2] tuples, according to
                    a distance criterium
                "residue_link_to" {}: Dict r1:r2 for pairs of linked residues in the N-term to C-term direction
                "cis_backbone_list" []: List of cis backbone links as [r1,r2,dihedral] tuples
                "lowtrans_backbone_list" []: List of trans backbone links with reduced dihedral as [r1,r2,dihedral] tuples

        """

        self.input_format = ''
        self.model_type = ''
        self.num_ats = 0
        self.nmodels = 0
        self.chain_ids = []
        self.modified = False
        self.all_residues = []
        self.biounit=False
        self.file_format=file_format
        if "pdb:"in input_pdb_path:
            pdbl = MMBPDBList(pdb=cache_dir, server=pdb_server) # MMBPDBList child defaults to Bio.PDB.PDBList if MMB server is not selected
            try:
                if '.' in input_pdb_path:
                    [pdbid,biounit]=input_pdb_path.split('.')
                    input_pdb_path = pdbid[4:].upper()
                    if pdb_server != 'mmb':
                        print ("Error: Biounits supported only on mmb server", file=sys.stderr)
                        sys.exit(1)
                    real_pdb_path = pdbl.retrieve_pdb_file(input_pdb_path, file_format='pdb', biounit=biounit)
                    self.biounit=biounit
                else:
                    input_pdb_path = input_pdb_path[4:].upper()
                    real_pdb_path = pdbl.retrieve_pdb_file(input_pdb_path, file_format=self.file_format)
                    if file_format=='pdb': # change file name to id.pdb
                        os.rename(real_pdb_path,input_pdb_path+".pdb")
                        real_pdb_path = input_pdb_path + ".pdb"
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
            print ('ERROR: unknown filetype', file=sys.stderr)
            sys.exit(2)
        try:
            warnings.simplefilter('ignore', BiopythonWarning)
            self.st = parser.get_structure('st', real_pdb_path)
            if self.input_format == 'pdb':
                self.headers = parse_pdb_header(real_pdb_path)
            else:
                self.headers = MMCIF2Dict(real_pdb_path)

        except OSError:
            print ("#ERROR: parsing PDB", file=sys.stderr)
            sys.exit(2)
        # Add .index field for correlative, unique numbering of residues
        self.residue_renumbering()
    #Atom renumbering for mmCIF, PDB uses atom number in file
        if self.input_format == 'cif' or self.biounit:
            self.atom_renumbering()

    #checking models type according to RMS among models
        self.nmodels = len(self.st)
        if self.nmodels > 1:
            self.models_type = mu.guess_models_type(self.st)
        else:
            self.models_type = 0
    #Identify chains and guess type (protein, dna, ...)
        self.set_chain_ids()
    #Calc general stats (num residues, atoms, etc)
        self.calc_stats()
    # Guess HETATM types
        self.guess_hetatm()

    def residue_renumbering(self):
        """Sets the Bio.PDB.Residue.index attribute to residues for a unique,
        consecutive, residue number, and generates the corresponding list
        in **all_residues**
        """
        i = 1
        self.all_residues=[]
        for r in self.st.get_residues():
            r.index = i
            if type(r).__name__ == 'DisorderedResidue':
                for ch_r in r.child_dict:
                    r.child_dict[ch_r].index = i
            self.all_residues.append(r)
            i += 1

    def atom_renumbering(self):
        """Sets  Bio.PDB.Atom.serial_number for all atoms in the structure, overriding original if any
        """
        i = 1
        for at in self.st.get_atoms():
            at.serial_number = i
            if hasattr(at, 'selected_child'):
                at.selected_child.serial_number = i
            i += 1

    def guess_hetatm(self):
        """ Guesses HETATM type as modified res, metal, wat, organic
        """
        self.hetatm={}
        for ty  in [mu.UNKNOWN, mu.MODRES, mu.METAL, mu.ORGANIC, mu.COVORGANIC, mu.WAT]:
            self.hetatm[ty]=[]
        for r in self.get_structure().get_residues():
            if not mu.is_hetatm(r):
                continue
            if mu.is_wat(r):
                self.hetatm[mu.WAT].append(r)
            elif len(r) == 1:
                self.hetatm[mu.METAL].append(r)
            elif 'N' in r or 'C' in r: # modified aminoacid candidate, TODO check connectivity with n-1 or n+1
                self.hetatm[mu.MODRES].append(r)
            #TODO check modified nucleotides
            else:
                self.hetatm[mu.ORGANIC].append(r)

    def calc_stats(self):
        """Calculates general statistics about the structure, and guesses whether
        it is a CA-only structure
        """
        nr = 0
        na = 0
        hi =     0
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
            na += len(r.get_list())
        self.num_res = nr
        self.num_ats = na
        self.res_insc = insi
        self.res_hetats = hi
        self.res_ligands = hi-wi
        self.num_wat = wi
        # Detecting whether it is a CA-only structure
        # num_ats should be much larger than num_res
        # waters removed
        # Taking polyGly as a lower limit
        self.ca_only = na - wi < (nr-wi)*4

    def get_ins_codes(self):
        """Makes a list with residues having insertion codes"""

        self.ins_codes_list=[]
        for r in self.st.get_residues():
            if mu.has_ins_code(r):
                self.ins_codes_list.append(r)

    def check_missing_atoms(self, valid_codes, residue_data):
        """ Makes a **list of missing atoms** in the structure

            Args:
                **valid_codes**: Valid residue 3-letter codes
                **residue_data**: Atom composition per residue taken
                    from residue library
            Returns:
                List of residues with missing atoms, as a tuples
                ["r",{"backbone":[atoms_list],"side":[atoms_list]}]
        """

        miss_at_list = []
        for r in self.st.get_residues():
            if r.get_resname() in valid_codes and not mu.is_hetatm(r):
                miss_at = mu.check_all_at_in_r(r, residue_data[r.get_resname().replace(' ','')])
                if self.is_C_term(r) and 'OXT' not in r:
                    if not 'backbone' in miss_at:
                        miss_at['backbone']=[]
                    miss_at['backbone'].append('OXT')
                if len(miss_at) > 0:
                    miss_at_list.append([r,miss_at])
        return miss_at_list

    def get_missing_side_chain_atoms(self, valid_codes, residue_data):
        """
            returns list of residues with complete backbone but with missing
            side chain atoms

            Args:
                **valid_codes**: list of valid 3-letter aminoacid codes
                **residue_data**: atom composition per residue taken
                    from residue library
            Returns:
                List of residues with missing atoms, as a tuples like
                ["r",[atoms_list]]
        """
        miss_side = []
        for res in self.check_missing_atoms(valid_codes, residue_data):
            [r,at_list]=res
            if 'side' in at_list and 'N' in r and 'CA' in r and 'C' in r:
                miss_side.append([r,at_list['side']])
        return miss_side

    def get_ion_res_list(self, ion_res, hydrogen_lists):
        """
            returns list of residues with potencial selection on adding H

            Args:
                **ion_res**: list of residue codes with potential selection
                **hydrogen_lists**: hydrogen atom names from residue library
            Returns:
                List of residues that require selection on adding H
                ["r",[atom_list]]
        """
        ion_res_list = []
        for res in self.all_residues:
            rcode=res.get_resname()
            if rcode in ion_res:
                ion_res_list.append([res,hydrogen_lists[rcode]])
        return ion_res_list

    def get_missing_backbone_atoms(self, valid_codes, residue_data):
        """
            returns list of residues with missing backbone atoms

            Args:
                **valid_codes**: list of valid 3-letter aminoacid codes
                **residue_data**: atom composition per residue taken
                    from residue library
            Returns:
                List of residues with missing atoms, as tuples
                ["r",[atoms_list]]
        """
        miss_bck = []
        for res in self.check_missing_atoms(valid_codes, residue_data):
            [r,at_list]=res
            if at_list['backbone']:
                miss_bck.append([r,at_list['backbone']])
        return miss_bck

    def get_backbone_breaks(self):
        """
            Determines several backbone anomalies
            1. Modified residues
            2. Backbone breaks
            3. Backbone links not corresponding to sequence
            4. Too long residue links
        """
        self.bck_breaks_list=[]
        self.wrong_link_list=[]
        self.not_link_seq_list=[]
        self.modified_residue_list=[]
        for i in range(0,len(self.all_residues)-1):
            r1 = self.all_residues[i]
            r2 = self.all_residues[i+1]
            if not mu.same_chain(r1,r2):
                continue
            if mu.is_hetatm(r1) or mu.is_hetatm(r2):
                if r1 in self.next_residue:
                    if self.next_residue[r1] == r2:
                        if mu.is_hetatm(r1):
                            self.modified_residue_list.append(r1)
                else:
                    continue
            #Skip NA #TODO include NA Backbone
            if self.chain_ids[r1.get_parent().id] != mu.PROTEIN:
                continue

            if r1 not in self.next_residue:
                self.bck_breaks_list.append([r1,r2])
                if mu.seq_consecutive(r1,r2):
                    d=0.
                    if 'N' in r1 and 'C' in r2:
                        d = r1['N']-r2['C']
                    self.not_link_seq_list.append([r1,r2,d])


            else:
                if r2 != self.next_residue[r1]:
                    self.wrong_link_list.append([r1,self.next_residue[r1],r2])

    def check_backbone_connect(self, backbone_atoms, COVLNK):
        """
        Determines backbone links usign a distance criterium and produces a dict with
        link relationships in the N-term to C-term direction

        Args:
            backbone_atoms: atoms to be considered as backbone
            COVLNK: Threshold distance for a covalent bond
        """
        if not hasattr(self,'backbone_links'):
            self.backbone_links = mu.get_backbone_links(self.get_structure(), backbone_atoms, COVLNK)
        self.next_residue = {}
        self.prev_residue = {}
        for lnk in self.backbone_links:
            [at1, at2] = lnk
            r1 = at1.get_parent()
            r2 = at2.get_parent()
            self.prev_residue[r2]=r1
            self.next_residue[r1]=r2

    def check_cis_backbone(self, COVLNK, check_models=True):
        """
        Determines omega dihedrals for two bound residues and classifies them
        as normal trans, low angle trans, and cis

        Args:
            backbone_atoms: atoms to be considered as backbone
            COVLNK: Distance threshold for a covalent bond
            check_models: Consider models as independent molecules

        Internal parameters:
            CISTHRES (20): Max dihedral for cis bonds
            TRANSTHRES (160): Min dihedral for trans bonds
        """
        CISTHRES = 20  #TODO check vaules withpdb checking
        TRANSTHRES = 160
        backbone_atoms = ['N','C']
        if not hasattr(self,'backbone_links'):
            self.backbone_links = mu.get_backbone_links(self.get_structure(), backbone_atoms, COVLNK)
        self.cis_backbone_list = []
        self.lowtrans_backbone_list = []
        for lnk in self.backbone_links:
            [at1, at2] = lnk
            r1 = at1.get_parent()
            r2 = at2.get_parent()
            if 'CA' in r1 and 'C' in r1 and 'CA' in r2 and 'N' in r2:
                dih =  mu.calc_bond_dihedral(r1['CA'],r1['C'],r2['N'],r2['CA'])
                if abs(dih) < CISTHRES:
                    self.cis_backbone_list.append([r1,r2,dih])
                elif abs(dih) < TRANSTHRES:
                    self.lowtrans_backbone_list.append([r1,r2,dih])

    def get_stats(self):
        """
         Returns a dict with calculates statistics

         Returns:
            Dict as {}
        """
        return {
            'nmodels': self.nmodels,
            'models_type': self.models_type,
            'nchains': len(self.chain_ids),
            'chain_ids': self.chain_ids,
            'num_res': self.num_res,
            'num_ats': self.num_ats,
            'res_insc': self.res_insc,
            'res_hetats': self.res_hetats,
            'res_ligands': self.res_ligands,
            'num_wat': self.num_wat,
            'ca_only': self.ca_only,
            'biounit': self.biounit
        }
    def print_headers(self):
        """
        Prints selected components from structure headers
        """
        self.get_headers()
        if 'entry_id'in self.meta:
            print (' PDB id: {}'.format(self.meta['entry_id']))
        print (' Title: {}'.format(self.meta['title']))
        print (' Experimental method: {}'.format( self.meta['method']))
        if 'keywords' in self.meta:
            print (' Keywords: {}'.format(self.meta['keywords']))
        if 'resolution' in self.meta:
            print (' Resolution: {} A'.format(self.meta['resolution']))
        if self.biounit:
            print (' Biounit no. {}'. format(self.meta['biounit']))

    def get_headers(self):
        """
        Extract selected components from structure headers
        """
        self.meta={}
        if self.input_format== 'cif':
            self.meta['entry_id'] = self.headers['_entry.id']
            self.meta['title']= self.headers['_struct.title']
            self.meta['method']=self.headers['_exptl.method']
            self.meta['keywords'] = self.headers['_struct_keywords.pdbx_keywords']
            if '_refine_hist.d_res_high' in self.headers:
                self.meta['resolution']=self.headers['_refine_hist.d_res_high']
        else:
            self.meta['title']= self.headers['name']
            self.meta['method']= self.headers['structure_method']
            if 'keywords' in self.headers:
                self.meta['keywords']=self.headers['keywords']
            if 'resolution' in self.headers:
                self.meta['resolution'] = self.headers['resolution']
        if self.biounit:
            self.meta['biounit'] =self.biounit


    def print_stats(self, prefix=''):
        """
        Prints statistics to stdout

        Args:
            prefix: Text prefix to prepend to printed data
        """
        stats = self.get_stats()
        if stats['nmodels'] > 1:
            print ('{} Num. models: {} (type: {}, {:8.3f} A)'.format(
                prefix,
                stats['nmodels'],
                mu.model_type_labels[stats['models_type']['type']],
                stats['models_type']['rmsd'])
            )
        else:
            print ('{} Num. models: {}'.format(prefix, stats['nmodels']))
        chids=[]
        for ch_id in sorted(stats['chain_ids']):
            if isinstance(stats['chain_ids'][ch_id],list):
                chids.append('{}: Unknown '.format(ch_id))
            else:
                chids.append('{}: {}'.format(ch_id, mu.chain_type_labels[stats['chain_ids'][ch_id]]))
        print ('{} Num. chains: {} ({})'.format(prefix, stats['nchains'], ', '.join(chids)))
        print ('{} Num. residues:  {}'.format(prefix, stats['num_res']))
        print ('{} Num. residues with ins. codes:  {}'.format(prefix, stats['res_insc']))
        print ('{} Num. HETATM residues:  {}'.format(prefix, stats['res_hetats']))
        print ('{} Num. ligands or modified residues:  {}'.format(prefix, stats['res_ligands']))
        print ('{} Num. water mol.:  {}'.format(prefix, stats['num_wat']))
        print ('{} Num. atoms:  {}'.format(prefix, stats['num_ats']))
        if stats['ca_only']:
            print ('Possible CA-Only structure')
        if len(self.hetatm[mu.MODRES]):
            print ('Modified residues found')
            for r in self.hetatm[mu.MODRES]:
                print (mu.residue_id(r))
        if len(self.hetatm[mu.METAL]):
            print ('Metal/Ion residues found')
            for r in self.hetatm[mu.METAL]:
                print (mu.residue_id(r))
        if len(self.hetatm[mu.ORGANIC]):
            print ('Small mol ligands found')
            for r in self.hetatm[mu.ORGANIC]:
                print (mu.residue_id(r))


    def get_structure(self):
        """
        Method pointing to the enclosed Bio.PDB.Structure object

        Returns:
            a Bio.PDB.Structure object
        """
        return self.st

    def save_structure(self, output_pdb_path):
        """
        Saves structure on disk in PDB format

        Args:
            output_pdb_path: OS path to the output file

        Errors:
            OSError: Error saving the file
        """
        if output_pdb_path:
            pdbio = PDBIO()
            pdbio.set_structure(self.st)
            try:
                pdbio.save(output_pdb_path)

            except OSError:
                sys.stderr.write ("#ERROR: unable to save PDB data on " + output_pdb_path)
        else:
            sys.stderr.write ("Error: output_pdb_path not provided \n")
            sys.exit(1)

# Methods to modify structure
    def select_model(self, nm):
        """
        Selects one model and delete the others from the structure. Model is
        renumbered to model 1

        Args:
            nm: Model number to keep
        """
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
        self.residue_renumbering()
        self.atom_renumbering()
        self.set_chain_ids()
        self.guess_hetatm()
        self.modified=True

    def has_models(self):
        """
        Shotcut method to check whether the structure has more than one model

        Returns: Boolean
        """
        return self.nmodels>1

    def has_superimp_models(self):
        """
        Shotcut method to check whether the structure has superimposed models (i.e. NMR or ensemble)

        Returns: Boolean
        """
        return self.models_type and self.models_type['type']==mu.ENSM

    def is_biounit(self):
        """Shortcut to check whether the structure has been loaded from a biounit
        """

    def set_chain_ids(self):
        """
        Identifies and sets the chain ids, guessing its nature (protein, dna, rna, ...)
        """
        self.chain_ids = {}
        for ch in self.st.get_chains():
            if not self.biounit and ch.get_parent().id > 0:
                continue
            self.chain_ids[ch.id] = mu.guess_chain_type(ch)
        self.chains_ids = sorted(self.chain_ids)

    def select_chains(self, select_chains):
        """
        Select one or more chains and remove the remaining.
        Args:
            select_chains: Comma separated chain ids
        """
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
        self.atom_renumbering()
        self.residue_renumbering()
        self.guess_hetatm()
        self.modified=True

    def select_altloc_residue(self, res, to_fix):
        """
        Selects one alternative conformation when altloc exists. Selection is
        done on the occupancy basis or from the conformation id.
        All relevant atoms in the residue **res** at modified in the same way.
        Triggers **modified** flag.

        Args:
            res: residue (as Bio.PDB.Residue)
            to_fix: atoms in residue to fix as
                {'ats':[atom_list],'select':'occupancy|conf_id'
        """
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
        """
        Removes residue **r** from the structure. Triggers **modified** flag
        and atom and residue renumbering
        """
        htm = mu.is_hetatm(r)
        mu.remove_residue(r)
        self.atom_renumbering()
        self.residue_renumbering()
        if htm:
            self.guess_hetatm()
        self.modified=True

    def fix_side_chain(self,r_at, res_library):
        """
        Fix missing side chain atoms in given residue. Triggers **modified** flag

        Args:
            **r_at**: tuple as [Bio.PDB.Residue, [list of atom ids]]
            **res_library**: Residue Library as structure_manager.residue_lib_manager.ResidueLib
        """
        for at_id in r_at[1]:
            print ("  Adding new atom " + at_id)
            if at_id == 'CB':
                coords = mu.build_coords_CB(r_at[0])
            else:
                coords = mu.build_coords_from_lib(r_at[0], res_library, r_at[0].get_resname(), at_id)
            mu.add_new_atom_to_residue(r_at[0], at_id, coords)
        self.atom_renumbering()
        self.modified=True

    def fix_backbone_atoms(self,r_at):
        """Adding missing backbone atoms not affecting main-chain like O and OXT
                Args:
            **r_at**: tuple as [Bio.PDB.Residue, [list of atom ids]]
            **res_library**: Residue Library as structure_manager.residue_lib_manager.ResidueLib
        """
        [r,at_list] = r_at
        print (mu.residue_id(r))
        if not 'C' in r:
            print ("Warning: not enough backbone to reconstruct missing atoms at "+ mu.residue_id(r), file=sys.stderr)
            print ("  Warning: not enough backbone to reconstruct missing atoms")
            return False
        if len(at_list) == 2 or at_list==['O']:
            if not 'CA' in r or not 'N' in r or not 'C' in r:
                print ("Warning: not enough backbone atoms to build O", file=sys.stderr)
                return 1
            print ("  Adding new atom O")
            mu.add_new_atom_to_residue(r,'O',mu.build_coords_O(r))
        if 'OXT' in at_list:
            if not 'CA' in r or not 'C' in r or not 'O' in r:
                print ("Warning: not enough backbone atoms to build OXT", file=sys.stderr)
                return 1
            print ("  Adding new atom OXT")
            mu.add_new_atom_to_residue(r,'OXT',mu.build_coords_SP2(1.229, r['C'], r['CA'], r['O']))

        self.atom_renumbering()
        self.modified=True
        return True

    def add_hydrogens(self,r_at_list, hydrogens_list, res_library, backbone_atoms, COVLNK, addH_rules, remove_H=True):
        """
        Add hydrogens according to selection in r_at_list

        Args:
           **r_at_list**: tuple as [Bio.PDB.Residue, Tauromeric Option]
           **hydrogens_list**: Hydrogen atom names per residue type
        """
        done_side = set()

        if remove_H:
            for r in self.all_residues:
                mu.remove_H_from_r (r, verbose=False)
        # residues with alternative forms
        for r_at in r_at_list:
            [r,opt] = r_at
        # Skip residues without addH rules
            if r.get_resname() not in addH_rules:
                print ("Warning: addH rules not available for residue type ", r.get_resname(),file=sys.stderr)
                continue
            if mu.is_hetatm(r):
                continue
            rcode=r.get_resname()
            self.add_hydrogens_side(r, res_library, opt, addH_rules[rcode][opt])
            r.resname=opt
            done_side.add(r)

        for r in self.all_residues:
            if mu.is_hetatm(r):
                continue
            rcode=r.get_resname()
            if not r in self.prev_residue:
                prev_residue=None
            else:
                prev_residue=self.prev_residue[r]
            self.add_hydrogens_backbone(r, prev_residue, res_library)

            if r not in done_side and rcode != 'GLY':
                if rcode not in addH_rules:
                    print ("Warning: addH rules not available for residue type ", rcode, file=sys.stderr)
                    continue
                self.add_hydrogens_side(r, res_library, rcode, addH_rules[rcode])

        self.residue_renumbering()
        self.atom_renumbering()
        self.modified=True

    def add_hydrogens_backbone(self, r, r_1, res_library):
        if 'N' not in r or 'CA'not in r or 'C' not in r:
            print ("Warning: Incomplete backbone in "+ mu.residue_id(r), file=sys.stderr)
            return 1

        if r_1 == None:
            # Nterminal TODO  Neutral NTerm
            if r.get_resname() == 'PRO':
                mu.add_new_atom_to_residue(r,'H',mu.build_coords_SP2(1.08, r['N'],r['CA'],r['CD']))
            else:
                crs = mu.build_coords_3xSP3(1.010, r['N'], r['CA'], r['C'])
                mu.add_new_atom_to_residue(r,'H1',crs[0])
                mu.add_new_atom_to_residue(r,'H2',crs[1])
                mu.add_new_atom_to_residue(r,'H3',crs[2])
        elif r.get_resname() != 'PRO':
            mu.add_new_atom_to_residue(r,'H',mu.build_coords_SP2(1.08, r['N'],r['CA'],r_1['C']))

        if r.get_resname() == 'GLY':
            crs = mu.build_coords_2xSP3(1.010, r['CA'],r['N'],r['C'])
            mu.add_new_atom_to_residue(r, 'HA2',crs[0])
            mu.add_new_atom_to_residue(r, 'HA3',crs[1])
        else:
            if 'CB' not in r:
                print ("Warning: Incomplete residue (CB atom) in "+ mu.residue_id(r),file=sys.stderr)
                return 1
            mu.add_new_atom_to_residue(r, 'HA',mu.build_coords_1xSP3(1.08,r['CA'],r['N'],r['C'],r['CB']))

    def add_hydrogens_side(self, r, res_library, opt, rules):
        if 'N' not in r or 'CA' not in r or 'C' not in r:
            print ("Warning: Incomplete backbone in "+ mu.residue_id(r),file=sys.stderr)
            return 1
        rcode=r.get_resname()
        for kr in rules.keys():
            rule=rules[kr]
            if rule['mode'] == 'B2':
                crs= mu.build_coords_2xSP3(
                    rule['dist'],
                    r[kr],
                    r[rule['ref_ats'][0]],
                    r[rule['ref_ats'][1]]
                )
                mu.add_new_atom_to_residue(r, rule['ats'][0],crs[0])
                mu.add_new_atom_to_residue(r, rule['ats'][1],crs[1])
            elif rule['mode'] == "B1":
                crs = mu.build_coords_1xSP3(
                    rule['dist'],
                    r[kr],
                    r[rule['ref_ats'][0]],
                    r[rule['ref_ats'][1]],
                    r[rule['ref_ats'][2]]
                )
                mu.add_new_atom_to_residue(r, rule['ats'][0],crs)
            elif rule['mode'] == 'S2':
                crs = mu.build_coords_SP2(
                    rule['dist'],
                    r[kr],
                    r[rule['ref_ats'][0]],
                    r[rule['ref_ats'][1]],
                )
                mu.add_new_atom_to_residue(r, rule['ats'][0],crs)
            elif rule['mode'] == 'L':
                for at_id in rule['ats']:
                    crs = mu.build_coords_from_lib(r, res_library, opt, at_id)
                    mu.add_new_atom_to_residue(r, at_id,crs)

    def is_N_term(self,r):
        return r not in self.prev_residue

    def is_C_term(self,r):
        return r not in self.next_residue
