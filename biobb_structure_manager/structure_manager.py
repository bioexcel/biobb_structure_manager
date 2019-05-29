"""Module to manage structure, based on BioPython Bio.PDB
"""
import warnings
import os
import sys
from Bio import BiopythonWarning
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.parse_pdb_header import parse_pdb_header
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Superimposer import Superimposer
from Bio.Seq import Seq, IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from biobb_structure_manager.mmb_server import MMBPDBList
from biobb_structure_manager.mutation_manager import MutationManager
from biobb_structure_manager.data_lib_manager import DataLibManager
from biobb_structure_manager.residue_lib_manager import ResidueLib

import biobb_structure_manager.model_utils as mu

CISTHRES = 20  #TODO check vaules with pdb checking
TRANSTHRES = 160

ALL_CONTACT_TYPES = [
    'severe',
    'apolar',
    'polar_acceptor',
    'polar_donor',
    'positive',
    'negative'
]
AMIDE_CONTACT_TYPES = [
    'polar_acceptor',
    'polar_donor',
]

class StructureManager():
    """Main Class wrapping Bio.PDB structure object
    """

    def __init__(
            self,
            input_pdb_path,
            data_library_path,
            res_library_path,
            pdb_server='ftp://ftp.wwpdb.org',
            cache_dir='tmpPDB',
            file_format='mmCif',
            fasta_sequence_path=''):

        """Class constructor. Sets an empty object and loads a structure
        according to parameters

        Args:
            **input_pdb_path** (str): path to input structure either in pdb or
            mmCIF format. Format is taken from file extension.
            Alternatively **pdb:pdbId** fetches the mmCIF file from RCSB
            **data_library_path**: Path to json data library
            **res_library_path**: Path to residue library
            **pdb_server** (str) :  **default** for Bio.PDB defaults (RCSB),
                                    **mmb** for MMB PDB API
            **cache_dir**: path to temporary dir to store downloaded structures
            **file_format**: structure file format to use
            **fasta_sequence_path**: path to canonical sequence file (needed for PDB input)

        Object structure:
            {
                "input_format" (str): Format of input file pdb|cif
                "models_type" (int): Guessed model type when num. models > 1 NMR|Bunit
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
                "modified_residue_list" []: List of residues being connected
                                            HETATM as PDB.Bio.Residue
                "backbone_links" []: List of found canonical backbone links
                                    as [at1, at2] tuples, according to a distance criterium

        """

        self.chain_ids = {}

        self.backbone_links = []
        self.modified_residue_list = []

        self.hetatm = {}
        self.num_res = 0
        self.num_ats = 0
        self.res_hetats = 0
        self.num_wat = 0
        self.res_insc = 0
        self.res_ligands = 0
        self.ca_only = False

        self.ss_bonds = []

        self.meta = {}

        self.all_residues = []
        self.next_residue = {}
        self.prev_residue = {}

        self.sequences = {}

        self.modified = False
        self.biounit = False
        self.fixed_side = False
        self.file_format = file_format

        self.data_library = DataLibManager(data_library_path)
        self.res_library = ResidueLib(res_library_path)

        self.input_format = self._load_structure_file(
            input_pdb_path, cache_dir, pdb_server, file_format
        )
        self.fasta = []

        if fasta_sequence_path:
            for record in SeqIO.parse(fasta_sequence_path, 'fasta'):
                self.fasta.append(record)

        #checking models type according to RMS among models
        self.nmodels = len(self.st)
        self.models_type = mu.guess_models_type(self.st) if self.nmodels > 1 else 0

        # Calc internal data
        self.update_internals()

    def _load_structure_file(self, input_pdb_path, cache_dir, pdb_server, file_format):
        if "pdb:"in input_pdb_path:
            # MMBPDBList child defaults to Bio.PDB.PDBList if MMB server is not selected
            pdbl = MMBPDBList(pdb=cache_dir, server=pdb_server)
            if '.' in input_pdb_path:
                [pdbid, biounit] = input_pdb_path.split('.')
                input_pdb_path = pdbid[4:].upper()
                if pdb_server != 'mmb':
                    raise WrongServerError
                real_pdb_path = pdbl.retrieve_pdb_file(
                    input_pdb_path, file_format='pdb', biounit=biounit
                )
                self.biounit = biounit
            else:
                input_pdb_path = input_pdb_path[4:].upper()
                real_pdb_path = pdbl.retrieve_pdb_file(
                    input_pdb_path, file_format=self.file_format
                )
                if file_format == 'pdb':
                    # change file name to id.pdb
                    os.rename(real_pdb_path, input_pdb_path + ".pdb")
                    real_pdb_path = input_pdb_path + ".pdb"
        else:
            real_pdb_path = input_pdb_path

        if '.pdb' in real_pdb_path:
            parser = PDBParser(PERMISSIVE=1)
            input_format = 'pdb'
        elif '.cif' in real_pdb_path:
            parser = MMCIFParser()
            input_format = 'cif'
        else:
            raise UnknownFileTypeError(input_pdb_path)

        warnings.simplefilter('ignore', BiopythonWarning)

        self.st = parser.get_structure('st', real_pdb_path)

        if input_format == 'pdb':
            self.headers = parse_pdb_header(real_pdb_path)
        else:
            self.headers = MMCIF2Dict(real_pdb_path)

        return input_format

    def _get_sequences(self):
        """ Extract sequences from structure, requires mmCIF or external fasta, only protein"""
        self.sequences = {}

        if self.fasta:
            chids = []
            seqs = []
            for rec in self.fasta:
                chids.append(rec.id.split('_')[1])
                seqs.append(str(rec.seq))
        else:
            if self.input_format != 'cif':
                print("Warning: sequence features only available in mmCIF" +\
                    " format or with external fasta input")
                return 1
            #TODO check for NA

            if not isinstance(self.headers['_entity_poly.pdbx_strand_id'], list):
                chids = [self.headers['_entity_poly.pdbx_strand_id']]
                seqs = [self.headers['_entity_poly.pdbx_seq_one_letter_code_can']]
            else:
                chids = self.headers['_entity_poly.pdbx_strand_id']
                seqs = self.headers['_entity_poly.pdbx_seq_one_letter_code_can']

        for i in range(0, len(chids)):
            for ch_id in chids[i].split(','):
                self.sequences[ch_id] = {
                    'can' : SeqRecord(
                        Seq(seqs[i].replace('\n', ''), IUPAC.protein),
                        'csq_' + ch_id,
                        'csq_' + ch_id,
                        'canonical sequence chain ' + ch_id
                    ),
                    'chains': chids[i].split(','),
                    'pdb': {}
                }
                self.sequences[ch_id]['can'].features.append(
                    SeqFeature(FeatureLocation(1, len(seqs[i])))
                )

        # PDB extrated sequences
        for mod in self.st:
            ppb = PPBuilder()
            for chn in mod.get_chains():
                seqs = []
                #self.sequences[ch_id]['pdb'][mod.id] = [1]
                ch_id = chn.id
                for frag in ppb.build_peptides(chn):
                    start = frag[0].get_id()[1]
                    end = frag[-1].get_id()[1]
                    frid = '{}:{}-{}'.format(ch_id, start, end)
                    sqr = SeqRecord(
                        frag.get_sequence(),
                        'pdbsq_' + frid,
                        'pdbsq_' + frid,
                        'PDB sequence chain ' + frid
                    )
                    sqr.features.append(SeqFeature(FeatureLocation(start, end)))
                    seqs.append(sqr)
                self.sequences[ch_id]['pdb'][mod.id] = seqs
        return 0

    def update_internals(self):
        """ Update internal data when structure is modified """

        # get canonical and structure sequences
        self._get_sequences()

        # Add .index field for correlative, unique numbering of residues
        self.residue_renumbering()
        #Atom renumbering for mmCIF, PDB uses atom number in file
        self.atom_renumbering()
        self.set_chain_ids()
        self.calc_stats()
        self.guess_hetatm()
        #Pre_calc rr distances with separated models
        self.rr_dist = mu.get_all_r2r_distances(
            self.st,
            'all',
            self.data_library.distances['R_R_CUTOFF'],
            join_models=False
        )
        #Precalc backbone . TODO Nucleic Acids
        self.check_backbone_connect(
            ('N', 'C'),
            self.data_library.distances['COVLNK']
        )

    def residue_renumbering(self):
        """Sets the Bio.PDB.Residue.index attribute to residues for a unique,
        consecutive, residue number, and generates the corresponding list
        in **all_residues**
        """
        i = 1
        self.all_residues = []
        for res in self.st.get_residues():
            res.index = i
            if type(res).__name__ == 'DisorderedResidue':
                for ch_r in res.child_dict:
                    res.child_dict[ch_r].index = i
            self.all_residues.append(res)
            i += 1

    def atom_renumbering(self):
        """ Sets  Bio.PDB.Atom.serial_number for all atoms in the structure,
            overriding original if any.
        """
        i = 1
        for atm in self.st.get_atoms():
            atm.serial_number = i
            if hasattr(atm, 'selected_child'):
                atm.selected_child.serial_number = i
            i += 1

    def guess_hetatm(self):
        """ Guesses HETATM type as modified res, metal, wat, organic
        """
        self.hetatm = {}
        for typ  in [mu.UNKNOWN, mu.MODRES, mu.METAL, mu.ORGANIC, mu.COVORGANIC, mu.WAT]:
            self.hetatm[typ] = []
        for res in self.st.get_residues():
            if not mu.is_hetatm(res):
                continue
            if mu.is_wat(res):
                self.hetatm[mu.WAT].append(res)
            elif len(res) == 1:
                self.hetatm[mu.METAL].append(res)
            elif 'N' in res or 'C' in res:
                # modified aminoacid candidate, TODO check connectivity with n-1 or n+1
                self.hetatm[mu.MODRES].append(res)
                #TODO check modified nucleotides
            else:
                self.hetatm[mu.ORGANIC].append(res)

    def calc_stats(self):
        """ Calculates general statistics about the structure, and guesses whether
            it is a CA-only structure
        """
        self.num_res = 0
        self.num_ats = 0
        self.res_hetats = 0
        self.num_wat = 0
        self.res_insc = 0
        for res in self.st.get_residues():
            self.num_res += 1
            if mu.is_wat(res):
                self.num_wat += 1
            if mu.is_hetatm(res):
                self.res_hetats += 1
            if mu.has_ins_code(res):
                self.res_insc += 1
            self.num_ats += len(res.get_list())
        self.res_ligands = self.res_hetats - self.num_wat
        # Detecting whether it is a CA-only structure
        # num_ats should be much larger than num_res
        # waters removed
        # Taking polyGly as a lower limit
        self.ca_only = self.num_ats - self.num_wat < (self.num_res - self.num_wat)*4

    def get_ins_codes(self):
        """Makes a list with residues having insertion codes"""
        return [
            res
            for res in self.st.get_residues()
            if mu.has_ins_code(res)
        ]

    def get_metal_atoms(self):
        """ Makes a list of possible metal atoms"""
        return mu.get_metal_atoms(self.st, self.data_library.atom_data['metal_atoms'])

    def get_SS_bonds(self):
        """ Stores and returns possible SS Bonds by distance"""
        self.ss_bonds = mu.get_all_at2at_distances(
            self.st,
            'SG',
            self.data_library.distances['SS_DIST'],
            not self.has_superimp_models()
        )
        return self.ss_bonds

    def check_chiral_sides(self):
        """ Returns a list of wrong chiral side chains"""
        chiral_res = self.data_library.get_chiral_data()
        chiral_list = [
            res
            for res in self.st.get_residues()
            if res.get_resname() in chiral_res
        ]

        if not chiral_list:
            return {}

        return {
            'list' : chiral_list,
            'res_to_fix' : [
                res
                for res in chiral_list
                if not mu.check_chiral_residue(res, chiral_res)
            ]
        }

    def get_chiral_bck_list(self):
        """ Returns a list of residues with chiral CAs"""
        prot_chains = 0
        chiral_bck_list = []
        for chn in self.st.get_chains():
            if self.chain_ids[chn.id] == mu.PROTEIN:
                prot_chains += 1
                for res in chn.get_residues():
                    if res.get_resname() != 'GLY' and not mu.is_hetatm(res):
                        chiral_bck_list.append(res)

        if not prot_chains:
            print("No protein chains detected, skipping")
            return {}

        if not chiral_bck_list:
            print("No residues with chiral CA found, skipping")
            return {'chiral_bck_list': []}

        return {
            'list': chiral_bck_list,
            'res_to_fix' : [res for res in chiral_bck_list if not mu.check_chiral_ca(res)]
        }

    def check_r_list_clashes(self, residue_list, contact_types):
        """ Checks clashes originated by a list of residues"""
        return mu.check_r_list_clashes(
            residue_list,
            self.rr_dist,
            self.data_library.distances['CLASH_DIST'],
            self.data_library.get_atom_lists(contact_types),
            not self.has_superimp_models()
        )

    def check_rr_clashes(self, res1, res2, contact_types):
        """ Checks all clashes between two residues"""
        return mu.check_rr_clashes(
            res1,
            res2,
            self.data_library.distances['CLASH_DIST'],
            self.data_library.get_atom_lists(contact_types)
        )

    def check_missing_atoms(self):
        """ Makes a **list of missing atoms** in the structure

            Returns:
                List of residues with missing atoms, as a tuples
                ["r",{"backbone":[atoms_list],"side":[atoms_list]}]
        """
        #TODO Nucleic acids
        valid_codes = self.data_library.get_valid_codes('protein')
        residue_data = self.data_library.get_all_atom_lists()
        miss_at_list = []
        for res in self.st.get_residues():
            if res.get_resname() in valid_codes and not mu.is_hetatm(res):
                miss_at = mu.check_all_at_in_r(
                    res, residue_data[res.get_resname().replace(' ', '')]
                )
                if self.is_C_term(res) and 'OXT' not in res:
                    if 'backbone' not in miss_at:
                        miss_at['backbone'] = []
                    miss_at['backbone'].append('OXT')
                if miss_at:
                    miss_at_list.append((res, miss_at))
        return miss_at_list

    def check_extra_atoms(self):
        """ Makes a **list of extra atoms** in the structure

            Returns:
                List of residues with extra atoms, as a tuples
                ["r",atoms_list]
        """
        #TODO Nucleic acids
        valid_codes = self.data_library.get_valid_codes('protein')
        residue_data = self.data_library.get_all_atom_lists()
        extra_at_list = []
        for res in self.st.get_residues():
            if res.get_resname() in valid_codes and not mu.is_hetatm(res):
                extra_ats = mu.check_unk_at_in_r(
                    res, residue_data[res.get_resname().replace(' ', '')]
                )
                if extra_ats:
                    extra_at_list.append((res, extra_ats))
        return extra_at_list

    def get_missing_atoms(self, fragment):
        """ Returns list of residues with missing atoms

            Returns:
                List of residues with missing atoms, as a tuples like
                ["r",[atoms_list]]
        """
        miss_ats = []
        for res_at in self.check_missing_atoms():
            res, at_list = res_at
            if fragment == 'side':
                if 'side' in at_list and 'N' in res and 'CA' in res and 'C' in res:
                    miss_ats.append((res, at_list['side']))
            else:
                if at_list['backbone']:
                    miss_ats.append((res, at_list['backbone']))
        return miss_ats

    def get_ion_res_list(self):
        """
            returns list of residues with potencial selection on adding H

            Returns:
                List of residues that require selection on adding H
                ["r",[atom_list]]
        """
        ion_res = self.data_library.ion_res
        hydrogen_lists = self.data_library.get_hydrogen_atoms()

        ion_res_list = []
        for res in self.all_residues:
            rcode = res.get_resname()
            if rcode in ion_res:
                ion_res_list.append((res, hydrogen_lists[rcode]))
        return ion_res_list

    def get_backbone_breaks(self):
        """
            Determines several backbone anomalies
            1. Modified residues
            2. Backbone breaks
            3. Backbone links not corresponding to sequence
            4. Too long residue links
        """
        bck_breaks_list = []
        wrong_link_list = []
        not_link_seq_list = []
        self.modified_residue_list = []
        for i in range(0, len(self.all_residues)-1):
            res1 = self.all_residues[i]
            res2 = self.all_residues[i+1]
            if not mu.same_chain(res1, res2):
                continue
            if mu.is_hetatm(res1) or mu.is_hetatm(res2):
                if res1 in self.next_residue:
                    if self.next_residue[res1] == res2:
                        if mu.is_hetatm(res1):
                            self.modified_residue_list.append(res1)
                else:
                    continue
            #Skip NA #TODO include NA Backbone
            if self.chain_ids[res1.get_parent().id] != mu.PROTEIN:
                continue

            if res1 not in self.next_residue:
                bck_breaks_list.append([res1, res2])
                if mu.seq_consecutive(res1, res2):
                    dist = 0.
                    if 'N' in res1 and 'C' in res2:
                        dist = res1['N'] - res2['C']
                    not_link_seq_list.append([res1, res2, dist])

            else:
                if res2 != self.next_residue[res1]:
                    wrong_link_list.append([res1, self.next_residue[res1], res2])
        return {
            'bck_breaks_list': bck_breaks_list,
            'wrong_link_list': wrong_link_list,
            'not_link_seq_list': not_link_seq_list
        }

    def check_backbone_connect(self, backbone_atoms, covlnk):
        """
        Determines backbone links usign a distance criterium and produces a dict with
        link relationships in the N-term to C-term direction

        Args:
            backbone_atoms: atoms to be considered as backbone
            COVLNK: Threshold distance for a covalent bond
        """
        self.backbone_links = mu.get_backbone_links(
            self.st, backbone_atoms, covlnk
        )
        self.next_residue = {}
        self.prev_residue = {}
        for lnk in self.backbone_links:
            [at1, at2] = lnk
            res1 = at1.get_parent()
            res2 = at2.get_parent()
            self.prev_residue[res2] = res1
            self.next_residue[res1] = res2

    def check_cis_backbone(self):
        """
        Determines omega dihedrals for two bound residues and classifies them
        as normal trans, low angle trans, and cis

        Args:
            backbone_atoms: atoms to be considered as backbone
            covlnk: Distance threshold for a covalent bond

        """
        cis_backbone_list = []
        lowtrans_backbone_list = []
        for lnk in self.backbone_links:
            [at1, at2] = lnk
            res1 = at1.get_parent()
            res2 = at2.get_parent()
            if 'CA' in res1 and 'C' in res1 and 'CA' in res2 and 'N' in res2:
                dih = mu.calc_bond_dihedral(res1['CA'], res1['C'], res2['N'], res2['CA'])
                if abs(dih) < CISTHRES:
                    cis_backbone_list.append((res1, res2, dih))
                elif abs(dih) < TRANSTHRES:
                    lowtrans_backbone_list.append((res1, res2, dih))
        return cis_backbone_list, lowtrans_backbone_list

    def check_amide_contacts(self):
        """ Check close contacts involving amide atoms """
        amide_res, amide_atoms = self.data_library.get_amide_data()

        amide_list = set(
            res for res in self.st.get_residues()
            if res.get_resname() in amide_res
        )

        if not amide_list:
            return {}

        c_list = self.check_r_list_clashes(
            amide_list,
            AMIDE_CONTACT_TYPES
        )
        amide_res_to_fix = []
        amide_cont_list = []
        for cls in c_list:
            for rkey in c_list[cls]:
                [at1, at2] = c_list[cls][rkey][0:2]
                res1 = at1.get_parent()
                res2 = at2.get_parent()
                add_pair = False
                if at1.id in amide_atoms and res1 in amide_list:
                    amide_res_to_fix.append(res1)
                    add_pair = True
                if at2.id in amide_atoms and res2 in amide_list:
                    amide_res_to_fix.append(res2)
                    add_pair = True
                if add_pair:
                    amide_cont_list.append(c_list[cls][rkey])
        return {
            'list':  amide_list,
            'res_to_fix': amide_res_to_fix,
            'cont_list': amide_cont_list,
        }

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
            print(' PDB id: {}'.format(self.meta['entry_id']))
        print(' Title: {}'.format(self.meta['title']))
        print(' Experimental method: {}'.format(self.meta['method']))
        if 'keywords' in self.meta:
            print(' Keywords: {}'.format(self.meta['keywords']))
        if 'resolution' in self.meta:
            print(' Resolution: {} A'.format(self.meta['resolution']))
        if self.biounit:
            print(' Biounit no. {}'. format(self.meta['biounit']))

    def get_headers(self):
        """
        Extract selected components from structure headers
        """
        self.meta = {}
        if self.input_format == 'cif':
            self.meta['entry_id'] = self.headers['_entry.id']
            self.meta['title'] = self.headers['_struct.title']
            self.meta['method'] = self.headers['_exptl.method']
            self.meta['keywords'] = self.headers['_struct_keywords.pdbx_keywords']
            if '_refine_hist.d_res_high' in self.headers:
                self.meta['resolution'] = self.headers['_refine_hist.d_res_high']
        else:
            self.meta['title'] = self.headers['name']
            self.meta['method'] = self.headers['structure_method']
            if 'keywords' in self.headers:
                self.meta['keywords'] = self.headers['keywords']
            if 'resolution' in self.headers:
                self.meta['resolution'] = self.headers['resolution']
        if self.biounit:
            self.meta['biounit'] = self.biounit

    def print_model_stats(self, prefix=''):
        """ Print stats """
        if self.nmodels > 1:
            print(
                '{} Num. models: {} (type: {}, {:8.3f} A)'.format(
                    prefix,
                    self.nmodels,
                    mu.MODEL_TYPE_LABELS[self.models_type['type']],
                    self.models_type['rmsd']
                )
            )
        else:
            print('{} Num. models: {}'.format(prefix, self.nmodels))

    def print_chain_stats(self, prefix=''):
        """ Print chains info """
        chids = []
        for ch_id in sorted(self.chain_ids):
            if isinstance(self.chain_ids[ch_id], list):
                chids.append('{}: Unknown '.format(ch_id))
            else:
                chids.append(
                    '{}: {}'.format(
                        ch_id, mu.CHAIN_TYPE_LABELS[self.chain_ids[ch_id]]
                    )
                )
        print('{} Num. chains: {} ({})'.format(prefix, len(self.chain_ids), ', '.join(chids)))


    def print_stats(self, prefix=''):
        """
        Prints statistics to stdout

        Args:
            prefix: Text prefix to prepend to printed data
        """
        stats = self.get_stats()
        self.print_model_stats(prefix)
        self.print_chain_stats(prefix)

        print('{} Num. residues:  {}'.format(prefix, stats['num_res']))
        print('{} Num. residues with ins. codes:  {}'.format(prefix, stats['res_insc']))
        print('{} Num. HETATM residues:  {}'.format(prefix, stats['res_hetats']))
        print('{} Num. ligands or modified residues:  {}'.format(prefix, stats['res_ligands']))
        print('{} Num. water mol.:  {}'.format(prefix, stats['num_wat']))
        print('{} Num. atoms:  {}'.format(prefix, stats['num_ats']))
        if stats['ca_only']:
            print('Possible CA-Only structure')
        if self.hetatm[mu.MODRES]:
            print('Modified residues found')
            for res in self.hetatm[mu.MODRES]:
                print(mu.residue_id(res))
        if self.hetatm[mu.METAL]:
            print('Metal/Ion residues found')
            for res in self.hetatm[mu.METAL]:
                print(mu.residue_id(res))
        if self.hetatm[mu.ORGANIC]:
            print('Small mol ligands found')
            for res in self.hetatm[mu.ORGANIC]:
                print(mu.residue_id(res))

    def save_structure(self, output_pdb_path, mod_id=None):
        """
        Saves structure on disk in PDB format

        Args:
            output_pdb_path: OS path to the output file
            mod_id (optional): model to write

        Errors:
            OSError: Error saving the file
        """
        if not output_pdb_path:
            raise OutputPathNotProvidedError
        pdbio = PDBIO()

        if mod_id is None:
            pdbio.set_structure(self.st)
            pdbio.save(output_pdb_path)
        else:
            pdbio.set_structure(self.st[mod_id])
            pdbio.save(output_pdb_path)


    def get_all_r2r_distances(self, res_group, join_models):
        """ Determine residue pairs within a given Cutoff distance
            calculated from the first atom available
            Args:
                res_group: list of residues to check | 'all'
                join_models: consider all models as separated molecules
            Output:
                List of tupes (r1,r2,dist)

        """
        return mu.get_all_r2r_distances(
            self.st,
            res_group,
            self.data_library.distances['R_R_CUTOFF'],
            join_models=join_models
        )

    def get_altloc_residues(self):
        return mu.get_altloc_residues(self.st)

# Methods to modify structure
    def select_model(self, keep_model):
        """ Selects one model and delete the others from the structure. Model is
            renumbered to model 1

            Args:
                keep_model: Model number to keep
        """
        ids = [nmodel.id for nmodel in self.st.get_models()]
        for ind, model_id in enumerate(ids):
            if ind != keep_model - 1:
                self.st.detach_child(model_id)
        if keep_model != 1:
            self.st[keep_model-1].id = 0
        self.nmodels = 1
        self.models_type = 0
        # Update internal data
        self.update_internals()
        self.modified = True

    def has_models(self):
        """ Shotcut method to check whether the structure has more than one model

            Returns: Boolean
        """
        return self.nmodels > 1

    def has_superimp_models(self):
        """ Shotcut method to check whether the structure has superimposed
            models (i.e. NMR or ensemble)

            Returns: Boolean
        """
        return self.models_type and self.models_type['type'] == mu.ENSM

    def set_chain_ids(self):
        """
        Identifies and sets the chain ids, guessing its nature (protein, dna, rna, ...)
        """
        self.chain_ids = {}
        for chn in self.st.get_chains():
            if not self.biounit and chn.get_parent().id > 0:
                continue
            self.chain_ids[chn.id] = mu.guess_chain_type(chn)

    def select_chains(self, select_chains):
        """
        Select one or more chains and remove the remaining.
        Args:
            select_chains: Comma separated chain ids
        """
        if not self.chain_ids:
            self.set_chain_ids()
        ch_ok = select_chains.split(',')
        for chn in ch_ok:
            if not chn in self.chain_ids:
                print('Warning: skipping unknown chain', chn)
        for mod in self.st:
            for chn in self.chain_ids:
                if chn not in ch_ok:
                    self.st[mod.id].detach_child(chn)
        # Update internal data
        self.update_internals()
        self.modified = True

    def select_altloc_residue(self, res, to_fix):
        """ Selects one alternative conformation when altloc exists. Selection is
            done on the occupancy basis or from the conformation id.
            All relevant atoms in the residue **res** at modified in the same way.
            Triggers **modified** flag.

            Args:
                res: residue (as Bio.PDB.Residue)
                to_fix: atoms in residue to fix as
                    {'ats':[atom_list],'select':'occupancy|conf_id'
        """
        for atm in to_fix['ats']:
            if to_fix['select'].lower() == 'occupancy':
                newat = atm.selected_child
            else:
                if to_fix['select'] in atm.child_dict.keys():
                    newat = atm.child_dict[to_fix['select']]
                else:
                    print(
                        'Warning: unknown alternative {} in {}'.format(
                            to_fix['select'], mu.atom_id(atm)
                        )
                    )
                    continue
            newat.disordered_flag = 0
            newat.altloc = ' '
            res.detach_child(atm.id)
            res.add(newat)
        res.disordered = 0
        self.modified = True

    def remove_residue(self, res, update_int=True):
        """
        Removes residue **r** from the structure. Triggers **modified** flag
        and atom and residue renumbering
        """
        mu.remove_residue(res)

        # Update internal data
        if update_int:
            self.update_internals()

        self.modified = True

    def fix_side_chain(self, r_at):
        """
        Fix missing side chain atoms in given residue. Triggers **modified** flag

        Args:
            **r_at**: tuple as [Bio.PDB.Residue, [list of atom ids]]
        """
        print(mu.residue_id(r_at[0]))
        for at_id in r_at[1]:
            print("  Adding new atom " + at_id)
            if at_id == 'CB':
                coords = mu.build_coords_CB(r_at[0])
            else:
                coords = mu.build_coords_from_lib(
                    r_at[0],
                    self.res_library,
                    r_at[0].get_resname(),
                    at_id
                )
            mu.add_new_atom_to_residue(r_at[0], at_id, coords)
        self.atom_renumbering()
        self.modified = True

    def fix_backbone_chain(self, brk_list, key_modeller=''):
        """ Fixes backbone breaks using Modeller """
        if key_modeller:
            os.environ['KEY_MODELLER9v21']=key_modeller
        from biobb_structure_manager.modeller_manager import ModellerManager

        ch_to_fix = set()
        for brk in brk_list:
            ch_to_fix.add(brk[0].get_parent().id)

        mod_mgr = ModellerManager()
        mod_mgr.seqs = self.sequences

        fixed_segments = []
        for mod in self.st:
            if self.has_models():
                print('Processing Model {}'.format(mod.id + 1))
                self.save_structure('{}/templ.pdb'.format(mod_mgr.tmpdir), mod.id)
            else:
                self.save_structure('{}/templ.pdb'.format(mod_mgr.tmpdir))
            for ch_id in self.chain_ids:
                if ch_id not in ch_to_fix:
                    continue
                print("Fixing backbone of chain " + ch_id)
                model_pdb = mod_mgr.build(mod.id, ch_id)
                parser = PDBParser(PERMISSIVE=1)
                model_st = parser.get_structure(
                    'model_st',
                    mod_mgr.tmpdir + "/" + model_pdb['name']
                )
                fixed_gaps = self.merge_structure(
                    model_st,
                    mod.id,
                    ch_id,
                    self.sequences[ch_id]['pdb'][mod.id][0].features[0].location.start
                )
                fixed_segments += fixed_gaps

        self.update_internals()

        return fixed_segments


    def merge_structure(self, new_st, mod_id, ch_id, offset):
        spimp = Superimposer()
        fixed_ats = [atm for atm in self.st[mod_id][ch_id].get_atoms() if atm.id == 'CA']
        moving_ats = []
        for atm in fixed_ats:
            moving_ats.append(new_st[0][' '][atm.get_parent().id[1] - offset + 1]['CA'])
        spimp.set_atoms(fixed_ats, moving_ats)
        spimp.apply(new_st.get_atoms())

        list_res = self.st[mod_id][ch_id].get_list()
        fixed_gaps = []
        for i in range(0, len(self.sequences[ch_id]['pdb'][mod_id])-1):
            gap_start = self.sequences[ch_id]['pdb'][mod_id][i].features[0].location.end
            gap_end = self.sequences[ch_id]['pdb'][mod_id][i+1].features[0].location.start
            pos = 0
            while pos < len(list_res) and self.st[mod_id][ch_id].child_list[pos].id[1] != gap_start:
                pos += 1
            self.remove_residue(self.st[mod_id][ch_id][gap_start], update_int=False)
            self.remove_residue(self.st[mod_id][ch_id][gap_end], update_int=False)
            for nres in range(gap_start, gap_end + 1):
                res = new_st[0][' '][nres - offset + 1].copy()
                res.id = (' ', nres, ' ')
                self.st[mod_id][ch_id].insert(pos, res)
                pos += 1
                print("Adding " + mu.residue_id(res))
            fixed_gaps.append(
                '{}{}-{}{}/{}'.format(ch_id, gap_start, ch_id, gap_end, mod_id + 1)
            )
            print()
        return fixed_gaps

    def fix_backbone_O_atoms(self, r_at):
        """Adding missing backbone atoms not affecting main-chain like O and OXT
                Args:
            **r_at**: tuple as [Bio.PDB.Residue, [list of atom ids]]
        """
        res, at_list = r_at
        print(mu.residue_id(res))
        if not 'C' in res:
            raise NotEnoughAtomsError
        if len(at_list) == 2 or at_list == ['O']:
            if 'CA' not in res or 'N' not in res or 'C' not in res:
                raise NotEnoughAtomsError
            print("  Adding new atom O")
            mu.add_new_atom_to_residue(res, 'O', mu.build_coords_O(res))
        if 'OXT' in at_list:
            if 'CA' not in res or 'C' not in res or 'O' not in res:
                raise NotEnoughAtomsError
            print("  Adding new atom OXT")
            mu.add_new_atom_to_residue(
                res,
                'OXT',
                mu.build_coords_SP2(mu.OINTERNALS[0], res['C'], res['CA'], res['O'])
            )

        self.atom_renumbering()
        self.modified = True
        return True

    def add_hydrogens(self, ion_res_list, remove_h=True):
        """
        Add hydrogens considering selections in ion_res_list

        Args:
           **r_at_list**: dict as Bio.PDB.Residue: Tauromeric Option
           **remove_h**: Remove Hydrogen atom before adding new ones
        """
        add_h_rules = self.data_library.get_add_h_rules()

        for res in self.all_residues:
            if mu.is_hetatm(res):
                continue

            if remove_h:
                mu.remove_H_from_r(res, verbose=False)

            if res not in self.prev_residue:
                prev_residue = None
            else:
                prev_residue = self.prev_residue[res]

            error_msg = mu.add_hydrogens_backbone(res, prev_residue)

            rcode = res.get_resname()

            if rcode == 'GLY':
                continue

            if rcode not in add_h_rules:
                print(NotAValidResidueError(rcode).message)
                continue

            if res in ion_res_list:
                if rcode != ion_res_list[res]:
                    print(
                        'Replacing {} by {}'.format(
                            mu.residue_id(res), ion_res_list[res]
                        )
                    )
                error_msg = mu.add_hydrogens_side(
                    res,
                    self.res_library,
                    ion_res_list[res],
                    add_h_rules[rcode][ion_res_list[res]]
                )
                res.resname = ion_res_list[res]
            else:
                error_msg = mu.add_hydrogens_side(res, self.res_library, rcode, add_h_rules[rcode])

            if error_msg:
                print(error_msg, mu.residue_id(res))

        self.residue_renumbering()
        self.atom_renumbering()
        self.modified = True

    def is_N_term(self, res):
        """ Detects whether it is N terminal residue."""
        return res not in self.prev_residue

    def is_C_term(self, res):
        """ Detects whether it is C terminal residue."""
        return res not in self.next_residue

    def prepare_mutations(self, mut_list):
        """ Fins residues to mutate from mut_list"""
        mutations = MutationManager(mut_list)
        mutations.prepare_mutations(self.st)
        return mutations

    def apply_mutations(self, mutations):
        """ Perform mutations """
        mutated_res = mutations.apply_mutations(
            self.data_library.get_mutation_map(),
            self.res_library
        )
        self.residue_renumbering()
        self.atom_renumbering()
        self.modified = True
        return mutated_res

    def invert_amide_atoms(self, res):
        """ Fix sidechains with incorrect amide assignments"""
        amide_res = self.data_library.get_amide_data()[0]
        res_type = res.get_resname()
        if res_type not in amide_res:
            raise NotAValidResidueError(res_type)
        mu.swap_atoms(
            res[amide_res[res_type][0]],
            res[amide_res[res_type][1]]
        )

    def fix_chiral_chains(self, res):
        """ Fix sidechains with chiral errors"""
        chiral_res = self.data_library.get_chiral_data()
        res_type = res.get_resname()
        mu.swap_atoms(
            res[chiral_res[res_type][0]],
            res[chiral_res[res_type][1]]
        )
        if res_type == 'ILE':
            mu.delete_atom(res, 'CD1')
            mu.build_atom(res, 'CD1', self.res_library, 'ILE')

#===============================================================================

class Error(Exception):
    """ Base class """
    pass

class WrongServerError(Error):
    def __init__(self):
        self.message = 'ERROR: Biounits supported only on MMB server'

class UnknownFileTypeError(Error):
    def __init__(self, type):
        self.message = 'ERROR: unknown filetype ({})'.format(type)

class OutputPathNotProvidedError(Error):
    def __init__(self):
        self.message = 'ERROR: output PDB path not provided'

class NotAValidResidueError(Error):
    def __init__(self, res):
        self.message = 'Warning: {} is not a valid residue in this context'.format(res)

class NotEnoughAtomsError(Error):
    def __init__(self):
        self.message = 'Warning: not enough backbone to build missing atoms'
