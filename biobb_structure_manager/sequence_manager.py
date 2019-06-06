
""" Module to manage sequence information for structures """


class SequenceData():
    def __init__(self):
        self.seqs = {}
        self.canonical = False
        self.fasta = []
        
        
    def add_empty_chain(self, ch_id):
        self.seqs[ch_id] = {'can':None , 'chains': [], 'pdb':{}}
        

    def load_sequence_from_fasta(self, fasta_sequence_path):
        """ Loads canonical sequence from external FASTA file"""
        self.fasta = []
        if fasta_sequence_path:
            try:
                for record in SeqIO.parse(fasta_sequence_path, 'fasta'):
                    self.fasta.append(record)
            except IOError:
                sys.exit("Error loading FASTA")

    def read_sequences(self, st, clean=True):
        """ Extracts sequences"""
        if clean:
            self.seqs = {}
            self.canonical = False
        
        if not self.canonical:
            self.read_canonical_seqs()
            
        self.read_structure_seqs()
    
    def read_canonical_seqs(self, st, ):
        """ Prepare canonical sequences """
        
        if not self.chain_ids:
            self.set_chain_ids()
            
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
                if ch_id not in self.chain_ids:
                    continue
                if ch_id not in self.sequences:
                    self.sequences[ch_id] = {'can':None , 'chains': [], 'pdb':{}}
                self.sequences[ch_id]['can'] = SeqRecord(
                    Seq(seqs[i].replace('\n', ''), IUPAC.protein),
                    'csq_' + ch_id,
                    'csq_' + ch_id,
                    'canonical sequence chain ' + ch_id
                )
                self.sequences[ch_id]['can'].features.append(
                    SeqFeature(FeatureLocation(1, len(seqs[i])))
                )

                for chn in chids[i].split(','):
                    if chn in self.chain_ids:
                        self.sequences[ch_id]['chains'].append(chn)

        self.canonical_sequence = True
        return 0
    
    def read_structure_seqs(self):
        """ Extracts sequences from structure"""
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
                if ch_id not in self.sequences:
                    self.sequences[ch_id] = {'can':None, 'chains':[], 'pdb':{}}
                self.sequences[ch_id]['pdb'][mod.id] = seqs

    def update_internals(self):
        """ Update internal data when structure is modified """

        
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
        # get canonical and structure sequences
        self._get_sequences(clean=True)

