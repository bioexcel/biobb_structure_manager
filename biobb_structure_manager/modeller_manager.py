""" Module to handle an interface to modeller """

import sys
import os
import uuid
from Bio import SeqIO
from Bio.Seq import Seq, IUPAC
from Bio.SeqRecord import SeqRecord

try:
    from modeller import *
    from modeller.automodel import *
except:
    sys.exit("Error importing modeller")

tmp_base_dir = '/tmp'
database_file = '/home/gelpi/DEVEL/BioExcel/biobb/biobb_structure_checking/test/test_modeller/pdb_95.pir'
class ModellerManager():
    def __init__(self, ch_id, seqs):
#        self.tmpdir = tmp_base_dir + "/mod" + str(uuid.uuid4())
        self.tmpdir = tmp_base_dir + "/modtest"
        self.ch_id = ch_id
        self.seqs = seqs
#        try:
#            os.mkdir(self.tmpdir)
#        except IOError as e:
#            sys.exit(e)
        self.env = environ()
        self.env.io.atom_files_directory = [self.tmpdir]        
        
    def run(self):
        tgt_seq = self.seqs['can'].seq
        # Check N-term
        pdb_seq = self.seqs['pdb'][0].seq
        nt_pos = tgt_seq.find(pdb_seq)
        tgt_seq = tgt_seq[nt_pos:]
        
        for i in range(1, len(self.seqs['pdb'])):
            gap_len = self.seqs['pdb'][i].features[0].location.start\
                - self.seqs['pdb'][i-1].features[0].location.end - 1 
            pdb_seq += '-'*gap_len
            pdb_seq += self.seqs['pdb'][i].seq
        #pdb_seq += '-'*(len(tgt_seq) - len(pdb_seq))
        tgt_seq = tgt_seq[0:len(pdb_seq)]
        alin = [
            SeqRecord(
                tgt_seq, 
                'target', 
                '', 
                'sequence:target:::::::0.00: 0.00'
            ),
            SeqRecord(
                pdb_seq, 
                'templ', 
                '', 
                'structureX:templ.pdb:{}:{}:{}:{}:::-1.00: -1.00'.format(
                    self.seqs['pdb'][0].features[0].location.start, 
                    self.ch_id,
                    self.seqs['pdb'][-1].features[0].location.end, 
                    self.ch_id
                )
            )
        ]
        
        SeqIO.write(alin, self.tmpdir + "/alin.pir", 'pir')

        a = automodel(
            self.env, 
            alnfile=self.tmpdir + "/alin.pir",
            knowns='templ',
            sequence='target',
            assess_methods=(assess.DOPE,assess.GA341)
        )
        a.starting_model = 1
        a.ending_model = 1
        orig_dir = os.getcwd()
        os.chdir(self.tmpdir)
        a.make()
        os.chdir(orig_dir)
        return a.outputs[0]

