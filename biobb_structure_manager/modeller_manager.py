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
    def __init__(self):
#        self.tmpdir = tmp_base_dir + "/mod" + str(uuid.uuid4())
        self.tmpdir = tmp_base_dir + "/modtest"
        self.ch_id = ''
        self.seqs = ''
#        try:
#            os.mkdir(self.tmpdir)
#        except IOError as e:
#            sys.exit(e)
        self.env = environ()
        self.env.io.atom_files_directory = [self.tmpdir]        
        
    def build(self, target_chain):
        tgt_seq = self.seqs[target_chain]['can'].seq
        templs = []
        knowns = []
        # Check N-term
        for ch_id in self.seqs[target_chain]['chains']:
            pdb_seq = self.seqs[ch_id]['pdb'][0].seq
            
            if ch_id == target_chain:
                nt_pos = tgt_seq.find(pdb_seq)
                tgt_seq = tgt_seq[nt_pos:]
        
            for i in range(1, len(self.seqs[ch_id]['pdb'])):
                gap_len = self.seqs[ch_id]['pdb'][i].features[0].location.start\
                    - self.seqs[ch_id]['pdb'][i-1].features[0].location.end - 1 
                pdb_seq += '-'*gap_len
                pdb_seq += self.seqs[ch_id]['pdb'][i].seq
        
            templs.append(
                SeqRecord(
                    pdb_seq, 
                    'templ' + ch_id, 
                    '', 
                    'structureX:templ.pdb:{}:{}:{}:{}:::-1.00: -1.00'.format(
                        self.seqs[ch_id]['pdb'][0].features[0].location.start, 
                        ch_id,
                        self.seqs[ch_id]['pdb'][-1].features[0].location.end, 
                        ch_id
                    )
                )
            )
            knowns.append('templ' + ch_id)
            if ch_id == target_chain:
                tgt_seq = tgt_seq[0:len(pdb_seq)]
        SeqIO.write (
            [
                SeqRecord(
                    tgt_seq, 
                    'target', 
                    '', 
                    'sequence:target:::::::0.00: 0.00'
                )
            ] + templs, self.tmpdir + "/alin.pir", 'pir')

        a = automodel(
            self.env, 
            alnfile=self.tmpdir + "/alin.pir",
            knowns=knowns,
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

