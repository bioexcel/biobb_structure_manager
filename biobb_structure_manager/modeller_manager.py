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
    def __init__(self, ch_id, target_seq, ):
#        self.tmpdir = tmp_base_dir + "/mod" + str(uuid.uuid4())
        self.tmpdir = tmp_base_dir + "/modtest"
        self.target = None
        self.ch_id = ch_id
        self.target = target_seq
#        try:
#            os.mkdir(self.tmpdir)
#        except IOError as e:
#            sys.exit(e)
        self.env = environ()
        self.env.io.atom_files_directory = [self.tmpdir]
        
    def run(self):
        tgt = self.target
        tgt.id = 'target'
        tgt.name = ''
        tgt.description = 'sequence:target:::::::0.00: 0.00'
        SeqIO.write(tgt, self.tmpdir + "/seqs.pir", 'pir')
        
        aln = alignment(self.env)
        mdl = model(
            self.env, 
            file='templ', 
            model_segment=('FIRST:' + self.ch_id,'LAST:' + self.ch_id)
        )
        aln.append_model(mdl, align_codes='templ', atom_files='templ.pdb')
        aln.append(file=self.tmpdir + '/seqs.pir', align_codes='target')
        aln.align2d()
        aln.write(file=self.tmpdir + '/alin.pir', alignment_format='PIR')
        
        a = automodel(
            self.env, 
            alnfile=self.tmpdir + "/alin.pir",
            knowns='templ',
            sequence='target',
            assess_methods=(assess.DOPE,assess.GA341)
        )
        a.starting_model = 1
        a.ending_model = 5
        a.make()
        sys.exit()
