""" Module to handle an interface to modeller """

import sys
import os
import uuid
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

try:
    from modeller import *
    from modeller.automodel import *
except:
    sys.exit("Error importing modeller")

tmp_base_dir = '/tmp'
class ModellerManager():
    def __init__(self):
        self.tmpdir = tmp_base_dir + "/mod" + str(uuid.uuid4())
        try:
            os.mkdir(self.tmpdir)
        except IOError as e:
            sys.exit(e)
        env = environ()
        env.io.atom_files_directory = [self.tmpdir]
        
    def prepare(self):
        print(vars(self))
        alin_file = self.tmpdir + "/alin.pir"
        target_sqr = SeqRecord(self.target_seq)
        target_sqr.id = "target"
        SeqIO.write(target_sqr, alin_file, "pir")
        SeqIO.write(SeqRecord(self.template_seq), alin_file, "pir")
    def run(self):
        pass
