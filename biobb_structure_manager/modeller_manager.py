""" Module to handle an interface to modeller """

import sys
import os
import uuid
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

try:
    from modeller import *
    from modeller.automodel import *
except:
    sys.exit("Error importing modeller")

TMP_BASE_DIR = '/tmp'

class ModellerManager():
    """ Class to handle Modeller calculations """
    def __init__(self):
        self.tmpdir = TMP_BASE_DIR + "/mod" + str(uuid.uuid4())
        self.tmpdir = "/tmp/modtest"
        self.ch_id = ''
        self.seqs = ''
        try:
            os.mkdir(self.tmpdir)
        except IOError as err:
            sys.exit(err)
        self.env = environ()
        self.env.io.atom_files_directory = [self.tmpdir]
        log.none()

    def build(self, target_model, target_chain):
        """ Prepares Modeller input and builds the model """
        tgt_seq = self.seqs[target_chain]['can'].seq
        templs = []
        knowns = []
        # Check N-term
        for ch_id in self.seqs[target_chain]['chains']:
            pdb_seq = self.seqs[ch_id]['pdb'][target_model][0].seq

            if ch_id == target_chain:
                nt_pos = tgt_seq.find(pdb_seq)
                tgt_seq = tgt_seq[nt_pos:]

            for i in range(1, len(self.seqs[ch_id]['pdb'][target_model])):
                gap_len = self.seqs[ch_id]['pdb'][target_model][i].features[0].location.start\
                    - self.seqs[ch_id]['pdb'][target_model][i-1].features[0].location.end - 1
                pdb_seq += '-'*gap_len
                pdb_seq += self.seqs[ch_id]['pdb'][target_model][i].seq

            templs.append(
                SeqRecord(
                    pdb_seq,
                    'templ' + ch_id,
                    '',
                    'structureX:templ.pdb:{}:{}:{}:{}:::-1.00: -1.00'.format(
                        self.seqs[ch_id]['pdb'][target_model][0].features[0].location.start,
                        ch_id,
                        self.seqs[ch_id]['pdb'][target_model][-1].features[0].location.end,
                        ch_id
                    )
                )
            )
            knowns.append('templ' + ch_id)
            if ch_id == target_chain:
                tgt_seq = tgt_seq[0:len(pdb_seq)]
        SeqIO.write(
            [
                SeqRecord(
                    tgt_seq,
                    'target',
                    '',
                    'sequence:target:::::::0.00: 0.00'
                )
            ] + templs, self.tmpdir + "/alin.pir", 'pir')

        amdl = automodel(
            self.env,
            alnfile=self.tmpdir + "/alin.pir",
            knowns=knowns,
            sequence='target',
            assess_methods=(assess.DOPE, assess.GA341)
        )
        amdl.starting_model = 1
        amdl.ending_model = 1

        orig_dir = os.getcwd()
        os.chdir(self.tmpdir)
        amdl.make()
        os.chdir(orig_dir)

        return amdl.outputs[0]

    def __del__(self):
        shutil.rmtree(self.tmpdir)
