"""
  module to load data library json file
"""

import json
import sys

class DataLibManager():

    def __init__(self, file_path):
        try:
            fh = open (file_path, "r")
            json_map = json.load(fh)
            self.residue_codes = json_map['data_library']['residue_codes']
            self.atom_data = json_map['data_library']['atom_data']
            self.residue_data = json_map['data_library']['residue_data']
            self.distances = json_map['data_library']['distances']

        except IOError:
            print ("ERROR: unable to open data library "+ file_path, file=sys.stderr)
            sys.exit(2)

    def get_valid_codes (self, mol_type):
        if mol_type == 'na':
            codes = self.residue_codes['dna']+ self.residue_codes['rna']
        else:
            codes = self.residue_codes[mol_type]
        return codes

    def get_all_atom_lists(self):
        ats_lists = {}
        for aa in self.residue_codes['protein']:
            ats_lists[aa] = {'backbone': self.residue_data['*']['bck_atoms'], 'side': self.residue_data[aa]['side_atoms']}
        return ats_lists

    def get_metal_atoms(self):
        return self.atom_data['metal_atoms']

    def get_atom_feature_list(self,feature):
        f_list = {}
        for rcode in self.residue_data:
            if feature in self.residue_data[rcode]:
                f_list[rcode]=self.residue_data[rcode][feature]
        if not '*' in f_list:
            f_list['*'] = []
        return f_list

    def get_atom_lists(self, contact_types):
        atom_lists={}
        for cls in contact_types:
            if cls != 'severe':
                atom_lists[cls]= self.get_atom_feature_list(cls + '_atoms')
        return atom_lists

    def get_mutation_rules(self,aa_in,aa_out,rule_group):
        return self.residue_data[aa_in]['mutation_rules'][aa_out][rule_group]

    def get_distances(self, distid):
        return self.distances[distid]

    def get_amide_data(self):
        alist=[]
        rlist={}
        for rcode in self.residue_data:
            if 'amide_atoms' in self.residue_data[rcode]:
                alist += self.residue_data[rcode]['amide_atoms']
                rlist[rcode] = self.residue_data[rcode]['amide_atoms']
        return [rlist,alist]

    def get_chiral_data(self):
        rlist={}
        for rcode in self.residue_data:
            if 'chiral_atoms' in self.residue_data[rcode]:
                rlist[rcode]=self.residue_data[rcode]['chiral_atoms']
        return rlist

    def get_mutation_map(self):
        mut_rules={}
        for rcode in self.residue_data:
            if 'mutation_rules' in self.residue_data[rcode]:
                mut_rules[rcode]=self.residue_data[rcode]['mutation_rules']
        return(mut_rules)
