"""
  Module to load data library json file
"""

import json
import sys

class DataLibManager():
    """ Manages projects' global data file
    """

    def __init__(self, file_path):
        try:
            data_file_h = open(file_path, "r")
            json_map = json.load(data_file_h)
            self.residue_codes = json_map['data_library']['residue_codes']
            self.atom_data = json_map['data_library']['atom_data']
            self.residue_data = json_map['data_library']['residue_data']
            self.distances = json_map['data_library']['distances']
            self.ion_res = json_map['data_library']['addH_check_residues']
            self.std_ion = json_map['data_library']['addH_std_ion']

        except IOError:
            print("ERROR: unable to open data library " + file_path, file=sys.stderr)
            sys.exit(2)

    def get_valid_codes(self, mol_type):
        """ Obtain valid residue codes """
        if mol_type == 'na':
            codes = self.residue_codes['dna'] + self.residue_codes['rna']
        else:
            codes = self.residue_codes[mol_type]
        return codes

    def get_all_atom_lists(self):
        """ Obtains lists of atoms per protein residue. """
        ats_lists = {}
        for code in self.residue_codes['protein']:
            ats_lists[code] = {
                'backbone': self.residue_data['*']['bck_atoms'],
                'side': self.residue_data[code]['side_atoms']
            }
        return ats_lists

    def get_metal_atoms(self):
        """ Gets metal atom names. """
        return self.atom_data['metal_atoms']

    def get_atom_feature_list(self, feature):
        """ Gets a residue list with a specific section of data . """
        f_list = {}
        for rcode in self.residue_data:
            if feature in self.residue_data[rcode]:
                f_list[rcode] = self.residue_data[rcode][feature]
        if '*' not in f_list:
            f_list['*'] = []
        return f_list

    def get_atom_lists(self, contact_types):
        """ Gets a list of atoms organized per contact types. """
        atom_lists = {}
        for cls in contact_types:
            if cls != 'severe':
                atom_lists[cls] = self.get_atom_feature_list(cls + '_atoms')
        return atom_lists

    def get_mutation_rules(self, aa_in, aa_out, rule_group):
        """ Gets mutation rules for a specific mutation, and rule type """
        return self.residue_data[aa_in]['mutation_rules'][aa_out][rule_group]

    def get_distances(self, distid):
        """ Get distance constant. """
        return self.distances[distid]

    def get_amide_data(self):
        """ Gets data related to amide residues """
        alist = []
        rlist = {}
        for rcode in self.residue_data:
            if 'amide_atoms' in self.residue_data[rcode]:
                alist += self.residue_data[rcode]['amide_atoms']
                rlist[rcode] = self.residue_data[rcode]['amide_atoms']
        return [rlist, alist]

    def get_chiral_data(self):
        """ Gets data related to chiral residues. """
        rlist = {}
        for rcode in self.residue_data:
            if 'chiral_atoms' in self.residue_data[rcode]:
                rlist[rcode] = self.residue_data[rcode]['chiral_atoms']
        return rlist

    def get_mutation_map(self):
        """ Gets the complete map of mutation rules per residue."""
        mut_rules = {}
        for rcode in self.residue_data:
            if 'mutation_rules' in self.residue_data[rcode]:
                mut_rules[rcode] = self.residue_data[rcode]['mutation_rules']
                mut_rules[rcode]['side_atoms'] = self.residue_data[rcode]['side_atoms']
        return mut_rules

    def get_hydrogen_atoms(self):
        """ Gets list of hydrogen atoms per residue. """
        h_atom_list = {}
        for rcode in self.residue_data:
            if 'hydrogen_atoms' in self.residue_data[rcode]:
                h_atom_list[rcode] = self.residue_data[rcode]['hydrogen_atoms']
        return h_atom_list

    def get_ion_data(self):
        """ Gets residues data related to ionization states. """
        return self.std_ion

    def get_add_h_rules(self):
        """ Gets rules for adding Hydrogen atoms to residues. """
        rules_list = {}
        for rcode in self.residue_data:
            if 'addH_rules' in self.residue_data[rcode]:
                rules_list[rcode] = self.residue_data[rcode]['addH_rules']
        return rules_list
