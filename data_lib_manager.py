"""
  module to load data library json file
"""

import json

class DataLib():

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
            
    def get_atom_list(self,aa):
        return {'backbone': self.residue_data['*']['bck_atoms'], 'side': self.residue_data[aa]['side_atoms']}

    def get_atom_feature_list(self,feature):
        f_list = {}
        for rcode in self.residue_data:
            if feature in self.residue_data[rcode]:
                f_list[rcode]=self.residue_data[rcode][feature]
        return f_list
    
    def get_mutation_rules(self,aa_in,aa_out,rule_group):
        return self.residue_data[aa_in]['mutation_rules'][aa_out][rule_group]
    
    def get_distance(self, distid):
        return self.distances[distid]
    


