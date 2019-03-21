"""
  Module to manage mutations
"""

import re
import sys

import structure_manager.model_utils as mu

ADD = 'Add'
DEL = 'Del'
MOV = 'Mov'

class MutationManager():

    def __init__(self, id_list, debug=False):
        self.mutation_list = []

        if 'file:' in id_list:
            #Load from file
            id_list = id_list.replace ('file:', '')
            print ('Reading mutation list from file {}'.format(id_list))
            self.id_list=[]
            for line in open(id_list, 'r'):
                self.id_list.append(line.replace('\n', '').replace('\r', ''))

        else:
            self.id_list = id_list.replace(' ', '').split(',')

        #convert to list ids to Mutation instances and make the list
        self.mutation_list = list(map (Mutation, self.id_list))
    # Check mutation_list and unroll chains/models
    def prepare_mutations (self, st, debug=False, stop_on_error=True):
        for mut in self.mutation_list:
            mut.prepare(st, debug, stop_on_error)

    # Perform the modifications according to rules
    def apply_mutations (self, st, mutation_map, residue_lib, remove_H='mut', debug=False):
        mutated_res=[]
        for mut in self.mutation_list:
            mutated_res += mut.apply(st, mutation_map, residue_lib, remove_H, debug)
        return mutated_res

    def __str__(self):
        return ','.join(self.id_list)

#===============================================================================

class Mutation():

    def __init__(self, mut_id):
        if ':' not in mut_id:
            mut_id = '*:' + mut_id

        self.id = mut_id.upper()

        [self.chain, mut] = mut_id.split(':')

        mut_comps = re.match('([A-z]*)([0-9]*)([A-z]*)', mut)

        self.old_id = mu.protein_residue_check(mut_comps.group(1))
        self.new_id = mu.protein_residue_check(mut_comps.group(3))
        self.res_num = mut_comps.group(2)

        self.id = self.chain + ":" + self.old_id + self.res_num + self.new_id

    def prepare (self, st, debug=False, stop_on_error=True): # Check which mutations are possible
        if debug:
            print ('#DEBUG: Checking ' + self.id)

        self.mutations = []
        ok = 0
        for model in st.get_models():
            if self.chain == '*':
                chains = []
                for ch in model.get_list():
                    chains.append(ch.get_id())
            else:
                chains = [self.chain]

            for ch in chains:
                if int(self.res_num) in model[ch]:
                    r = model[ch][int(self.res_num)]
                    if r.get_resname() == self.old_id:
                        self.mutations.append(
                            {
                                'model':model.get_id(),
                                'chain':ch,
                                'residue':r.get_id(),
                                'new_id':self.new_id
                            })
                        ok += 1
                    else:
                        print ('#WARNING: Unknown residue {}:{}{}'.format(ch, self.old_id,self.res_num))
                else:
                    print ('#WARNING: Unknown residue {}:{}{}',format(ch,self.old_id,self.res_num))

        if ok == 0 and stop_on_error:
            print ('#ERROR: no mutations available for {}'.format(self.id), file=sys.stderr)
            sys.exit(1)

    def apply(self, st, mut_map, res_lib, remove_H, debug=False):
        if debug:
            print (self.mutations)
            print ("Mutation Rules")
            print (mut_map[self.old_id][self.new_id])
        mutated_res=[]
        for m in self.mutations:
            r = st[m['model']][m['chain']][m['residue']]
            rname = r.get_resname().replace(' ','')
# Deleting H
            if remove_H == 'mut':
                mu.remove_H_from_r(r, verbose= True)
            # checking side chain
            side_atoms=[]
            bck_atoms=[]
            for at in r.get_atoms():
                if at.id in mut_map[rname]['side_atoms']:
                    side_atoms.append(at.id)
                else:
                    bck_atoms.append(at.id)
            for at_id in ['N','CA','C']:
                if at_id not in bck_atoms:
                    print ('#ERROR: Backbone atoms missing for {}, aborting'.format(mu.residue_id(r)), file=sys.stderr)
                    sys.exit(1)
            missing_ats = []
            for at_id in mut_map[rname]['side_atoms']:
                if at_id not in side_atoms:
                    missing_ats.append(at_id)
            print ("Replacing " + mu.residue_id(r) + " into " + self.new_id)
            in_rules = []
            extra_adds=[]
# Renaming ats
            for rule in mut_map[rname][self.new_id][MOV]:
                [old_at, new_at] = rule.split("-")
                print ('  Renaming {} to {}'.format(old_at,new_at))
                if old_at in side_atoms:
                    mu.rename_atom(r, old_at, new_at)
                else:
                    print ('#WARNING: atom {} missing in {}'.format(old_at, mu.residue_id(r)))
                    extra_adds.append(new_at)
                in_rules.append(old_at)

# Deleting atoms
            for at_id in mut_map[rname][self.new_id][DEL]:
                print ('  Deleting {}'.format(at_id))
                if at_id in side_atoms:
                    mu.delete_atom(r, at_id)
                else:
                    print ('#WARNING: atom {} already missing in {}'.format(at_id,mu.residue_id(r)))
                in_rules.append(at_id)


# Adding atoms (new_id required as r.resname is still the original)
# Adding missing atoms that keep name
            for at_id in mut_map[rname]['side_atoms']:
                if at_id not in in_rules and at_id in missing_ats:
                    print ('  Adding missing atom {}'.format(at_id))
                    mu.build_atom(r, at_id, res_lib, self.new_id)
            for at_id in extra_adds:
                print ('  Adding new atom {}'.format(at_id))
                mu.build_atom(r, at_id, res_lib, self.new_id)
            for at_id in mut_map[rname][self.new_id][ADD]:
                print ('  Adding new atom {}'.format(at_id))
                mu.build_atom(r, at_id, res_lib, self.new_id)
#Renaming residue
            r.resname = self.new_id
            mutated_res.append(r)
        print("")
        return mutated_res



    def __str__(self):
      return self.id

