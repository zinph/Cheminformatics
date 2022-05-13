#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      pzin
#
# Created:     13/05/2022
# Copyright:   (c) pzin 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from itertools import combinations, permutations
from molvs import standardize_smiles

class AmideCoupling:
    def __init__(self):
        pass

    def pairwise(self, given_list):
        '''
        Generate all possible pairs based on a given list.
        s -> (s0,s1), (s1,s2), (s2, s3), ...
        '''
        res = list(permutations(given_list, 2))
        return res


    def make_chunks(self, given_list, n_blocks_per_container):
        '''
        Generate all possible pairs based on a given list.
        s -> (s0,s1), (s1,s2), (s2, s3), ...
        '''
        res = list(permutations(given_list, n_blocks_per_container))
        return res


    def deprotect_fmoc(self, target_smile):
        FMOC_smirks = '[#6:2][N:1]C(=O)OCC1c2ccccc2-c2ccccc12>>[#6:2][N:1]'
        fmoc_rxn = AllChem.ReactionFromSmarts(FMOC_smirks)
        reactant = Chem.MolFromSmiles(target_smile)
        products = fmoc_rxn.RunReactants((reactant,))
        try:
            deprotected_mol = products[0][0]
            deprotected_smile = Chem.MolToSmiles(deprotected_mol)
        except:
            deprotected_smile = target_smile
        return deprotected_smile


    def amide_coupling(self, smile1, smile2):
        deprotect_first_smile  = self.deprotect_fmoc(smile1)
        deprotect_second_smile = self.deprotect_fmoc(smile2)
        mol1 = Chem.MolFromSmiles(deprotect_first_smile)
        mol2 = Chem.MolFromSmiles(deprotect_second_smile)
        smarts = "[C:1](=[O:2])O.[N:3] >> [C:1](=[O:2])[N:3]"
        rxn3 = AllChem.ReactionFromSmarts (smarts)
        products = rxn3.RunReactants ([mol1, mol2])
        resulting_smile_list = []
        try:
            for i in range(len(products)):
                resulting_smile = Chem.MolToSmiles(products[i][0])
                resulting_smile_list.append(resulting_smile)
        except:
            pass
        return resulting_smile_list


    def multiple_amide_coupling(self, list_of_smiles):
        first_smile = list_of_smiles[0]
        for i in range(len(list_of_smiles)-1):
            try:
                first_smile = self.amide_coupling(first_smile, list_of_smiles[i+1])[0]
            except:
                first_smile = ''
        return first_smile


    def make_libary(self, possible_smile_blocks):
        enumerated_library    = []
        for smiles_to_piece in possible_smile_blocks:
            resulting_cpd = self.multiple_amide_coupling(list(smiles_to_piece))
            try:
                canonicalized_resulting_smile = standardize_smiles(resulting_cpd)
                if canonicalized_resulting_smile not in enumerated_library:
                    enumerated_library.append(canonicalized_resulting_smile)
            except:
                pass
        return enumerated_library


def main():
    s1 = "O=C(O)CCC1CCN(C(=O)OCC2c3ccccc3-c3ccccc32)CC1"
    s2 = "CN(CCCCC(=O)O)C(=O)OCC1c2ccccc2-c2ccccc21"
    s3 = "O=C(NC1CCN(c2ccc(C(=O)O)cn2)CC1)OCC1c2ccccc2-c2ccccc21"
    s4 = "O=C(NC1CCC(C(=O)O)CC1)OCC1c2ccccc2-c2ccccc21"
    s5 = "O=C(O)Cn1c(-c2ccccc2)cc2ccccc21"
    s6 = "N#CCn1nccc1C(=O)O"

    list_of_smiles = [s1,s2,s3,s4,s5,s6]
    desired_blocks_per_molecule = 3
    amide_coupling_object = AmideCoupling()
    possible_smile_blocks = amide_coupling_object.make_chunks(list_of_smiles, desired_blocks_per_molecule)

    print(amide_coupling_object.make_libary(possible_smile_blocks))

if __name__ == '__main__':
    main()
