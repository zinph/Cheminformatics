# -*- coding: utf-8 -*-
"""
Created on Wed May 20 22:11:12 2020

@author: phyophyokyawzin
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys

class MACCS:
    def __init__(self, smiles):
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles
        
    def compute_MACCS(self, name):
        MACCS_list = []
        header = ['bit' + str(i) for i in range(167)]
        for i in range(len(self.mols)):
            ds = list(MACCSkeys.GenMACCSKeys(self.mols[i]).ToBitString())
            MACCS_list.append(ds)
        df = pd.DataFrame(MACCS_list,columns=header)
        df.insert(loc=0, column='smiles', value=self.smiles)
        df.to_csv(name[:-4]+'_MACCS.csv', index=False)