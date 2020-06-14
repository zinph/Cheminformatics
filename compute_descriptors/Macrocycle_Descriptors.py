# -*- coding: utf-8 -*-
"""
Created on Wed May 20 22:04:24 2020

@author: phyophyokyawzin
"""

import itertools
import pandas as pd
from rdkit import Chem
from mordred import Calculator, descriptors
from mordred.RingCount import RingCount

class Macrocycle_Descriptors:

    def __init__(self, smiles):
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles
        self.mordred = None


    def compute_ringsize(self, mol):
        '''
        check for macrolides of RS 3 to 99, return a  list of ring counts.
        [RS3,RS4,.....,RS99]
        [0,0,0,...,1,...,0]
        '''
        RS_3_99 = [i+3 for i in range(97)]
        RS_count = []
        for j in RS_3_99:
            RS = RingCount(order=j)(mol)
            RS_count.append(RS)
        return RS_count

    def macrolide_ring_info(self):
        headers = ['n'+str(i+13)+'Ring' for i in range(87)]+['SmallestRS','LargestRS']
        # up to nR12 is already with mordred, start with nR13 to nR99
        ring_sizes = []
        for i in range(len(self.mols)):
            RS = self.compute_ringsize(self.mols[i])  # nR3 to nR99
            RS_12_99 = RS[9:]    # start with nR12 up to nR99
            ring_indices = [i for i,x in enumerate(RS_12_99) if x!=0]  # get index if item isn't equal to 0
            # if there is a particular ring present, the frequency won't be zero. Find those indexes. 
			if ring_indices:
                # Add 12 (starting ring count) to get up to the actual ring size
                smallest_RS = ring_indices[0]+12     # Retrieve the first index (for the smallest core RS - note the list is in ascending order)
                largest_RS = ring_indices[-1]+12	 # Retrieve the last index (for the largest core RS)
                RS_12_99.append(smallest_RS)  # Smallest RS
                RS_12_99.append(largest_RS)  # Largest RS
            else:
                RS_12_99.extend(['',''])
            ring_sizes.append(RS_12_99[1:]) # up to nR12 is already with mordred, start with nR13 to nR99
        df = pd.DataFrame(ring_sizes, columns=headers)
        return df

    def sugar_count(self):
        sugar_patterns = [
        '[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)]',
        '[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]',
        '[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C(O)@C1)]',
        '[OX2;$([r5]1@C(!@[OX2H1])@C@C@C1),$([r6]1@C(!@[OX2H1])@C@C@C@C1)]',
        '[OX2;$([r5]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]',
        '[OX2;$([r5]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]',
        ]
        sugar_mols = [Chem.MolFromSmarts(i) for i in sugar_patterns]
        sugar_counts = []
        for i in self.mols:
            matches_total = []
            for s_mol in sugar_mols:
                raw_matches = i.GetSubstructMatches(s_mol)
                matches = list(sum(raw_matches, ()))
                if matches not in matches_total and len(matches) !=0:
                    matches_total.append(matches)
            sugar_indices = set((list(itertools.chain(*matches_total))))
            count = len(sugar_indices)
            sugar_counts.append(count)
        df = pd.DataFrame(sugar_counts, columns=['nSugars'])
        return df

    def core_ester_count(self):
        '''
        Returns pandas frame containing the count of esters in core rings of >=12 membered macrocycles.
        '''
        ester_smarts = '[CX3](=[OX1])O@[r;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r11]'
        core_ester = []
        ester_mol = Chem.MolFromSmarts(ester_smarts)
        for i in self.mols:
            ester_count = len(i.GetSubstructMatches(ester_mol))
            core_ester.append(ester_count)
        df = pd.DataFrame(core_ester, columns=['core_ester'])
        return df

    def mordred_compute(self, name):
        calc = Calculator(descriptors, ignore_3D=True)
        df = calc.pandas(self.mols)
        self.mordred = df
        df.insert(loc=0, column='smiles', value=self.smiles)
        df.to_csv(name[:-4]+'_mordred.csv', index=False)

    def compute_mordred_macrocycle(self, name):
        if not isinstance(self.mordred, pd.DataFrame):
            self.mordred = self.mordred_compute(name)
        ring_df = self.macrolide_ring_info()
        sugar_df = self.sugar_count()
        ester_df = self.core_ester_count()
#        self.mrc = pd.concat([ring_df,sugar_df, ester_df], axis=1)
        mordred_mrc = pd.concat([self.mordred, ring_df,sugar_df, ester_df], axis=1)
        mordred_mrc.to_csv(name[:-4]+'_mordred_mrc.csv', index=False)

