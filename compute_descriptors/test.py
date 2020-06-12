# -*- coding: utf-8 -*-
"""
Created on Wed May 20 22:03:52 2020

@author: phyophyokyawzin
"""


import pandas as pd
from molvs import standardize_smiles
from ECFP6 import *


def main():

    filename = 'data/macrolides_smiles.csv'
    df = pd.read_csv(filename)
    smiles = [standardize_smiles(i) for i in df['smiles'].values]

    output_filename = 'desc' + filename[4:]

    ## Compute ECFP6 Fingerprints and export file.
    maccs_descriptor = ECFP6(smiles)
    maccs_descriptor.compute_ECFP6(output_filename)


if __name__ == '__main__':
    main()