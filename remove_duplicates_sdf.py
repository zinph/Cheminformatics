#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:     Remove duplicate compounds from sdf file based on Inchi keys.
#
# Author:      zinph
#
# Created:     14/11/2019
# Copyright:   (c) zinph 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd
from tkinter.filedialog import askopenfilename
import os

def main():
    # Choose the sdf file when prompted

    filename = askopenfilename()
    df = PandasTools.LoadSDF(filename, molColName='structure')
    num_original = df.shape[0]
    suppl = Chem.SDMolSupplier(filename)
    mols = [x for x in suppl]
    inchis = [Chem.inchi.MolToInchi(x) for x in mols]
    df['inchis'] = inchis

    df = df.drop_duplicates(subset=['inchis'])
    num_current = df.shape[0]
    num_removed = num_original - num_current
    print(str(num_removed) + ' duplicate compounds were removed.')

    df = df.loc[:, df.columns != 'inchis']
    PandasTools.WriteSDF(df,filename[:-4]+'_DuplicatesRemoved.sdf',molColName='structure',properties=list(df.columns))

if __name__ == '__main__':
    main()
