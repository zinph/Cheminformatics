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
    df = PandasTools.LoadSDF(filename, molColName='structure')  # Load sdf file to pandas
    num_original = df.shape[0]                             # get number of original compounds

    suppl = Chem.SDMolSupplier(filename)                   # load sdf file with rdkit
    mols = [x for x in suppl]                              # get mols from sdf file
    inchis = [Chem.inchi.MolToInchi(x) for x in mols]      # generate inchi keys from mols
    df['inchis'] = inchis                                  # add "inchis" column to pandas frame

    df = df.drop_duplicates(subset=['inchis'])             # drop duplicate rows based on "inchis" column
    num_current = df.shape[0]                              # check number of compounds after duplicate removal
    num_removed = num_original - num_current               # check how many were removed
    print(str(num_removed) + ' duplicate compounds were removed.')

    df = df.loc[:, df.columns != 'inchis']                 # remove the temporarily added column of "inchis"
    PandasTools.WriteSDF(df,filename[:-4]+'_DuplicatesRemoved.sdf',molColName='structure',properties=list(df.columns))  # export _DuplicatesRemoved sdf file

if __name__ == '__main__':
    main()
