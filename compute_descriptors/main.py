# -*- coding: utf-8 -*-
"""
Created on Wed May 20 22:03:52 2020

@author: phyophyokyawzin
"""


import time, os
import pandas as pd
from molvs import standardize_smiles
from rdkit import Chem
from MACCS import *
from RDKit_2D import *
from ECFP6 import *
from Macrocycle_Descriptors import *


def convert_time(second):
    day = second/86400
    hour = (day - int(day))*24
    minute = (hour - int(hour))*60
    second = round((minute - int(minute))*60,4)
    return(str(int(day)) + ' DAYS: '+ str(int(hour)) + ' HOURS: '+ str(int(minute)) + ' MINUTES: ' + str(second) + ' SECONDS')

def create_directory():
    while not os.path.exists('desc'):
        os.mkdir('desc')


def main():

    filename = 'data/macrolides_smiles.csv'
    create_directory()
    df = pd.read_csv(filename)
    smiles = [standardize_smiles(i) for i in df['smiles'].values]

    start_time = time.time()
    output_filename = 'desc' + filename[4:]

    ### Compute ECFP6 Fingerprints and export file.
    ecfps_descriptor = ECFP6(smiles)
    ecfps_descriptor.compute_ECFP6(output_filename)

    ## Compute MACCS Fingerprints and export file.
    maccs_descriptor = MACCS(smiles)
    maccs_descriptor.compute_MACCS(output_filename)

    ## Compute RDKit 2D Descriptors and export file.
    rdk_descriptor = RDKit_2D(smiles)
    rdk_descriptor.compute_2Drdkit(output_filename)

    ## Compute mordred_mrc Descriptors and export file.
    mrc_descriptor = Macrocycle_Descriptors(smiles)
    mrc_descriptor.mordred_compute(output_filename)
    mrc_descriptor.compute_mordred_macrocycle(output_filename)


    duration = convert_time(time.time()-start_time)
    print(duration)


if __name__ == '__main__':
    main()