import pandas as pd
from rdkit import Chem
from tkinter.filedialog import askopenfilename

def main():
    #     test_smi = 'C:/ProgramData/Simulations Plus, Inc/ADMET_Predictor9.5/Demo2D.smi'
    filename = askopenfilename()
    #     df = pd.read_csv(test_smi, sep ='\t')
    df = pd.read_csv(filename, sep = '\t')

    num_original = df.shape[0]  

    mols = [Chem.MolFromSmiles(x) for x in df['Canonical SMILES'].tolist()]   # get mols

    inchi = []
    for i in mols:
        if i != None:
            inchi.append(Chem.inchi.MolToInchi(i))
        else:
            inchi.append(None)

    df['inchis'] = inchi
    df[df.duplicated(subset=['inchis'], keep=False)].to_csv(filename[:-4]+'_Duplicates.smi',sep = '\t', index=False)
    df = df.drop_duplicates(subset=['inchis'], keep='last')

    num_current = df.shape[0]                              # check number of compounds after duplicate removal
    num_removed = num_original - num_current               # check how many were removed
    print(str(num_removed) + ' duplicate compounds were removed.')

    df = df.loc[:, df.columns != 'inchis']                 # remove the temporarily added column of "inchis"
    df.to_csv(filename[:-4]+'_DuplicatesRemoved.smi',sep = '\t', index=False)  # export _DuplicatesRemoved sdf file)

if __name__ == '__main__':
    main()
