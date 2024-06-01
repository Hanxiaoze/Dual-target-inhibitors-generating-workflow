import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import ConvertToNumpyArray
import os

input_folder ='./mol2_Q_to_smi_warm_up'
input_files = [os.path.join(input_folder, filename) for filename in os.listdir(input_folder) if filename.endswith('.smi')]

for i in input_files:
    print(i)
    f = open(i)
    mol = f.readlines()[0].split()[0]
    try:
        mol = Chem.MolFromSmiles(mol)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    except:
        print(i, '#########')



def read_mol2_files(folder):
    """读取指定文件夹中所有的smi文件，并返回分子列表"""
    mol_list = []
    for filename in os.listdir(folder):
        if filename.endswith(".smi"):
            filepath = os.path.join(folder, filename)
            f = open(filepath)
            mol = f.readlines()[0].split()[0]
            mol = Chem.MolFromSmiles(mol)
            if mol is not None:
                mol_list.append(mol)
    return mol_list



# 计算分子指纹数组
def compute_fingerprints(molecules):
    fingerprints = []
    for mol in molecules:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
        arr = np.zeros((1,))
        ConvertToNumpyArray(fp, arr)
        fingerprints.append(arr)
    return np.array(fingerprints)