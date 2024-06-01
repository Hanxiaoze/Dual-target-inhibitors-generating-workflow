import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import ConvertToNumpyArray
import os

from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'

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

# 读取分子数据并生成分子对象列表
def read_molecules(smiles_list):
    molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    return molecules

# 计算分子指纹数组
def compute_fingerprints(molecules):
    fingerprints = []
    for mol in molecules:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
        arr = np.zeros((1,))
        ConvertToNumpyArray(fp, arr)
        fingerprints.append(arr)
    return np.array(fingerprints)

# 绘制分子指纹热图
def plot_fingerprints(fingerprints):
    plt.matshow(fingerprints, cmap='viridis', aspect='auto')
    plt.xlabel('Bit Index')
    plt.ylabel('Molecule ID')
    plt.title('Molecular Fingerprints')
    plt.colorbar()
    plt.show()

# 读取分子数据
mol_list = read_mol2_files('./mol2_Q_to_smi_warm_up')

# # 示例分子 SMILES 列表
# smiles_list = ['CCO', 'CCN', 'CCC']

# # 读取分子数据
# molecules = read_molecules(smiles_list)

# 计算分子指纹数组
fingerprints = compute_fingerprints(mol_list)

# 绘制分子指纹热图
plot_fingerprints(fingerprints)