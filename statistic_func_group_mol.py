import matplotlib.pyplot as plt
from collections import defaultdict
from rdkit import Chem
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

# 读取分子数据
mol_list = read_mol2_files('./mol2_Q_to_smi_warm_up')

# 定义要统计的官能团
functional_groups = {
    '-O-': 'O',
    '-OH': '[OH]',
    '-N=': 'N',
    '-COOH': 'C(=O)[OH]',    
    '-NH2': '[NH2]',
    'Secondary Amine group': '[NH]',
    'Tertiary Amine group': 'CN(C)C',
    'Benzene Ring': 'c1ccccc1',
    '-CHO': '[CH]=O',
    '-(C=O)OC-': 'C(=O)OC',
    '-(C=O)-': 'C(=O)',
    '-(C=O)N-': 'C(=O)N',
    '-C(C=O)C-': 'C(C=O)C',
    'Pyridine group': 'n1ccccc1',
    'Pyrrole group': 'n1cccc1',
    'Furan group': 'c1ccoc1',
    'Imidazole group': 'c1cnc[nH]1',
    '-C=C-': 'C=C',
    '-C≡C-': 'C#C',
    '-C=N-': 'C=N',
    '-F': 'F',
    '-Cl': 'Cl',
    '-C≡N': 'C#N'
}

# 初始化官能团计数器
functional_group_counts = defaultdict(int)

# 遍历SMILES分子数据集
for mol in mol_list:

    if mol is not None:
        # 遍历每种官能团
        for fg_name, pattern in functional_groups.items():
            substruct = Chem.MolFromSmarts(pattern)
            if substruct is not None:
                # 检查分子中是否存在该官能团
                if mol.HasSubstructMatch(substruct):
                    functional_group_counts[fg_name] += 1

# 按照 functional_group_counts 中的顺序对键值对进行排序
sorted_functional_groups = sorted(functional_group_counts.items(), key=lambda x: list(functional_groups.keys()).index(x[0]))

# 提取排序后的键和值
sorted_keys = [item[0] for item in sorted_functional_groups]
sorted_values = [item[1] for item in sorted_functional_groups]

# 绘制柱状图
bars = plt.bar(sorted_keys, sorted_values, color='blue')

# 添加数值标签
for bar, count in zip(bars, sorted_values):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1, str(count), ha='center', color='black')
plt.xlabel('Functional Group')
plt.ylabel('Count molecules')
# plt.title('Functional Group Counts')
plt.xticks(rotation=45, ha='right')
plt.show()