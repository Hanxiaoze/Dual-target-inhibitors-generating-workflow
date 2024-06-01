import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger
from sklearn.model_selection import train_test_split
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.manifold import TSNE
# from mlinsights.mlmodel import PredictableTSNE

from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'

from hdbscan import HDBSCAN
 
sns.set_context('poster')
sns.set_style('white')
sns.set_color_codes()
plot_kwds = {'alpha' : 0.5, 's' : 80, 'linewidths':0}
RDLogger.DisableLog('rdApp.*')
seed = 794

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

fps = []


# for mol in mol_list:
#     fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)      # 2 是 Morgan 指纹的半径，它决定了在生成指纹时要考虑的邻居范围的大小。较大的半径会考虑更远的邻居环境。
#     arr = np.zeros((1,))   # 初始化一个 1 维的 0 数组          # fp 是一个稀疏的二进制向量表示生成的 Morgan 指纹
#     DataStructs.ConvertToNumpyArray(fp, arr)
#     flag = True
#     for a in fps:
#         if np.array_equal(a, arr):
#             flag = False
#             break
#     if flag:
#         fps.append(arr)        


for mol in mol_list:
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)      # 2 是 Morgan 指纹的半径，它决定了在生成指纹时要考虑的邻居范围的大小。较大的半径会考虑更远的邻居环境。
    arr = np.zeros((1,))   # 初始化一个 1 维的 0 数组          # fp 是一个稀疏的二进制向量表示生成的 Morgan 指纹
    DataStructs.ConvertToNumpyArray(fp, arr)
    fps.append(arr)        


fps = np.array(fps)
print(fps.shape)

# clusterer = HDBSCAN(algorithm='generic', 
#                     min_cluster_size=7,
#                     # allow_single_cluster=True,
#                     min_samples=3, 
#                     metric='euclidean',
#                     # cluster_selection_method='leaf'
#                     )


# 实例化HDBSCAN类 
clusterer = HDBSCAN(algorithm='generic', 
                    min_cluster_size=6,
                    # allow_single_cluster=True,
                    min_samples=1, 
                    metric='euclidean',
                    # cluster_selection_method='leaf'
                    )
clusterer.fit(fps)    # 聚类

palette = sns.color_palette()


cluster_colors = [sns.desaturate(palette[min(col, len(palette) - 1)] if col >= 0 else (0.5, 0.5, 0.5), sat)
                  for col, sat in zip(clusterer.labels_, clusterer.probabilities_)]
 

# tsne = TSNE(random_state=seed)
tsne = TSNE(n_components=3, random_state=seed)
res = tsne.fit_transform(fps)      # 降维



# plt.scatter(res[:,0], res[:,1], c=cluster_colors, s=800, alpha=0.2, linewidths=0)     # 绘制二维



# 创建 3D 坐标轴
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 绘制散点图
ax.scatter(res[:,0], res[:,1], res[:,2], c=cluster_colors, s=100, alpha=1, linewidths=0)   # 绘制三维
ax.set_xlabel('t-SNE Component 1', fontsize=12)
ax.set_ylabel('t-SNE Component 2', fontsize=12)
ax.set_zlabel('t-SNE Component 3', fontsize=12)
# 设置坐标轴刻度的字体大小
ax.tick_params(axis='both', which='major', labelsize=12)

plt.show()