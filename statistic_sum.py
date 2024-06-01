import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'

# 读取文件
file_path_0 = './E_M_gen_affinity_warm_up.txt'  # 请替换为实际文件路径
file_path ='./E_M_6W_affinity.txt'

data = pd.read_csv(file_path, delim_whitespace=True, header=None, names=['Name', 'Value1', 'Value2'])
data_0 = pd.read_csv(file_path_0, delim_whitespace=True, header=None, names=['Name', 'Value1', 'Value2'])

# 计算后两列之和
data['Sum'] = -data['Value1'] - data['Value2']
data_0['Sum'] = -data_0['Value1'] - data_0['Value2']

# 计算大于1, 2, 3, 4的概率
thresholds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
probabilities = [(data['Sum'] > threshold).mean() for threshold in thresholds]
probabilities_0 = [(data_0['Sum'] > threshold).mean() for threshold in thresholds]

# 绘制折线图
plt.figure(figsize=(12, 4))
plt.plot(thresholds, probabilities, marker='o', color='Blue', label='Training molecule set')
plt.plot(thresholds, probabilities_0, marker='x', color='Orange', label='Generated molecule set')

# 标记出每个数据点的y值
# for i in range(len(probabilities)):
#     plt.text(thresholds[i], probabilities[i], f'{round(probabilities[i],6)}', color='Blue', ha='center', va='bottom')
#     plt.text(thresholds[i], probabilities_0[i]+0.05, f'{round(probabilities_0[i],6)}', color='Orange', ha='center', va='bottom')

# plt.title('Probability of Sum > Threshold')
plt.xlabel('Param i of  \'SUM gate criterion\' ((\u2212Epro) + (\u2212Mpro) > i) of the two binding affinity values')
plt.ylabel('Probability')
plt.xticks(np.arange(1, 19, 1))
plt.legend()
plt.show()