import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'

# 读取文件

file_path ='./E_M_6W_affinity.txt'
file_path_0 = './E_M_gen_affinity_warm_up.txt'  

data = pd.read_csv(file_path, delim_whitespace=True, header=None, names=['Name', 'Value1', 'Value2'])
data_0 = pd.read_csv(file_path_0, delim_whitespace=True, header=None, names=['Name', 'Value1', 'Value2'])


prob_list = []
prob_list_0 = []
for i in np.arange(0.5, 8.5, 0.5):
    data['and_'+str(i)] = (-data['Value1'] > i+0) & (-data['Value2'] > i)            #计算训练集与门列 
    data_0['and_'+str(i)] = (-data_0['Value1'] > i+0) & (-data_0['Value2'] > i)      #计算生成集与门列 

    prob_list.append((data['and_'+str(i)].mean()))
    prob_list_0.append((data_0['and_'+str(i)].mean()))

plt.figure(figsize=(12, 4))
plt.plot(np.arange(0.5, 8.5, 0.5), prob_list, marker='o', color='Blue', label='Training molecule set')
plt.plot(np.arange(0.5, 8.5, 0.5), prob_list_0, marker='x', color='Orange', label='Generated molecule set')
# 标记出每个数据点的y值
# for i in range(len(prob_list)):
#     plt.text(np.arange(1, 11, 0.5)[i], prob_list[i], f'{round(prob_list[i],6)}', color='Blue', ha='right', va='bottom')
#     plt.text(np.arange(1, 11, 0.5)[i], prob_list_0[i]+0.05, f'{round(prob_list_0[i],6)}', color='Orange', ha='right', va='bottom')

plt.xlabel('Param i of  \'AND gate criterion\' (\u2212Epro>i & \u2212Mpro>i) of the two binding affinity values')
plt.ylabel('Probability')
plt.xticks(np.arange(0, 8.5, 0.5))
plt.legend()
plt.show()