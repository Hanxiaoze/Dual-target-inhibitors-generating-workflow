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
    data['and_'+str(i)] = (-data['Value1'] > i+1) & (-data['Value2'] > i)            #计算训练集与门列 
    data_0['and_'+str(i)] = (-data_0['Value1'] > i+1) & (-data_0['Value2'] > i)      #计算生成集与门列 

    prob_list.append((data['and_'+str(i)].mean()))
    prob_list_0.append((data_0['and_'+str(i)].mean()))

plt.figure(figsize=(12, 4))
plt.plot(np.arange(0.5, 8.5, 0.5), prob_list, marker='o', color='Blue', label='Training molecule set')
plt.plot(np.arange(0.5, 8.5, 0.5), prob_list_0, marker='x', color='Orange', label='Generated molecule set')

plt.xlabel('Param i of  \'AND gate criterion\' (\u2212Epro>i+1 & \u2212Mpro>i) of the two binding affinity values')
plt.ylabel('Probability')
plt.xticks(np.arange(0, 9.0, 0.5))
plt.legend()
plt.show()