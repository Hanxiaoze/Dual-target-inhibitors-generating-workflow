import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'

# 读取txt文件并解析数据
file_path = './sorted_sum_E_M_gen_affinity_warm_up.txt' 
data = {'x_labels': [], 'data2': [], 'data3': [], 'data4': []}

c = 0

with open(file_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        c += 1
        if c <= 42:
            items = line.strip().split()
            x_label = items[0].replace('_log.txt', '')
            data['x_labels'].append(x_label)
            data['data2'].append(-float(items[1]))
            data['data3'].append(-float(items[2]))
            data['data4'].append(-float(items[3]))

# 绘制折线图
            
plt.figure(figsize=(12, 4))

plt.plot(data['x_labels'], data['data2'], linestyle=':', marker='^', label='Epro + Mpro')
plt.plot(data['x_labels'], data['data3'], linestyle=':', marker='x', label='Epro')
plt.plot(data['x_labels'], data['data4'], linestyle=':', marker='o', label='Mpro')

plt.axhline(y=17.0, color='r', linestyle='--')

# 在每个数据点的正上方标记数值
for i in range(len(data['x_labels'])):
    plt.text(data['x_labels'][i], data['data2'][i], f"{data['data2'][i]}", ha='center', va='bottom', rotation=45)
    plt.text(data['x_labels'][i], data['data3'][i]+0.5, f"{data['data3'][i]}", ha='center', va='bottom', rotation=45)
    plt.text(data['x_labels'][i], data['data4'][i]-0.5, f"{data['data4'][i]}", ha='center', va='top', rotation=45)

# 设置图形标题和标签
# plt.title('Data from txt file')
plt.xlabel('Molecule ID')
plt.ylabel('\u2212 Binding affinity values (kcal/mol)')

plt.ylim(6, 22)  
plt.yticks(np.arange(6, 22, 0.5))

# 添加图例
plt.legend()


# 显示图形
plt.xticks(rotation=90)  # 旋转x轴标签
plt.tight_layout()  # 自动调整布局，防止标签被裁剪
plt.show()