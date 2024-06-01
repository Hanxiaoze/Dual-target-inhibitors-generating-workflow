import sys

arg1 = str(sys.argv[1])
arg2 = str(sys.argv[2])

import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'


big_file_path = './E_M_6W_affinity.txt'
small_file_path = './E_M_gen_affinity_warm_up.txt'

# 手动设置 x 和 y 的范围
x_range = (-14, -2)


# 用于存储浮点数的列表
big_E_float_list = []
big_M_float_list = []

small_E_float_list = []
small_M_float_list = []

# 逐行读取文件
with open(big_file_path, 'r') as file:
    for line in file:
        # 按空格分割每一行
        fields = line.split()
        if len(fields) >= 3:
            # 尝试将第二个字段转换为浮点数，并添加到列表中
            try:                
                big_E_float_list.append(float(fields[1]))
                big_M_float_list.append(float(fields[2]))
            except ValueError:
                print(f"Unable to convert '{fields[1]} or {fields[2]}' to float in line: {line}")


with open(small_file_path, 'r') as small_file:
    for s_line in small_file:
        # 按空格分割每一行
        s_fields = s_line.split()
        if len(s_fields) >= 3:
            # 尝试将第二个字段转换为浮点数，并添加到列表中
            try:                
                small_E_float_list.append(float(s_fields[1]))
                small_M_float_list.append(float(s_fields[2]))
            except ValueError:
                print(f"Unable to convert '{s_fields[1]} or {s_fields[2]}' to float in line: {s_line}")


# 创建图形和子图
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))

# 计算概率分布
bins = np.linspace(-14, -2, 21)  # 指定 bin 的个数或宽度

hist_big_E, edges = np.histogram(big_E_float_list, bins=bins, range=x_range, density=True)
hist_big_M, _ = np.histogram(big_M_float_list, bins=bins, range=x_range, density=True)
hist_small_E, _ = np.histogram(small_E_float_list, bins=bins, range=x_range, density=True)
hist_small_M, _ = np.histogram(small_M_float_list, bins=bins, range=x_range, density=True)

# 计算折线图的中心点坐标
x_centers = (edges[:-1] + edges[1:]) / 2

# 绘制折线图
axes.plot(x_centers, hist_big_E, color='blue', marker='x', linestyle=':', label='Training molecules Epro binding affinity distribution')
axes.plot(x_centers, hist_big_M, color='orange', marker='x', linestyle=':', label='Training molecules Mpro binding affinity distribution')
axes.plot(x_centers, hist_small_E, color='blue', marker='o', linestyle='-', label='Generated molecules Epro binding affinity distribution')
axes.plot(x_centers, hist_small_M, color='orange', marker='o', linestyle='-', label='Generated molecules Mpro binding affinity distribution')

axes.legend()
# axes.set_title('Generated molecules with affinity conditions Epro= -{}.0 / Mpro= -{}.0 kcal/mol'.format(arg1, arg2))
axes.set_xlim([-14, -2])  # 设置 x 轴范围
axes.set_xlabel('Binding affinity value (kcal/mol)')  # 设置 x 轴标签
axes.set_ylabel('Probability density')  # 设置 y 轴标签

axes.set_ylim(0, 0.8)

# 显示图形
plt.tight_layout()
plt.show()
