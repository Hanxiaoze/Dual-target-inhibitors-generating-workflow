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
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 6), sharex=True)

# 绘制上面子图
axes[0].hist([big_E_float_list, big_M_float_list], bins=20, density=True, range=x_range, color=['blue', 'orange'], alpha=0.7, label=['Training Epro affinity distribution', 'Training Mpro affinity distribution'])
axes[0].legend()
axes[0].set_title('Affinity Distribution Probility Histogram for Training Epro and Training Mpro')
axes[0].set_xlim([-14, -2])  # 设置 x 轴范围
axes[0].set_xlabel('Affinity Value (kcal/mol)')  # 设置 x 轴标签
axes[0].set_ylabel('Probility')  # 设置 y 轴标签
axes[0].tick_params(axis='x', which='both', labelbottom=True)  # 显示 x 轴刻度标签

# 绘制下面子图
axes[1].hist([small_E_float_list, small_M_float_list], bins=20, density=True, range=x_range, color=['blue', 'orange'], alpha=0.7, label=['Generated Epro affinity distribution', 'Generated Mpro affinity distribution'])
axes[1].legend()
axes[1].set_title('Affinity Distribution Probility Histogram for Generated Epro and Generated Mpro')
axes[1].set_xlim([-14, -2])  # 设置 x 轴范围
axes[1].set_xlabel('Affinity Value (kcal/mol)')  # 设置 x 轴标签
axes[1].set_ylabel('Probility')  # 设置 y 轴标签


# 显示图形
plt.tight_layout()
plt.show()
