import matplotlib.pyplot as plt
import numpy as np


# 从文件中读取数据
file_path = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/E_M_sum_affinity.out'
data = np.loadtxt(file_path, skiprows=0, usecols=(1, 2))

# 创建3D图形
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 使用scatter()函数绘制3D散点图
hist = ax.scatter(data[:, 0], data[:, 1], c=data[:, 1], cmap='Blues', marker='o', s=5, alpha=0.5, zdir='z', depthshade=True)

# 添加标题和标签
ax.set_title('3D Frequency Histogram')
ax.set_xlabel('Epro_affinity_value')
ax.set_ylabel('Mpro_affinity_value')
ax.set_zlabel('Frequency')

# 添加颜色条
fig.colorbar(hist, ax=ax, label='Frequency')

# 限制Z轴范围在 0-0.2
ax.set_zlim(0, 0.002)

# 显示图形
plt.show()
