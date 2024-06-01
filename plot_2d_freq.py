import matplotlib.pyplot as plt
import numpy as np

# 从文件中读取数据
file_path = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/E_M_6W_affinity.txt'
data = np.loadtxt(file_path, skiprows=0, usecols=(1, 2))

# 使用plt.hist2d()函数绘制二维频率统计直方图
plt.hist2d(data[:, 0], data[:, 1], bins=40, cmap='Blues')

# 添加标题和标签
plt.title('2D Frequency Histogram')
plt.xlabel('Epro_affinity_value')
plt.ylabel('Mpro_affinity_value')

# 添加颜色条
plt.colorbar(label='Frequency')

# 显示图形
plt.show()