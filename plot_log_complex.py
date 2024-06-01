import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'

# 创建一个Pandas DataFrame来加载数据
data = pd.read_csv('./log.csv')

x = (data.index + 1).tolist()

# 提取各列数据
y1 = data['Time']
y2 = data['Learning rate']
y3 = data['Train loss']
y4 = data['Validation loss']
y5 = data['KLD_types']
y6 = data['KLD_dists']

# 创建画布
fig, ax1 = plt.subplots(figsize=(10, 6))  # 调整图的大小

# 绘制y1
#plt.plot(x, y1, label='Training Time')


# 绘制y3
ax1.plot(x, y3, label='Training loss', linewidth=3.5, alpha=0.5)

# 绘制y4
ax1.plot(x, y4, label='Validation loss', linewidth=3.5, alpha=0.5)

# 绘制y5
ax1.plot(x, y5, label='KLD_types')

# 绘制y6
ax1.plot(x, y6, label='KLD_dists')

ax1.set_xlabel('Training epochs', fontsize=12)
ax1.set_ylabel('Loss metric', fontsize=12)

#ax1.tick_params(axis='y', colors='blue')
#ax1.yaxis.set_tick_params(color='blue')
#ax1.spines['left'].set_edgecolor('blue')


# 添加图例
legend1 = ax1.legend(loc='upper right', fontsize=12)
ax1.add_artist(legend1)
# plt.legend()



############
############
# 绘制y2
ax2 = ax1.twinx()  # 创建次坐标轴
ax2.plot(x, y2, label='Learning rate', color='mediumvioletred', linewidth=1.5, alpha=0.2)
ax2.set_ylabel('Learning rate', color='mediumvioletred', fontsize=12)
ax2.tick_params(axis='y', colors='mediumvioletred')
ax2.yaxis.set_tick_params(color='mediumvioletred')
ax2.spines['right'].set_edgecolor('mediumvioletred')


legend2 = ax2.legend(loc='upper left', fontsize=12)



# 设置标题
# plt.title('Training Results', fontsize=14)

# 自适应y轴刻度范围
plt.autoscale(enable=True, axis='y')

# 显示图形
plt.show()
