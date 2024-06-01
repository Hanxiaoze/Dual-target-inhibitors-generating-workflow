import numpy as np
import sys

import argparse

def get_parser():
    main_parser = argparse.ArgumentParser()
    main_parser.add_argument('data_path',
                             help='Path to .npz file')
    return main_parser

parser = get_parser()
args = parser.parse_args()
data_path = args.data_path

# 更改输出的最大宽度
sys.setrecursionlimit(10000)
np.set_printoptions(threshold=np.inf, linewidth=300)

# 读取npz文件
data = np.load(data_path)
print(data.files)
ii = 0
for i in data['stat_heads']:    
    print(str(i)+": "+str(data['stats'][ii]))
    ii += 1

stat_heads=data['stat_heads']
stats=data['stats']

for s in range(len(stats[0].tolist())):
    print(f'The statistic result of 【 molecule {s+1} 】:')
    line =""
    for t in range(len(stat_heads)):
        line = line + str(stat_heads[t]) + ": " + str(int(stats[t][s])) + "   "
    print(line+"\n\n")
    

data.close()
