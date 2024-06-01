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
# data_path = '/home/zhou/sch/cG-SchNet-main/splits/4_comp_relenergy_split.npz'
# ['train_idx', 'val_idx', 'test_idx', 'invalid_idx']

# 更改输出的最大宽度
sys.setrecursionlimit(10000)
np.set_printoptions(threshold=np.inf, linewidth=300)

# 读取npz文件
data = np.load(data_path)
print(data.files)

for i in data['train_idx']: 
    print('train_idx:  ' + str(i))

for i in data['val_idx']: 
    print('val_idx:  ' + str(i))

for i in data['test_idx']: 
    print('test_idx:  ' + str(i))

for i in data['invalid_idx']: 
    print('invalid_idx:  ' + str(i))



data.close()