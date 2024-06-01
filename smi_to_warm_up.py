import os
import shutil
import re
import sys
import math

arg1 = float(sys.argv[1])

def get_numeric_value(filename):
    match = re.search(r'x_(\d+)\.smi', filename)
    if match:
        return int(match.group(1))
    return None

def copy_files_above_median(source_dir, target_dir):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # 获取所有符合格式的文件名
    files = [f for f in os.listdir(source_dir) if re.match(r'x_\d+\.smi', f)]
    
    # 按照文件名中的数字部分排序
    files.sort(key=get_numeric_value)

    # 计算中间索引值
    mid_index = math.floor(len(files) * arg1)

    # 复制数字大于中间值的文件到目标文件夹
    for file in files[mid_index:]:
        shutil.copy(os.path.join(source_dir, file), target_dir)
        print(f"Copied {file} to {target_dir}")

# 示例调用
source_directory = './mol2_Q_to_smi'
target_directory = './mol2_Q_to_smi_warm_up'
copy_files_above_median(source_directory, target_directory)