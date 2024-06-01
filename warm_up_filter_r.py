import sys
import math


arg1 = float(str(sys.argv[1]))


input_file = './E_M_gen_affinity.txt'
output_file = './E_M_gen_affinity_warm_up.txt'

rate = arg1

with open(input_file, 'r') as inf:
    lines = inf.readlines()                 # 读取文件内容到列表中
    total_num = len(lines)

    with open(output_file, 'w') as outf:
        for infl in lines:                  # 遍历列表中的每一行内容
            id = int(infl.split('_')[1])

            if id > math.floor(total_num * rate):
                outf.write(infl)