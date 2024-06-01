# 读取原始文件
input_file_path = './E_M_gen_affinity_warm_up.txt'
output_file_path = './sorted_sum_E_M_gen_affinity_warm_up.txt'

# 读取数据并计算和
data = []
with open(input_file_path, 'r') as file:
    for line in file:
        parts = line.split()
        if len(parts) >= 3:
            name = parts[0]
            value1 = float(parts[1])
            value2 = float(parts[2])
            total = value1 + value2
            data.append((name, total, value1, value2))

# 按和排序
sorted_data = sorted(data, key=lambda x: float(x[1]))

# 写入新文件
with open(output_file_path, 'w') as file:
    for entry in sorted_data:
        file.write(f"{entry[0]}  {entry[1]:.2f}  {entry[2]}  {entry[3]}\n")

print("Process completed.")
