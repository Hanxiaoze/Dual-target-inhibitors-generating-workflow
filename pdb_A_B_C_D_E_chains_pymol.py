import os
# 打开原始文件和新文件
with open('./target_protein_1.pdb', 'r') as input_file, open('./target_protein_1.pdb_tmp', 'w') as output_file:
    # 逐行读取原始文件
    lines = input_file.readlines()

    # 遍历每一行
    for i, line in enumerate(lines):
        # 检查是否在 第一条链 范围内
        if 2 <= i + 1 <= 1216:                            
            if len(line) >= 22:
                line = line[:21] + "A" + line[22:]

        # 检查是否在 第二条链 范围内
        if 1217 <= i + 1 <= 2431:                            
            if len(line) >= 22:
                line = line[:21] + "B" + line[22:]

        # 检查是否在 第三条链 范围内
        if 2432 <= i + 1 <= 3646:                            
            if len(line) >= 22:
                line = line[:21] + "C" + line[22:]

        # 检查是否在 第四条链 范围内
        if 3647 <= i + 1 <= 4861:                            
            if len(line) >= 22:
                line = line[:21] + "D" + line[22:]

        # 检查是否在 第五条链 范围内
        if 4862 <= i + 1 <= 6076:                            
            if len(line) >= 22:
                line = line[:21] + "E" + line[22:]

        # 写入修改后的行到新文件
        output_file.write(line)

os.system(f"mv ./target_protein_1.pdb_tmp ./target_protein_1.pdb")