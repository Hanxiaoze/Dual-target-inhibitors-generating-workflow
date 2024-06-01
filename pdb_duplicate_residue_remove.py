import os

# 打开原始文件和新文件
with open('./target_protein_2.pdb', 'r') as input_file, open('./target_protein_2.pdb_tmp', 'w') as output_file:
    # 逐行读取原始文件
    lines = input_file.readlines()

    # 遍历每一行
    for i, line in enumerate(lines):
        if len(line) >= 17:        
            if line[16] != "A" and line[16] != "B":
                output_file.write(line)
            if line[16] == "A":
                line = line[0 : 16] + " " + line[17 : ]
                output_file.write(line)
        else:
            output_file.write(line)

os.system(f"mv ./target_protein_2.pdb_tmp ./target_protein_2.pdb")