import os
# 设置输入pdb文件夹和输出xyz文件夹
#input_folder = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/generated_1_xyz
#output_folder = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/generated_1_pdb'
input_folder = '.'
output_folder = '.'

# 获取输入文件夹中的所有pdb文件
input_files = [os.path.join(input_folder, filename) for filename in os.listdir(input_folder) if filename.endswith('.pdb')]


for input_file in input_files:
    # 打开pdb文件进行读取
    with open(input_file, 'r') as pdb_file:
        lines = pdb_file.readlines()

    # 创建XYZ文件进行写入
    output_file = os.path.join(output_folder, os.path.basename(input_file).replace('.pdb', f'_re.xyz'))
    with open(output_file, 'w') as xyz_file:
        # 写入XYZ头部信息
        xyz_file.write(str(len(lines) - 2) + "\n")
        xyz_file.write("\n")

        atom_number = 1

        # 写入原子信息
        for line in lines[2:]:  # 跳过PDB文件的前两行
            atom_data = line.split()
            
            xyz_file.write(f"{atom_data[2]:4} {float(atom_data[5]):10.4f} {float(atom_data[6]):10.4f} {float(atom_data[7]):10.4f}\n")