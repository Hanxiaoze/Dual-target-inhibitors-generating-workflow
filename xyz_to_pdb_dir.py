import os
# 设置输入xyz文件夹和输出pdb文件夹
input_folder = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/generated_1_xyz'
output_folder = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/generated_1_pdb'
# input_folder = '.'
# output_folder = '.'

# 获取输入文件夹中的所有xyz文件
input_files = [os.path.join(input_folder, filename) for filename in os.listdir(input_folder) if filename.endswith('.xyz')]


for input_file in input_files:
    # 打开XYZ文件进行读取
    with open(input_file, 'r') as xyz_file:
        lines = xyz_file.readlines()

    # 创建PDB文件进行写入
    output_file = os.path.join(output_folder, os.path.basename(input_file).replace('.xyz', f'.pdb'))
    with open(output_file, 'w') as pdb_file:
        # 写入PDB头部信息
        pdb_file.write("REMARK XYZ to PDB conversion\n")
        pdb_file.write("TITLE Converted from XYZ file\n")

        atom_number = 1

        # 写入原子信息
        for line in lines[2:]:  # 跳过XYZ文件的前两行
            atom_data = line.split()
            element = atom_data[0]
            x, y, z = float(atom_data[1]), float(atom_data[2]), float(atom_data[3])
            pdb_file.write(f"ATOM  {atom_number:4} {element:3}  MOL     1    {x:8.3f}{y:8.3f}{z:8.3f}\n")
            atom_number += 1
