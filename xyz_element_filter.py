import os
import re

input_xyz_dir_path='./ZINC_in-man_xyz_with_prop'
filtered_xyz_dir_path='./ZINC_in-man_xyz_with_prop_elem_filter'

if os.path.exists(filtered_xyz_dir_path):
    os.system('rm -r {}'.format(filtered_xyz_dir_path))
os.makedirs(filtered_xyz_dir_path, exist_ok=True)


# 获取工作目录下所有xxx.xyz文件
xyz_files = [file for file in os.listdir(input_xyz_dir_path) if file.endswith('.xyz')]

# 循环处理每个mol2文件
for xyz_file in xyz_files:
    with open(str(input_xyz_dir_path)+'/'+str(xyz_file), 'r') as file:
        contents = file.read()

    # 使用正则表达式提取原子坐标
    atom_pattern = re.compile(r'([A-Z])(?:([a-z]))?\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s')
    atoms = atom_pattern.findall(contents)
    how_many_atoms = len(atoms)

    # 打印原子坐标
    print(f"原子坐标 ({xyz_file}):")
    flg = True
    for atom in atoms:
        element, t, x, y, z = atom
        print(f"{element}{t}: {x}, {y}, {z}")
        #if element+t not in ['H', 'B', 'C', 'N', 'O', 'F', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl']:
        #[1, 6, 7, 8, 9, 16, 17]  
        # if element+t not in ['H', 'C', 'N', 'O', 'F', 'S', 'Cl']:
        if element+t not in ['H', 'C', 'N', 'O', 'F']:
            flg = False
        if how_many_atoms > 55:
            flg = False

    print()

    if flg == True:
        # 写入xyz文件
        xyz_file_name=str(filtered_xyz_dir_path)+'/'+str(xyz_file)
        if os.path.exists(xyz_file_name):
            os.system('rm {}'.format(xyz_file_name))
        
        
        with open(xyz_file_name, "wt") as fout:
            fout.write(contents)
            

        print(f'Write 【{xyz_file_name}】 done !!!')

    else:
        continue

