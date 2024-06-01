import os
import re

mol2_dir_path='/home/zhou/Dock_working/5J89'
xyz_dir_path='/home/zhou/sch/cG-SchNet-main/data/world_small_55/xyz_55_HCNOF_raw'

os.system('rm {}/*'.format(xyz_dir_path))

# 获取工作目录下所有xxx.mol2文件
mol2_files = [file for file in os.listdir(mol2_dir_path) if file.endswith('.mol2')]

# 循环处理每个mol2文件
for mol2_file in mol2_files:
    with open(str(mol2_dir_path)+'/'+str(mol2_file), 'r') as file:
        contents = file.read()

    # 使用正则表达式提取原子坐标
    atom_pattern = re.compile(r'(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+([A-Z])(?:([a-z]))?')
    atoms = atom_pattern.findall(contents)
    how_many_atoms = len(atoms)

    # 打印原子坐标
    print(f"原子坐标 ({mol2_file}):")
    flg = True
    for atom in atoms:
        x, y, z, element, t = atom
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
        xyz_file_name=str(xyz_dir_path)+'/'+str(mol2_file)+'.xyz'
        os.system('rm {}'.format(xyz_file_name))
        
        
        with open(xyz_file_name, "wt") as fout:
            fout.write(str(how_many_atoms)+'\n')
            fout.write('Properties=species:S:1:pos:R:3 pbc="F F F"\n')
            for atom in atoms:
                x, y, z, element, t = atom
                fout.write('{}{}        {}         {}         {}\n'.format(element, t, x, y, z))

        print(f'Write 【{xyz_file_name}】 done !!!')

    else:
        continue

