import os
import re

xyz_no_property_path = '/home/zhou/sch/cG-SchNet-main/data/world_small_55/xyz_55_HCNOF_raw'

xyz_with_property_path = '/home/zhou/sch/cG-SchNet-main/data/world_small_55/xyz_55_HCNOF_with_property'

prop_1_path = '/home/zhou/Dock_working/world_E_alphafold_step7_1_1st_frame/in_site_last_information.txt'

prop_2_path = '/home/zhou/Dock_working/CoV2_Mpro_7ZB7/last_information.txt'

os.system('rm {}/*'.format(xyz_with_property_path))

# 获取工作目录下所有 xxx.xyz 文件
xyz_files = [file for file in os.listdir(xyz_no_property_path) if file.endswith('.xyz')]
# 循环处理每个 .xyz 文件
for no_prop_file in xyz_files:
    with open(str(xyz_no_property_path)+'/'+str(no_prop_file), 'r') as n_p_file:
        pppp = ''
        pp_1 = ''
        pp_2 = ''
        # 'world969.mol2.xyz'
        # 'world3303_log.txt   -9.8   ZINC000006716957'
        flg = str(no_prop_file).split('.')[0] + '_log.txt'
        with open(str(prop_1_path), 'r') as p1_file:
            p1_lines = p1_file.readlines()
            for n1 in range(len(p1_lines)):
                if flg in p1_lines[n1]:
                    pp_1 = p1_lines[n1].split()[1]
                
        with open(str(prop_2_path), 'r') as p2_file:
            p2_lines = p2_file.readlines()
            for n2 in range(len(p2_lines)):
                if flg in p2_lines[n2]:
                    pp_2 = p2_lines[n2].split()[1]
                
        pppp = pp_1 + '    ' + pp_2
        if (len(pppp.split()) == 2):
            n_p_lines = n_p_file.readlines()
            for i in range(len(n_p_lines)):
                if n_p_lines[i] == 'Properties=species:S:1:pos:R:3 pbc="F F F"\n':
                    n_p_lines[i] = str(pppp) + '\n'
            
            with open(str(xyz_with_property_path)+'/'+str(no_prop_file), 'w') as w_p_file:
                w_p_file.writelines(n_p_lines)
        else:
            continue





