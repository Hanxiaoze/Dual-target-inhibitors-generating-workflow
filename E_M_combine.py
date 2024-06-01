# E_file_path = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/Epro_sum_affinity.out'

# M_file_path = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/Mpro_sum_affinity.out'

# E_M_path_to_write = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/E_M_sum_affinity.out'



E_file_path = './Docking_dir_for_target_protein_1/sorted_final_energy.txt'

M_file_path = './Docking_dir_for_target_protein_2/sorted_final_energy.txt'

E_M_path_to_write = './E_M_gen_affinity.txt'


with open(M_file_path, 'r') as M_data:
    M_dict = {line.split()[0]: line.split()[1] for line in M_data}   #  一次性读取 M_file 加载成字典，提高查找性能

with open(E_file_path, 'r') as E_data, open(E_M_path_to_write, 'w') as E_M_file:  # 打开文件，'w'表示写入模式，如果文件存在则覆盖
    E_lines = E_data.readlines()   # 一次性将 E_file 读取进内存，避免循环中反复读取
    for E_line in E_lines:
        E_fields = E_line.split()
        mol_name = E_fields[0]
        if mol_name in M_dict:
            E_M_file.write(f"{mol_name}   {E_fields[1]}    {M_dict[mol_name]}\n")

