import os 

mol_save = []
save_dir_path = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/sorted_E_M_gen_affinity_new.txt'
remove_dir_path = '/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/generated_1_mol2_Q_to_smi'
with open(save_dir_path) as f:
    for l in f.readlines():
        mol_save.append(l.split('_log.txt')[0])

print(len(mol_save))
print(mol_save)


input_folder = remove_dir_path
input_files = [os.path.join(input_folder, filename) for filename in os.listdir(input_folder) if filename.endswith('.smi')]


for to_remove in input_files:
    if to_remove.split('.')[0].split('/')[-1] not in mol_save:
        os.system(f'rm {to_remove}')
        print(f'rm {to_remove}')
