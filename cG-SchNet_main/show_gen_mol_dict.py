from ase import Atoms
import ase.visualize as asv
import pickle

import argparse
def get_parser():
    main_parser = argparse.ArgumentParser()
    main_parser.add_argument('data_path',
                             help='Path to .npz file')
    return main_parser

parser = get_parser()
args = parser.parse_args()
data_path = args.data_path

# 打开二进制文件
file = open(data_path, "rb")

# 读取序列化对象
generated = pickle.load(file)
file.close()

print(generated)

ats = []
n_total_atoms = 0
n_molecules = 0
for key in generated:
    n = 0
    for i in range(len(generated[key]['_atomic_numbers'])):
        at = Atoms(generated[key]['_atomic_numbers'][i],
        positions=generated[key]['_positions'][i])
        ats += [at]
        n += 1
        n_molecules += 1
        n_total_atoms += n * key
asv.view(ats)
print(f'Total number of atoms placed: {n_total_atoms} '
    f'(avg {n_total_atoms / n_molecules:.2f})', flush=True)
