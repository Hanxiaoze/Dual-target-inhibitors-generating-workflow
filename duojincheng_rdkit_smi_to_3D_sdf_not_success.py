import os
from rdkit import Chem
from rdkit.Chem import AllChem
from multiprocessing import Pool

# 输入文件夹和输出文件夹
input_folder = '/home/zhou/ZINC_data_bank/zinc_in-man_222D_smi_files'  # 替换为包含SMILES文件的文件夹路径
output_folder = '/home/zhou/ZINC_data_bank/ZINC_in-man_smi_to_3D_sdf_files'  # 替换为保存SDF文件的文件夹路径

# 获取输入文件夹中的所有SMILES文件
input_files = [os.path.join(input_folder, filename) for filename in os.listdir(input_folder) if filename.endswith('.smi')]

def process_file(input_file):
    # 读取SMILES文件
    with open(input_file, 'r') as f:
        smiles_list = f.read().splitlines()

    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mol = Chem.AddHs(mol)
                #AllChem.EmbedMolecule(mol, randomSeed=42)  # 这样还是只能变成 2D sdf
                
                #params = AllChem.ETKDGv3()      # 这种方法也不行
                #AllChem.EmbedMultipleConfs(mol, numConfs = 1 , params = params)
                
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())

                output_file = os.path.join(output_folder, os.path.basename(input_file).replace('.smi', f'.sdf'))
                writer = Chem.SDWriter(output_file)
                writer.write(mol)
        except Exception as e:
            print(f"Error processing SMILES: {smiles} - {str(e)}")

if __name__ == '__main__':
    # 创建输出文件夹
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 创建进程池并并行处理文件
    with Pool(processes=os.cpu_count()) as pool:
        pool.map(process_file, input_files)

    print("Conversion completed.")
