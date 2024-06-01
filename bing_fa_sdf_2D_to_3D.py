import os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ProcessPoolExecutor
def generate_3d_sdf(input_sdf, output_sdf):
    suppl = Chem.SDMolSupplier(input_sdf)
    writer = Chem.SDWriter(output_sdf)
    for mol in suppl:
        if mol is not None:
            mol = Chem.AddHs(mol)   # 添加氢原子
            params = AllChem.ETKDGv3()
            params.useSmallRingTorsions = True
            AllChem.EmbedMultipleConfs(mol, numConfs = 1 , params = params)
            #AllChem.EmbedMolecule(mol, AllChem.ETKDGv2(), randomSeed=42)   # 采样生成 3D 结构
            writer.write(mol)
    writer.close()
if __name__ == "__main__":
    input_folder = '/home/zhou/ZINC_data_bank/ZINC_in-man_sdf_files_222D'  # 输入文件夹路径
    output_folder = '/home/zhou/ZINC_data_bank/ZINC_in-man_sdf_files_222D_to_333D'  # 输出文件夹路径
    num_workers = 72  # 设置并行工作核的数量
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for file_name in os.listdir(input_folder):
            if file_name.endswith('.sdf'):
                input_sdf_file = os.path.join(input_folder, file_name)
                output_sdf_file = os.path.join(output_folder, f'{file_name}')
                executor.submit(generate_3d_sdf, input_sdf_file, output_sdf_file)
    print("3D SDF文件已生成。")
