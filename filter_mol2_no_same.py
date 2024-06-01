import os
from rdkit import Chem
from shutil import copyfile

def mol2_to_smiles(file_path):
    """
    将mol2文件中的分子结构转换为SMILES表示
    """
    mol_supplier = Chem.MolFromMol2File(file_path, sanitize=False)
    if mol_supplier is not None:
        return Chem.MolToSmiles(mol_supplier, isomericSmiles=True)

def main():
    input_folder = "./generated_2_10kcal_10kcal/mol2_Q"
    output_folder = "./generated_2_10kcal_10kcal/mol2_Q_no_same"

    # 创建输出文件夹
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 存储SMILES和文件路径的字典
    smiles_dict = {}

    # 遍历文件夹中的所有.mol2文件
    for file_name in os.listdir(input_folder):
        if file_name.endswith('.mol2'):
            input_file_path = os.path.join(input_folder, file_name)

            # 将.mol2文件中的分子结构转换为SMILES
            smiles = mol2_to_smiles(input_file_path)

            if smiles in smiles_dict:
                print(os.path.basename(input_file_path), ' is same as ',os.path.basename(smiles_dict[smiles]))
            else:
                # 将SMILES和文件路径添加到字典中
                smiles_dict[smiles] = input_file_path

    # 检查SMILES是否有重复，如果有，则复制不同的.mol2文件到输出文件夹
    for smiles, file_path in smiles_dict.items():        
        output_file_path = os.path.join(output_folder, os.path.basename(file_path))
        copyfile(file_path, output_file_path)

if __name__ == "__main__":
    main()
