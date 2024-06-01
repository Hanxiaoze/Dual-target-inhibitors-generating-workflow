import os
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel as ob
from openbabel import pybel

# 指定包含xyz文件的目录
input_directory = "/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/generated_1_xyz"
output_directory = "/home/zhou/sch/cG-SchNet-main/models/ZINC_in-man_no_H_small_55_2/generated/generated_1_corr_sdf"

# 获取输入目录中的所有xyz文件
input_files = [os.path.join(input_directory, file) for file in os.listdir(input_directory) if file.endswith(".xyz")]

for xyz_file in input_files:
    with open(xyz_file, 'r') as xyzf:
        lines = xyzf.readlines()
        # 解析 XYZ 文件内容
        num_atoms = int(lines[0])
        atoms_block = lines[2:2 + num_atoms]
        #print(atoms_block)
        # 创建 RDKit 分子对象
        #mol = Chem.MolFromMolBlock('\n'.join(f'{line.strip()}' for line in atoms_block))
        #mol = Chem.MolFromMol2File(xyz_file)
        
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("xyz", "smi")
        mol = ob.OBMol()
        conv.ReadFile(mol, xyz_file)
        smiles = conv.WriteString(mol)      

        
        #smiles = Chem.MolToSmiles(mol)


        try:
            # 创建包含错误结构的分子
            mol_1 = Chem.MolFromSmiles(smiles)

            # 修正键信息
            Chem.SanitizeMol(mol_1)

            # 修正原子杂化
            Chem.Kekulize(mol_1)

            # 添加氢原子
            mol_with_hydrogens = Chem.AddHs(mol_1)

            # 生成和优化3D几何构象
            AllChem.EmbedMolecule(mol_with_hydrogens)
            AllChem.UFFOptimizeMolecule(mol_with_hydrogens)

            # 输出纠正后的分子结构
            print("Successed in Corrected Structure:", Chem.MolToSmiles(mol_with_hydrogens))

            sdf_file = os.path.join(output_directory, os.path.basename(xyz_file).replace('.xyz', '.sdf'))
            

            # 创建一个 SDWriter 对象，指定输出的文件名
            sdf_writer = Chem.SDWriter(sdf_file)

            # 写入分子到 SDF 文件
            sdf_writer.write(mol_with_hydrogens)
        except:
            print(f"Failed  !!!!!!!!!!!!!!!!!!!!!!!!!!!! in {xyz_file}")
