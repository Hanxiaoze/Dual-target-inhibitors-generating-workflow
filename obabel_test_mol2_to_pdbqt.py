#from openbabel import openbabel
from openbabel import pybel

# 输入文件和输出文件
input_file = '/home/zhou/ZINC_data_bank/ZINC_in-man___3D_add_H___mol2___108314_files_lib/3d_add_H_ZINC000085590202.mol2'  # 替换为包含 mol2 的输入文件路径
output_file = '/home/zhou/ZINC_data_bank/3d_add_H_ZINC000085590202.pdbqt'  # 替换为输出 pdbqt 文件路径

# 创建Open Babel分子对象
#obConversion = openbabel.OBConversion()
#obConversion.SetInFormat("smi")  # 设置输入格式为SMILES
#obConversion.SetOutFormat("sdf")  # 设置输出格式为SDF


# 逐个读取SMILES分子并将其转换为SDF格式
for mol in pybel.readfile("mol2", input_file):
        mol.addh()  # 添加氢原子
        mol.make3D()  # 生成3D结构
        
        output_stream = pybel.Outputfile('pdbqt', output_file, overwrite=True)  # 创建输出流        
        output_stream.write(mol)  # 将分子写入 SDF 文件
        output_stream.close()  # 关闭输出文件

    
print("pdbqt file with 3D coordinates generated.")
