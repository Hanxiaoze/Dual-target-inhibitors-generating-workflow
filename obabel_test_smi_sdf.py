#from openbabel import openbabel
from openbabel import pybel

# 输入文件和输出文件
input_file = '/home/zhou/ZINC_data_bank/zinc_in-man_222D_smi_files/ZINC000000899435.smi'  # 替换为包含SMILES的输入文件路径
output_file = '/home/zhou/ZINC_data_bank/output.sdf'  # 替换为输出SDF文件路径

# 创建Open Babel分子对象
#obConversion = openbabel.OBConversion()
#obConversion.SetInFormat("smi")  # 设置输入格式为SMILES
#obConversion.SetOutFormat("sdf")  # 设置输出格式为SDF


# 逐个读取SMILES分子并将其转换为SDF格式
for mol in pybel.readfile("smi", input_file):
        mol.addh()  # 添加氢原子
        mol.make3D()  # 生成3D结构
        
        output_stream = pybel.Outputfile('sdf', output_file, overwrite=True)  # 创建输出流        
        output_stream.write(mol)  # 将分子写入 SDF 文件
        output_stream.close()  # 关闭输出文件

    
print("SDF file with 3D coordinates generated.")
