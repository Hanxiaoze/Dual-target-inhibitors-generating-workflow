import os
from openbabel import pybel

# 输入文件夹和输出文件夹
input_folder = './sdf'
output_folder = './mol2_Q'

# 获取输入文件夹中的所有  文件
input_files = [os.path.join(input_folder, filename) for filename in os.listdir(input_folder) if filename.endswith('.sdf')]

for input_file in input_files:
    output_file = os.path.join(output_folder, os.path.basename(input_file).replace('.sdf', f'.mol2'))
    for m in pybel.readfile("sdf", input_file):   # 读取此  文件
        m.addh()  # 添加氢原子
        m.make3D()  # 生成3D结构
        #m.calccharges(charge_model="gasteiger")   # 计算 GASTEIGER 电荷
        #m.OBMol.AssignGasteigerPartialCharges()
        #m.calccharges("gasteiger")
        output_stream = pybel.Outputfile('mol2', output_file, overwrite=True)  # 创建输出流
        output_stream.write(m)  # 将分子写入 pdb 文件
        output_stream.close()  # 关闭输出文件
        print(f"{output_file} mol2 file with 【Gasteiger Q】 generated successfully!")





