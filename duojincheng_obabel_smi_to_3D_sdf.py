import os
from openbabel import pybel
from multiprocessing import Pool

# 输入文件夹和输出文件夹
input_folder = './zinc_in-man_2D_smi_files'  # 替换为包含SMILES文件的文件夹路径
output_folder = './ZINC_in-man_smi_to_3D_sdf_files'  # 替换为保存SDF文件的文件夹路径

# 获取输入文件夹中的所有SMILES文件
input_files = [os.path.join(input_folder, filename) for filename in os.listdir(input_folder) if filename.endswith('.smi')]

def process_file(input_file):    
        
    output_file = os.path.join(output_folder, os.path.basename(input_file).replace('.smi', f'.sdf'))

    # 逐个读取SMILES分子并将其转换为SDF格式
    for mol in pybel.readfile("smi", input_file):
            mol.addh()  # 添加氢原子
            mol.make3D()  # 生成3D结构
            
            output_stream = pybel.Outputfile('sdf', output_file, overwrite=True)  # 创建输出流        
            output_stream.write(mol)  # 将分子写入 SDF 文件
            output_stream.close()  # 关闭输出文件
    
    print(f"{input_file} SDF file with 3D coordinates generated successfully!")
    
        

if __name__ == '__main__':
    # 创建输出文件夹
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 创建进程池并并行处理文件
    with Pool(processes=os.cpu_count()) as pool:
        pool.map(process_file, input_files)

    print("Conversion completed.")
