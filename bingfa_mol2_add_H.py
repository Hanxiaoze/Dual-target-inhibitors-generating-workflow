import os
import multiprocessing
from openbabel import openbabel

# 指定包含MOL2文件的目录
input_directory = "/home/zhou/ZINC_data_bank/ZINC_in-man_mol2_files"
output_directory = "/home/zhou/ZINC_data_bank/ZINC_in-man_mol2_add_H_files"

# 获取输入目录中的所有MOL2文件
input_files = [os.path.join(input_directory, file) for file in os.listdir(input_directory) if file.endswith(".mol2")]

# 定义一个函数来处理单个MOL2文件
def process_mol2_file(input_file):
    # 构建输出文件的路径
    output_file = os.path.join(output_directory, os.path.basename(input_file))
    
    # 创建Open Babel的分子对象
    mol = openbabel.OBMol()
    
    # 读取MOL2文件
    conv = openbabel.OBConversion()
    conv.SetInFormat("MOL2")
    conv.ReadFile(mol, input_file)
    
    # 添加氢原子
    mol.AddHydrogens()
    
    # 保存处理后的分子到输出文件
    conv.SetOutFormat("MOL2")
    conv.WriteFile(mol, output_file)
    
    print(f"Processed {input_file}")

if __name__ == "__main__":
    # 设置并行进程的数量
    num_processes = multiprocessing.cpu_count()
    
    # 创建进程池
    with multiprocessing.Pool(num_processes) as pool:
        pool.map(process_mol2_file, input_files)