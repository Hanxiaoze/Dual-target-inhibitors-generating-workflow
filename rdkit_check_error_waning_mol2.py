import os
from shutil import copyfile
from rdkit import Chem
from io import StringIO
import sys

def validate_mol2(input_folder, output_folder_E):
    # 创建输出文件夹
    if not os.path.exists(output_folder_E):
        os.makedirs(output_folder_E)

    # 遍历输入文件夹中的所有文件
    for file_name in os.listdir(input_folder):
        if file_name.endswith('.mol2'):
            input_file_path = os.path.join(input_folder, file_name)

            # 尝试读取mol2文件
            try:
                # 使用 StringIO 创建一个临时缓冲区
                buffer = StringIO()

                # 重定向标准输出到缓冲区
                sys.stdout = buffer

                # 尝试使用RDKit读取mol2文件
                mol = Chem.MolFromMol2File(input_file_path, sanitize=True)

                # 恢复标准输出
                sys.stdout = sys.__stdout__

                # 获取缓冲区中的输出
                output = buffer.getvalue()
                if "Error" in output or "Can't" in output:
                    print(f"【Error】 in output or 【Can't】 in Rdkit reading mol2 file '{input_file_path}'.")
                    copyfile(input_file_path, os.path.join(output_folder_E, file_name))            

                if mol is None:
                    print(f"【mol is None】in Rdkit reading mol2 file '{input_file_path}'.")
                    # 复制到输出文件夹E
                    copyfile(input_file_path, os.path.join(output_folder_E, file_name))

            except Exception as e:
                print(f"Error: {e}")
                # 复制到输出文件夹E
                copyfile(input_file_path, os.path.join(output_folder_E, file_name))

if __name__ == "__main__":
    input_folder = "./mol2_Q"
    output_folder_E = "./mol2_Q_E"
    validate_mol2(input_folder, output_folder_E)
