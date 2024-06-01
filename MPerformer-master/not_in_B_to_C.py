import os
import shutil

def copy_files_not_in_B_to_C(folder_A, folder_B, folder_C):    
    # 创建目标文件夹C
    if not os.path.exists(folder_C):
        os.makedirs(folder_C)

    # 获取文件夹A和B中的文件列表
    files_A = set([file for file in os.listdir(folder_A) if file.endswith('.xyz')])
    files_B = set([file[:-4] + '.xyz' for file in os.listdir(folder_B) if file.endswith('.sdf')])

    # 找出在A中但不在B中的.xyz文件，并将它们复制到文件夹C中
    for file_name in files_A - files_B:
        source_path = os.path.join(folder_A, file_name)
        destination_path = os.path.join(folder_C, file_name)
        shutil.copyfile(source_path, destination_path)
        print(f"Copied '{file_name}' from folder A to folder C.")

        

# 示例用法
folder_A = "./xyz"
folder_B = "./xxxx"
folder_C = "./xyz-xxxx"
copy_files_not_in_B_to_C(folder_A, folder_B, folder_C)
