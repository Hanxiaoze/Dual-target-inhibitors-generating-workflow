import os
import shutil

# 源文件夹路径
source_folder = "./ZINC_in-man_sdf_files"

# 目标文件夹路径
x_folder = "./ZINC_in-man_sdf_files_2D"
y_folder = "./ZINC_in-man_sdf_files_3D"

# 搜索的字符串
string_x = "RDKit          2D"
string_y = "RDKit          3D"

# 遍历源文件夹中的文件
for root, dirs, files in os.walk(source_folder):
    for filename in files:
        file_path = os.path.join(root, filename)
        
        with open(file_path, 'r', encoding='utf-8') as file:
            file_contents = file.read()
            # 检查文件中是否包含字符串 "xxxx"
            if string_x in file_contents:
                # 构建目标文件夹路径
                target_folder = x_folder
            # 检查文件中是否包含字符串 "yyyy"
            elif string_y in file_contents:
                # 构建目标文件夹路径
                target_folder = y_folder
            else:
                # 如果没有匹配的字符串，跳过此文件
                continue

        # 确保目标文件夹存在，如果不存在则创建
        os.makedirs(target_folder, exist_ok=True)

        # 构建目标文件的完整路径
        target_path = os.path.join(target_folder, filename)

        # 移动文件到目标文件夹
        shutil.move(file_path, target_path)

        print(f"Moved '{filename}' to '{target_folder}'")

