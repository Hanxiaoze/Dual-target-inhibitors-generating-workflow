import os
import shutil

# 源文件夹路径和目标文件夹路径
source_folder = '/home/zhou/ZINC_data_bank/ZINC_in-man_smi_to_3D_sdf_files'  # 替换为源文件夹的路径
destination_folder = '/home/zhou/ZINC_data_bank/ZINC_in-man___3D_add_H___sdf___files_lib_for_start'  # 替换为目标文件夹的路径

# 创建目标文件夹（如果它不存在）
if not os.path.exists(destination_folder):
    os.makedirs(destination_folder)

# 遍历源文件夹中的所有文件
for filename in os.listdir(source_folder):
    source_file = os.path.join(source_folder, filename)

    # 在新文件名前添加"xxxx"
    new_filename = "3d_add_H_" + filename

    # 构建目标文件的完整路径
    destination_file = os.path.join(destination_folder, new_filename)

    # 复制文件到目标文件夹
    shutil.copy(source_file, destination_file)

print("文件复制完成。")
