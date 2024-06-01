import os

def check_for_nan(input_folder):
    count = 1
    # 遍历输入文件夹中的所有文件
    for file_name in os.listdir(input_folder):
        if file_name.endswith('.mol2'):
            input_file_path = os.path.join(input_folder, file_name)

            # 打开文件并逐行检查是否存在'nan'
            with open(input_file_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'nan' in line:
                        print('\n', count)
                        count += 1
                        print(f"File '{input_file_path}' contains 'nan'.")
                        print(line)
                        # 如果存在'nan'，删除该文件
                        os.remove(input_file_path)
                        print(f"File '{input_file_path}' because having 'nan', has been deleted.")
                        # 退出当前文件的检查
                        break

if __name__ == "__main__":
    input_folder = "./ZINC_in-man_mol2"
    check_for_nan(input_folder)
