import os
import filecmp

def compare_folders_content(folder1, folder2):
    diff_files = []
    dcmp = filecmp.dircmp(folder1, folder2)
    
    for sub_diff in dcmp.diff_files:
        diff_files.append(sub_diff)
    
    for common_dir in dcmp.common_dirs:
        sub_diff_files = compare_folders_content(
            os.path.join(folder1, common_dir),
            os.path.join(folder2, common_dir)
        )
        diff_files.extend(sub_diff_files)
    
    return diff_files

if __name__ == '__main__':
    folder1 = '/home/zhou/ZINC_data_bank/mol2_to_pdbqt_files'
    folder2 = '/home/zhou/ZINC_data_bank/pdbqt_small_55_zinc_10_100000'
    
    diff_files = compare_folders_content(folder1, folder2)
    
    print("Files with different content:")
    for file in diff_files:
        print(file)
