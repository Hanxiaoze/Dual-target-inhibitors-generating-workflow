import os
import shutil

def copy_files_without_string(source_folder, destination_folder, search_string):
    os.makedirs(destination_folder, exist_ok=True)
    
    for file in [file for file in os.listdir(source_folder) if os.path.isfile(os.path.join(source_folder, file))]:
        
        source_file_path = os.path.join(source_folder, file)
        destination_file_path = os.path.join(destination_folder, file)
            
        with open(source_file_path, 'r') as source_file:
            file_contents = source_file.read()
            if search_string not in file_contents:
                shutil.copy(source_file_path, destination_file_path)
                print(f"File '{file}' copied to destination folder.")

if __name__ == '__main__':
    source_folder = './ZINC_in-man_pdbqt'
    destination_folder = './ZINC_in-man_pdbqt_55'
    search_string = 'ATOM     56' # 要搜索的字符串

    copy_files_without_string(source_folder, destination_folder, search_string)

    print("File copying completed.")
