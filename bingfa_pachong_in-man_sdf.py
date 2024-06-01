import os
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from concurrent.futures import ProcessPoolExecutor

# 设置下载路径
download_dir = 'ZINC_in-man_sdf_files'
os.makedirs(download_dir, exist_ok=True)

id_list = []

with open('ZINC_code_in-man.txt', 'r') as file:
    molecule_ids = file.readlines()
    for id in molecule_ids:
        id_list.append(id.strip())

def download_molecule(molecule_id):
    try:
        # 构建下载链接
        download_url = f'https://zinc.docking.org/substances/{molecule_id}.sdf'
        
        # 构建文件名
        file_name = f'{molecule_id}.sdf'
        file_path = os.path.join(download_dir, file_name)
        
        # 发起HTTP请求并下载文件
        response = requests.get(download_url)
        if response.status_code == 200:
            with open(file_path, 'wb') as f:
                f.write(response.content)
                print(f"Downloaded {file_name}")
            return f"Downloaded {file_name}"
        else:
            print(f"Failed to download {file_name}")
            return f"Failed to download {file_name}"
            
    except Exception as e:
        print(f"Error for {molecule_id}: {str(e)}")
        return f"Error for {molecule_id}: {str(e)}"

if __name__ == '__main__':
        
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(download_molecule, id_list))

        for result in results:
            print(result)

    print("Download completed.")
