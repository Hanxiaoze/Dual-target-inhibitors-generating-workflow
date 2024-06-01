import os
import re
import requests
from concurrent.futures import ProcessPoolExecutor

# 设置下载路径
download_dir = 'zinc_in-man_2D_smi_files'
os.makedirs(download_dir, exist_ok=True)

def download_molecule(molecule_id):
    try:
        # 构建下载链接
        download_url = f'https://zinc.docking.org/substances/{molecule_id}'
        
        # 构建文件名
        file_name = f'{molecule_id}'
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
    matching_parts = []

    # 遍历文件夹中的所有文件
    for filename in os.listdir("./ZINC_in-man_sdf_files_2D"):
        # 使用正则表达式来查找匹配的部分
        match = re.search(r'ZINC\d{12}', filename)
        if match:
            matching_part = match.group(0)
            matching_parts.append(matching_part)

    molecule_ids = [f'{i}.smi' for i in matching_parts]
    
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(download_molecule, molecule_ids))
        
        for result in results:
            print(result)

    print("Download completed.")
