import os
import subprocess
import time

input_folder = "./xyz-xxxx"
output_folder = "./sdf"
conda_env_path = "/home/zhou/anaconda3/envs/MPerformer"
python_script = "./predict.py"

# 确保输出文件夹存在
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# 遍历输入文件夹中的所有文件
for file_name in os.listdir(input_folder):
    if file_name.endswith('.xyz'):
        input_file_path = os.path.join(input_folder, file_name)
        

        # 构造调用conda环境的命令
        command = f"conda run -p {conda_env_path} python {python_script} --filename {input_file_path} --outputs_path {output_folder}"
        
        # 执行命令
        os.system(command)
        os.system('cp ./sdf/* ./xxxx/')
        # time.sleep(1)

        # 可选：休眠1秒
        # import time
        # time.sleep(1)
