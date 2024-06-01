#!/bin/bash

# conda init bash



for file in ./xyz/*.xyz; do
    # 检查当前文件夹是否存在文件
    conda run -p /home/zhou/anaconda3/envs/MPerformer python ./predict.py --filename $file --outputs_path ./sdf
    #sleep 1
done


