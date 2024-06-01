import torch

# 创建输入张量
input = torch.tensor([[1, 2],
                      [3, 4],
                      [5, 6]])

# 创建索引张量
index = torch.tensor([[0, 1],
                      [1, 0],
                      [0, 0]])

# 在第1维度上聚合操作
output = torch.gather(input, 1, index)

print(output)
