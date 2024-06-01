import torch
import torch.nn.functional as F

# 示例数据
batch_size = 2
n_atoms = 5
n_dims = 3
n_positions = 10

mols = torch.rand(batch_size, n_atoms, n_dims)  # 随机生成分子坐标数据
grid = torch.rand(n_positions, n_dims)  # 随机生成网格位置数据

print("batch_size = ", batch_size,"    n_atoms = ", n_atoms,"    n_dims = ",n_dims,"    n_positions = ",n_positions)

# 打印示例数据
print("示例数据 mols:")
print("mols = torch.rand(batch_size, n_atoms, n_dims)")
print(mols)
print(mols.shape)

print("\n示例数据 grid:")
print("grid = torch.rand(n_positions, n_dims)")
print(grid)
print(grid.shape)

# 如果 mols 张量的维度比 grid 张量的维度多1，添加一个额外的维度
if len(mols.size()) == len(grid.size()) + 1:
    grid = grid.unsqueeze(0)
    print("\n添加额外维度后的 grid:")
    print(grid)
    print(grid.shape)

print("\nmols[:, :, None, :]")
print(mols[:, :, None, :])
print(mols[:, :, None, :].shape)

print("\ngrid[:, None, :, :]")
print(grid[:, None, :, :])
print(grid[:, None, :, :].shape)

print("\nmols[:, :, None, :] - grid[:, None, :, :]")
print(mols[:, :, None, :] - grid[:, None, :, :])
print((mols[:, :, None, :] - grid[:, None, :, :]).shape)

print("\n(mols[:, :, None, :] - grid[:, None, :, :]).pow_(2)")
print((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2))

print("\ntorch.sum((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2), -1)")
print(torch.sum((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2), -1))
print((torch.sum((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2), -1)).shape)

print("\nF.relu(torch.sum((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2), -1), inplace=True)")
print(F.relu(torch.sum((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2), -1), inplace=True))
print((F.relu(torch.sum((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2), -1), inplace=True)).shape)

print("\nF.relu(torch.sum((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2), -1), inplace=True).sqrt_()")
print(F.relu(torch.sum((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2), -1), inplace=True).sqrt_())
print((F.relu(torch.sum((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2), -1), inplace=True).sqrt_()).shape)

# 计算每个分子的每个原子与每个位置之间的距离
result = F.relu(torch.sum((mols[:, :, None, :] - grid[:, None, :, :]).pow_(2), -1), inplace=True).sqrt_()

# 打印计算结果
print("\n计算结果:")
print(result)
