import torch

n_batch = 2
positions = torch.tensor([
    [[1.0, 2.0, 3.0],
     [4.0, 5.0, 6.0],
#     [7.0, 8.0, 9.0],
     [10.0, 11.0, 12.0]],

    [[20.0, 21.0, 22.0],
     [23.0, 24.0, 25.0],
     [26.0, 27.0, 28.0]]
])

neighbors = torch.tensor([
    [[0, 1],
     [2, 0],
#     [0, 3],
     [1, 2]],

    [[1, 2],
     [0, 1],
     [0, 2]]
])

idx_m = torch.arange(n_batch, device=positions.device, dtype=torch.long)[:, None, None]
pos_xyz = positions[idx_m, neighbors[:, :, :], :]
dist_vec = pos_xyz - positions[:, :, None, :]

print("idx_m:")
print(idx_m)
print('\n')
print("pos_xyz:")
print(pos_xyz)
print('\n')
print("positions[:, :, None, :]:")
print(positions[:, :, None, :])
print('\n')
print("dist_vec = pos_xyz - positions[:, :, None, :] :")
print(dist_vec)
