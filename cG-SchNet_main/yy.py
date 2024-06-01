import torch
import numpy as np

all_types = [1, 6, 7, 8, 9, 10]
all_types_tensor = torch.tensor(all_types)
type_labels = torch.eye(len(all_types))
print("torch.max(all_types_tensor)+1\n", torch.max(all_types_tensor)+1)
type_idc_converter = torch.zeros(torch.max(all_types_tensor)+1).long()
print("type_idc_converter = torch.zeros(torch.max(all_types_tensor)+1).long()\n", type_idc_converter)
type_idc_converter[all_types_tensor.long()] = torch.arange(len(all_types))
# type_idc_converter = tensor([0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5])
print("type_idc_converter[all_types_tensor.long()] = torch.arange(len(all_types))\n", type_idc_converter)

def print_dimensions(lst):
    print(len(lst), end='')
    if isinstance(lst[0], list):
        print('x', end='')
        print_dimensions(lst[0])
    else:
        print()

n_dims = 3  # make grid in 3d space
radial_limits = [0.2 , 1.5]
n_bins = 10
max_dist = 2

grid_max = radial_limits[1]
grid_steps = int(grid_max * 2 * 20) + 1  # gives steps of length 0.05
coords = np.linspace(-grid_max, grid_max, grid_steps)
print("coords\n", coords)

grid = np.meshgrid(*[coords for _ in range(n_dims)])
print("grid\n", grid)
print("len(grid)\n", len(grid))
print("grid[0].shape()\n", grid[0].shape)
grid = np.stack(grid, axis=-1)  # stack to array (instead of list)
print("grid\n", grid)
# reshape into 2d array of positions
shape_a0 = np.prod(grid.shape[:n_dims])
print("shape_a0\n", shape_a0)
grid = np.reshape(grid, (shape_a0, -1))
print("grid\n", grid)
# cut off cells that are out of the spherical limits
grid_dists = np.sqrt(np.sum(grid**2, axis=-1))
grid_mask = np.logical_and(grid_dists >= radial_limits[0],
                            grid_dists <= radial_limits[1])
print("grid_mask\n", grid_mask)
grid = grid[grid_mask]
# assemble special grid extending only in one direction (x-axis) for the first step
# (we don't need to populate a 3d grid due to rotational invariance at first step)
start_grid = np.zeros((n_bins, n_dims))
start_grid[:, 0] = np.linspace(0, max_dist, n_bins)  # only extend along x-axis

print("grid\n",grid)
print("start_grid\n",start_grid)


coords = [1,2,3,4,5]
n_dims = 3
grid = np.meshgrid(*[coords for _ in range(n_dims)])
print("*[coords for _ in range(n_dims)]\n",*[coords for _ in range(n_dims)])
print("grid\n",grid)


finish = torch.tensor([0, 0, 1, 0], dtype=torch.bool)
print("finish[finish]\n", finish[finish])

