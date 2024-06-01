import torch
import numpy as np
centers = torch.linspace(0, 1, 11)
print("centers: \n",centers)
target=torch.tensor([[0.1, 0.2, 0.3],[0.4, 0.5, 0.6],[0.7, 0.8, 0.9]])
print(target)
print("target.size: \n",target.size())
centers = centers.view(*[1 for _ in target.size()], -1)
print("centers.view(*[1 for _ in target.size()], -1) \n ",centers)
print("target.unsqueeze(-1)-centers\n",target.unsqueeze(-1)-centers)
width = 0.000001
n_bins = 11
max_size = 1
labels = torch.exp(-(1 / width) * (target.unsqueeze(-1)-centers) ** 2)
print("labels: \n", labels)
max_dist_label = torch.zeros(n_bins)
max_dist_label[-1] = 1.
print("max_dist_label: \n",max_dist_label)
labels = torch.where(target.unsqueeze(-1) <= max_size, labels, max_dist_label)
print("labels = torch.where(target.unsqueeze(-1) <= max_size, labels, max_dist_label)\n",labels)
labels = labels / torch.sum(labels, -1, keepdim=True)
print("labels = labels / torch.sum(labels, -1, keepdim=True)\n",labels)
target.unsqueeze(-1)
print("target.unsqueeze(-1) \n",target.unsqueeze(-1))
label_idcs = torch.arange(n_bins).reshape(1, 1, -1)
print("label_idcs = torch.arange(n_bins).reshape(1, 1, -1) \n", label_idcs)
t = tuple([slice(0, 3), slice(0, 4), slice(0, 5)])
print("t: \n", t)


# store all possible types in a tensor, build one-hot encoded label vectors
all_types = [1, 6, 7, 8, 9, 10]
all_types_tensor = torch.tensor(all_types)
print("all_types_tensor:\n",all_types_tensor)
type_labels = torch.eye(len(all_types))  # one-hot encoding of types
print("type_labels:\n",type_labels)
# build array that converts type to correct row index in the type labels
type_idc_converter = torch.zeros(torch.max(all_types_tensor)+1).long()
print("type_idc_converter\n",type_idc_converter)
print("all_types_tensor.long()\n",all_types_tensor.long())
type_idc_converter[all_types_tensor.long()] = torch.arange(len(all_types))
print("type_idc_converter\n",type_idc_converter)


n_tokens = 2
n_atoms = 4 + 2
neighbors = torch.ones(4, 4).long() * 8
print("neighbors: \n", neighbors)
# for i in range(n_tokens, 0, -1):
#             neighbors = \
#                 torch.cat((neighbors, torch.ones(n_atoms-i, 1).long()*n_atoms-i), 1)
#             print("neighbors: \n", neighbors)
            # neighbors = \
            #     torch.cat((neighbors, torch.arange(n_atoms-i).view(1, -1)), 0)
            # print("neighbors: \n", neighbors)

batch_size = 2
neighbors = np.arange(n_atoms-1)[None, None, :] * np.ones((batch_size, n_atoms, 1))
print("np.arange(n_atoms-1)[None, None, :]\n",np.arange(n_atoms-1)[None, None, :])
print("np.ones((batch_size, n_atoms, 1)\n",np.ones((batch_size, n_atoms, 1)))
print("neighbors = np.arange(n_atoms-1)[None, None, :] * np.ones((batch_size, n_atoms, 1))\n", neighbors)

labels = np.random.choice(100, [batch_size, n_atoms, 300])
print("np.random.choice(100, [batch_size, n_atoms, 300])\n", labels)
print("np.sum(labels, axis=-1, keepdims=True)\n", np.sum(labels, axis=-1, keepdims=True))
labels = labels / np.sum(labels, axis=-1, keepdims=True)
print("labels / np.sum(labels, axis=-1, keepdims=True)\n", labels)

dist_mask = np.ones([batch_size, n_atoms])
print("dist_mask \n",dist_mask)
dist_mask[np.random.uniform(0, 1, size=dist_mask.shape) <= 0.1] = 0
print("dist_mask[np.random.uniform(0, 1, size=dist_mask.shape) <= 0.1] = 0 \n",dist_mask)
