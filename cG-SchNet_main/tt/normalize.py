import torch

labels = torch.tensor([[[0.2, 0.3, 0.5, 0.2],
                        [0.1, 0.2, 0.7, 0.4],
                        [0.4, 0.4, 0.2, 0.6]],
                       [[0.7, 0.1, 0.1, 0.5],
                        [0.3, 0.4, 0.2, 0.3],
                        [0.5, 0.3, 0.3, 0.2]]])
print("\nlabels:")
print(labels)

print("\ntotal_sum = torch.sum(labels, -1, keepdim=True)")
total_sum = torch.sum(labels, -1, keepdim=True)
print("total_sum:")
print(total_sum)

normalized_labels = labels / total_sum
print("\nnormalized_labels = labels / total_sum")
print("normalized_labels:")
print(normalized_labels)
