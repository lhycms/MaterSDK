import numpy as np

# Define the original list
a = [3, 1, 4, 2, 5]

# Get the sorted list and its index mapping
sorted_idx = [i for i, _ in sorted(enumerate(a), key=lambda x: x)]
sorted_list = [a[i] for i in sorted_idx]
idx_map = {i: sorted_idx[i] for i in range(len(sorted_idx))}

print(sorted_list)  # Output: [1, 2, 3, 4, 5]
print(idx_map)      # Output: {0: 1, 1: 3, 2: 0, 3: 2, 4: 4}