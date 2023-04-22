import numpy as np

# Create example array
arr = np.array([[1, 0, 3],
                [0, 0, 0],
                [4, 0, 6]])

# Find columns where all values are 0
all_zeros = np.all(arr == 0, axis=0)
print(all_zeros)

# Remove columns where all values are 0
arr = arr[:, ~all_zeros]

print(arr)