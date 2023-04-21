import numpy as np
import itertools
from matersdk.io.publicLayer.structure import DStructure



r_cutoff = 3.2


structure = DStructure.from_file(file_format="pwmat",
                                 file_path="/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config")
#result = find_connected_atoms(struct=structure, tolerance=0.45)

num_atoms = len(structure.species)
fc = np.array(structure.frac_coords)    # fc.shape = (12, 3)
fc_copy = np.repeat(fc[:, :, np.newaxis], 27, axis=2)   # fc.shape = (12, 3, 27)
neighbors = np.array(list(itertools.product([0, 1, -1], [0, 1, -1], [0, 1, -1]))).T # neighbors.shape = (3, 27)
neighbors = np.repeat(neighbors[np.newaxis, :, :], 1, axis=0)   # neighbors.shape = (1, 3, 27)
fc_diff = fc_copy - neighbors   # fc_diff.shape = (12, 3, 27)


latmat = structure.lattice.matrix
adjacent_matrix = np.zeros((num_atoms, num_atoms))

for i in range(num_atoms):
    for j in range(i+1, num_atoms):
        frac_diff = fc_diff[j] - fc_diff[i]
        distance_ij = np.dot(latmat.T, frac_diff)   # 
        if sum(np.linalg.norm(distance_ij, axis=0) <= r_cutoff) > 0:
            adjacent_matrix[i, j] = 1
            adjacent_matrix[j, i] = 1
            
print(adjacent_matrix)