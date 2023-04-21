import numpy as np
import itertools
from matersdk.io.publicLayer.structure import DStructure



r_cutoff = 3.2

struct = DStructure.from_file(file_format="pwmat",
                                 file_path="/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config")
struct.reformat_elements_()
print(struct)
#result = find_connected_atoms(struct=structure, tolerance=0.45)



# pylint: disable=E1136
n_atoms = len(struct.species)
fc = np.array(struct.frac_coords)
# shape = (12, 3, 27)
fc_copy = np.repeat(fc[:, :, np.newaxis], 27, axis=2)
neighbors = np.array(list(itertools.product([0, 1, -1], [0, 1, -1], [0, 1, -1]))).T
# shape = (1, 3, 27)
neighbors = np.repeat(neighbors[np.newaxis, :, :], 1, axis=0)
fc_diff = fc_copy - neighbors
species = list(map(str, struct.species))
# in case of charged species
latmat = struct.lattice.matrix
connected_matrix = np.zeros((n_atoms, n_atoms))

for i in range(n_atoms):
    for j in range(i + 1, n_atoms):
        max_bond_length = 3.2
        # (3, 27)
        frac_diff = fc_diff[j] - fc_copy[i]
        print(frac_diff)
        distance_ij = np.dot(latmat.T, frac_diff)
        # print(np.linalg.norm(distance_ij,axis=0))
        if sum(np.linalg.norm(distance_ij, axis=0) < max_bond_length) > 0:
            connected_matrix[i, j] = 1
            connected_matrix[j, i] = 1

print(connected_matrix)