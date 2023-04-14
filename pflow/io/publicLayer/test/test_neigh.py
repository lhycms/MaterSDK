import unittest

# python3 -m pflow.io.publicLayer.test.test_neigh
from ..structure import DStructure
from ..neigh import StructureNeighbors


class NeighborsTest(unittest.TestCase):
    def test_all(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        scaling_matrix = [3, 3, 1]
        reformat_mark = True
        n_neighbors = 60
        algorithm = "ball_tree"
        coords_are_cartesian = True
        
        structure = DStructure.from_file(
                        file_format="pwmat", 
                        file_path=atom_config_path)
        neighbors = StructureNeighbors(
                        structure=structure,
                        scaling_matrix=scaling_matrix,
                        reformat_mark=reformat_mark,
                        n_neighbors=n_neighbors,
                        algorithm=algorithm
                        )
        
        print()
        print("Step 1. primitive_cell 中的原子在 supercell 中对应的index:", end="\t")
        print(neighbors._get_key_idxs(scaling_matrix))
        
        print()
        print("Step 2. primitive_cell 中原子的近邻原子情况:")
        key_nbr_species, key_nbr_distances, key_nbr_coords = \
                                neighbors._get_key_neighs_info(
                                        scaling_matrix=scaling_matrix,
                                        n_neighbors=n_neighbors,
                                        algorithm=algorithm,
                                        coords_are_cartesian=coords_are_cartesian)
        print("\t2.1. The number of atoms in primitive cell:\t", len(neighbors.structure.species))
        print("\t2.2. The shape of key_nbr_species:\t", key_nbr_species.shape)
        print("\t2.3. The shape of key_nbr_distances:\t", key_nbr_distances.shape)
        print("\t2.4. The shape of key_nbr_coords:\t", key_nbr_coords.shape)
        
        

    

if __name__ == "__main__":
    unittest.main()