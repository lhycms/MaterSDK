import unittest

# python3 -m matersdk.io.publicLayer.test.test_neigh
from ..structure import DStructure
from ..neigh import StructureNeighborsDescriptor


#class StructureNeighborsDescriptor()


class StructureNeighborsV1Test(unittest.TestCase):
    def all(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        scaling_matrix = [5, 5, 1]
        reformat_mark = True
        n_neighbors = 200
        algorithm = "ball_tree"
        coords_are_cartesian = True
        
        structure = DStructure.from_file(
                        file_format="pwmat", 
                        file_path=atom_config_path)
        #neighbors = StructureNeighborsV1(
        #                structure=structure,
        #                scaling_matrix=scaling_matrix,
        #                reformat_mark=reformat_mark,
        #                coords_are_cartesian=coords_are_cartesian,
        #                n_neighbors=n_neighbors,
        #                algorithm=algorithm
        #                )
        neighbors = StructureNeighborsDescriptor.create(
                        "v1",
                        structure=structure,
                        scaling_matrix=scaling_matrix,
                        reformat_mark=reformat_mark,
                        coords_are_cartesian=coords_are_cartesian,
                        n_neighbors=n_neighbors,
                        algorithm=algorithm)
        
        print()
        print("Step 1. primitive_cell 中原子的近邻原子情况:")
        key_nbr_species, key_nbr_distances, key_nbr_coords = \
                                neighbors._get_key_neighs_info(
                                        scaling_matrix=scaling_matrix,
                                        n_neighbors=n_neighbors,
                                        algorithm=algorithm,
                                        coords_are_cartesian=coords_are_cartesian)
        print("\t1.1. The number of atoms in primitive cell:\t", len(neighbors.structure.species))
        print("\t1.2. The shape of key_nbr_atomic_numbers:\t", key_nbr_species.shape)
        print("\t1.3. The shape of key_nbr_distances:\t", key_nbr_distances.shape)
        print("\t1.4. The shape of key_nbr_coords:\t", key_nbr_coords.shape)
        


class StructureNeighborsV2Test(unittest.TestCase):
    def test_all(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        scaling_matrix = [5, 5, 1]
        reformat_mark = True
        n_neighbors = 200
        coords_are_cartesian = True
        
        structure = DStructure.from_file(
                        file_format="pwmat", 
                        file_path=atom_config_path)
        #neighbors_v2 = StructureNeighborsV2(
        #                structure=structure,
        #                scaling_matrix=scaling_matrix,
        #                reformat_mark=reformat_mark,
        #                coords_are_cartesian=coords_are_cartesian,
        #                n_neighbors=n_neighbors)
        neighbors_v2 = StructureNeighborsDescriptor.create(
                        "v2",
                        structure=structure,
                        scaling_matrix=scaling_matrix,
                        reformat_mark=reformat_mark,
                        coords_are_cartesian=coords_are_cartesian,
                        n_neighbors=n_neighbors)
        
        
        ### Step 1.
        print()
        print("Step 1. primitive_cell 中原子的近邻原子情况:")
        
        key_nbr_species, key_nbr_distances, key_nbr_coords = \
                    neighbors_v2._get_key_neighs_info(
                                    scaling_matrix=scaling_matrix,
                                    n_neighbors=n_neighbors,
                                    coords_are_cartesian=coords_are_cartesian)
        
        print("\t1.1. The number of atoms in primitive cell:\t", len(neighbors_v2.structure.species))
        print("\t1.2. The shape of key_nbr_species:\t", key_nbr_species.shape)
        print("\t1.3. The shape of key_nbr_distances:\t", key_nbr_distances.shape)
        print("\t1.4. The shape of key_nbr_coords:\t", key_nbr_coords.shape)

    
if __name__ == "__main__":
    unittest.main()