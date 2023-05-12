import unittest

# python3 -m matersdk.io.publicLayer.test.test_neigh
from ..structure import DStructure
from ..neigh import StructureNeighborsDescriptor
from ..neigh import StructureNeighborUtils


#class StructureNeighborsDescriptor()


class StructureNeighborsV1Test(unittest.TestCase):
    def all(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
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
    def all(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
        scaling_matrix = [5, 5, 1]
        reformat_mark = True
        n_neighbors = 200
        coords_are_cartesian = True
        
        structure = DStructure.from_file(
                        file_format="pwmat", 
                        file_path=atom_config_path)
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


class StructureNeighborsV3Test(unittest.TestCase):
    def all(self):
        #atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/feature/movement/LiSi.config"
        #scaling_matrix = [5, 5, 1]
        scaling_matrix = [3, 3, 3]
        #rcut = 3.2
        rcut = 6.5
        reformat_mark = True
        coords_are_cartesian = True   
        
        structure = DStructure.from_file(
                        file_format="pwmat",
                        file_path=atom_config_path)
        neighbors_v3 = StructureNeighborsDescriptor.create(
                        "v3",
                        structure=structure,
                        scaling_matrix=scaling_matrix,
                        reformat_mark=reformat_mark,
                        coords_are_cartesian=coords_are_cartesian,
                        rcut=rcut)
        
        ### Step 1.
        print()
        print("Step 1. 在截断半径 {0} 内，最大近邻原子数为:".format(rcut), end="\t")
        print(neighbors_v3.get_max_num_nbrs(
                    scaling_matrix=scaling_matrix,
                    rcut=rcut,
                    coords_are_cartesian=coords_are_cartesian)
        )
        
        print()
        print("Step 2. primitive_cell 中原子的近邻原子情况:")
        key_nbr_species, key_nbr_distances, key_nbr_coords = \
                    neighbors_v3._get_key_neighs_info(
                                    scaling_matrix=scaling_matrix,
                                    rcut=rcut,
                                    coords_are_cartesian=coords_are_cartesian)
        print("\t1.1. The number of atoms in primitive cell:\t", len(neighbors_v3.structure.species))
        print("\t1.2. The shape of key_nbr_species:\t", key_nbr_species.shape)
        print("\t1.3. The shape of key_nbr_distances:\t", key_nbr_distances.shape)
        print("\t1.4. The shape of key_nbr_coords:\t", key_nbr_coords.shape)



class StructureNeighborUtilsTest(unittest.TestCase):
    def test_all(self):
        #atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/feature/movement/LiSi.config"
        #scaling_matrix = [5, 5, 1]
        scaling_matrix = [3, 3, 3]
        #rcut = 3.2
        rcut = 6.5
        reformat_mark = True
        coords_are_cartesian = True   
        
        structure = DStructure.from_file(
                        file_format="pwmat",
                        file_path=atom_config_path)
        
        ### Step 1.
        max_num_nbrs_real = StructureNeighborUtils.get_max_num_nbrs_real(
                                    structure=structure,
                                    scaling_matrix=scaling_matrix,
                                    rcut=rcut,
                                    coords_are_cartesian=coords_are_cartesian)
        print("max_num_nbrs_real (excluding center atom self) = {0}".format(max_num_nbrs_real))


    
if __name__ == "__main__":
    unittest.main()