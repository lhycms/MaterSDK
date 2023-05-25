import unittest

# python3 -m matersdk.io.publicLayer.test.test_neigh
from ..structure import DStructure
from ..neigh import StructureNeighborsDescriptor
from ..neigh import StructureNeighborUtils


#class StructureNeighborsDescriptor()


class StructureNeighborsV1Test(unittest.TestCase):
    def test_all(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
        #atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/feature/movement/LiSi.config"
        scaling_matrix = [5, 5, 1]
        #scaling_matrix = [3, 3, 3]
        rcut = 3.2
        #rcut = 6.5
        reformat_mark = True
        coords_are_cartesian = True   
        
        structure = DStructure.from_file(
                        file_format="pwmat",
                        file_path=atom_config_path)
        neighbors_v1 = StructureNeighborsDescriptor.create(
                        "v1",
                        structure=structure,
                        scaling_matrix=scaling_matrix,
                        reformat_mark=reformat_mark,
                        coords_are_cartesian=coords_are_cartesian,
                        rcut=rcut,
                        #max_nbrs_num=100,
                        )
        
        print()
        print("Step 1. primitive_cell 中原子的近邻原子情况:")
        print("\t1.1. The number of atoms in primitive cell:\t", len(neighbors_v1.structure.species))
        print("\t1.2. The shape of key_nbr_species:\t", neighbors_v1.key_nbr_atomic_numbers.shape)
        print("\t1.3. The shape of key_nbr_distances:\t", neighbors_v1.key_nbr_distances.shape)
        print("\t1.4. The shape of key_nbr_coords:\t", neighbors_v1.key_nbr_coords.shape)
        #print(neighbors_v1.key_nbr_atomic_numbers)



class StructureNeighborUtilsTest(unittest.TestCase):
    def all(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
        #atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/feature/movement/LiSi.config"
        scaling_matrix = [5, 5, 1]
        #scaling_matrix = [3, 3, 3]
        #rcut = 3.2
        rcut = 6.5
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
        print("Step 1. 此 atom.config 在截断半径 {0} 内，最大近邻原子数（不包括中心原子自身）为:".format(rcut), end="\t")
        print(max_num_nbrs_real)

    
if __name__ == "__main__":
    unittest.main()