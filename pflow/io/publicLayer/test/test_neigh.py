import unittest

# python3 -m pflow.io.publicLayer.test.test_neigh
from ..structure import DStructure
from ..neigh import StructureNeighbors
from ..neigh import AdjacentMatrix
from ..neigh import DpFeaturePair


class NeighborsTest(unittest.TestCase):
    def all(self):
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
        print("\t2.2. The shape of key_nbr_atomic_numbers:\t", key_nbr_species.shape)
        print("\t2.3. The shape of key_nbr_distances:\t", key_nbr_distances.shape)
        print("\t2.4. The shape of key_nbr_coords:\t", key_nbr_coords.shape)
        


class AdjacentMatrixTest(unittest.TestCase):
    def all(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        scaling_matrix = [3, 3, 3]
        structure = DStructure.from_file(
                            file_format="pwmat",
                            file_path=atom_config_path
                            )
        rcut = 3.2
        
        
        adjacent_matrix = AdjacentMatrix(
                                structure=structure,
                                rcut=rcut,
                                scaling_matrix=scaling_matrix
                                )
        
        ### Step 1. 
        print()
        print("Step 1. get_neigh_primitive_frac_coords:")
        #adjacent_matrix._get_neigh_primtive_cell_frac_coords()
        
        ### Step 2. 
        print()
        print("Step 2. The adjacent matrix (radius cutoff = {0})".format(rcut))
        adjacent_matrix.get_adjacent_matrix()



class DpFeatureTest(unittest.TestCase):
    def test_all(self):
        ### Step 0.1. 
        atom_config_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        scaling_matrix = [3, 3, 1]
        reformat_mark = True
        n_neighbors = 60    # 需要设置得大一些
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
                        algorithm=algorithm,
                        coords_are_cartesian=coords_are_cartesian)
        dp_feature = DpFeaturePair(structure_neighbors=neighbors)        
        
        
        ### Step 1. 抽取一对 "中心原子-近邻原子" 的 DpFeature
        print()
        print("Step 1. extract_feature:")
        center_atomic_number = 42
        nbr_atomic_number = 42
        rcut = 3.2
        max_num_nbrs = 10   # 需要设置的大一些
        
        dp_feature_pair_an, dp_feature_pair_d, dp_feature_pair_c = \
                    dp_feature.extract_feature_pair(
                                    center_atomic_number=center_atomic_number,
                                    nbr_atomic_number=nbr_atomic_number,
                                    rcut=rcut,
                                    max_num_nbrs=max_num_nbrs)
        print("1.1. Atomic number -- dp_feature_pair_an:")
        print(dp_feature_pair_an)
        print()
        print("1.2. Distance -- dp_feature_pair_d:")
        print(dp_feature_pair_d)
        print()
        print("1.3. Coords -- dp_feature_pair_c:")
        print(dp_feature_pair_c)
        
    
if __name__ == "__main__":
    unittest.main()