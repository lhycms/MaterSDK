import unittest

# python3 -m pflow.feature.deepmd.test.test_se
from ....io.publicLayer.structure import DStructure
from ....io.publicLayer.neigh import StructureNeighbors
from ..se import DeepmdSeR


class DeepmdSeRTest(unittest.TestCase):
    def test_all(self):
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
        
        ### Step 1. Print the attributions of DeepmdSeR
        center_atomic_number = 42
        nbr_atomic_number = 42
        rcut = 3.4
        rcut_smooth = 3.1
        max_num_nbrs = 10
        
        deepmd_se_r = DeepmdSeR(
                        structure_neighbors=neighbors,
                        center_atomic_number=center_atomic_number,
                        nbr_atomic_number=nbr_atomic_number,
                        rcut=rcut,
                        rcut_smooth=rcut_smooth,
                        max_num_nbrs=max_num_nbrs)
        print()
        print("Step 1. Print the attributions of DeepmdSeR:")
        print(deepmd_se_r.dp_feature_pair_an)
        print(deepmd_se_r.dp_feature_pair_d)
        print(deepmd_se_r.dp_feature_pair_rc)
        
        
        ### Step 2. Get smooth edition s_{ij}
        print()
        print("Step 2. Get segmented form of s in Deepmd:")
        #print(deepmd_se_r._get_s(rcut=rcut, rcut_smooth=rcut_smooth))
        print(deepmd_se_r.dp_feature_pair_s)
        
    

if __name__ == "__main__":
    unittest.main()