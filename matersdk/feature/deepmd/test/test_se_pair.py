import unittest

# python3 -m matersdk.feature.deepmd.test.test_se_pair
from ....io.publicLayer.structure import DStructure
from ....io.publicLayer.neigh import StructureNeighborsDescriptor
from ..se_pair import DpseTildeRPairDescriptor


class DpseTildeRPairTest(unittest.TestCase):
    def test_all_v1(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
        scaling_matrix = [5, 5, 1]
        reformat_mark = True
        coords_are_cartesian = True
        
        center_atomic_number = 42
        nbr_atomic_number = 42
        rcut = 3.2
        rcut_smooth = 3.0
        
        structure = DStructure.from_file(
                        file_format="pwmat",
                        file_path=atom_config_path)
        neighbors = StructureNeighborsDescriptor.create(
                        'v1',
                        structure=structure,
                        rcut=rcut,
                        scaling_matrix=scaling_matrix,
                        reformat_mark=reformat_mark,
                        coords_are_cartesian=coords_are_cartesian)
        
        ### Step 1. Print the attributions of DeepmdSeR
        dpse_tildeR_pair = DpseTildeRPairDescriptor.create(
                        'v1',
                        structure_neighbors=neighbors,
                        center_atomic_number=center_atomic_number,
                        nbr_atomic_number=nbr_atomic_number,
                        rcut=rcut,
                        rcut_smooth=rcut_smooth)
        print()
        print("Step 1. Print the attributions of DeepmdSeR:")
        print("\t1. deepmd_se_r.dp_feature_pair_an.shape = ", dpse_tildeR_pair.dp_feature_pair_an.shape)
        print("\t2. deepmd_se_r.dp_feature_pair_d.shape = ", dpse_tildeR_pair.dp_feature_pair_d.shape)
        print("\t3. deepmd_se_r.dp_feature_pair_rc.shape = ", dpse_tildeR_pair.dp_feature_pair_rc.shape)
        
        
        ### Step 2. Get smooth edition s_{ij}
        print()
        print("Step 2. Get segmented form of s in Deepmd:")
        print("\t1. s.shape = ", dpse_tildeR_pair._get_s(rcut=rcut, rcut_smooth=rcut_smooth).shape)
        
        
        ### Step 3.
        print()
        print("Step 3.")
        print("\t1. deepmd_se_r.dp_feature_pair_tildeR.shape = ", dpse_tildeR_pair.dp_feature_pair_tildeR.shape)


        ### Step 4.
        print()
        print("Step 4.")
        print("\t1. deepmd_se_r.get_tildeR(max_num_nbrs=24).shape = ", dpse_tildeR_pair.get_tildeR(max_num_nbrs=24).shape)
        print(dpse_tildeR_pair.get_tildeR(max_num_nbrs=24))


if __name__ == "__main__":
    unittest.main()