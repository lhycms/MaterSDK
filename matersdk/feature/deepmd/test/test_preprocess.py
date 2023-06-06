import unittest
import numpy as np
from matersdk.io.publicLayer.structure import DStructure
from matersdk.io.publicLayer.neigh import StructureNeighborsDescriptor
from matersdk.feature.deepmd.se_pair import DpseTildeRPairDescriptor

# python3 -m matersdk.feature.deepmd.test.test_preprocess
from ..preprocess import TildeRPairNormalizer


class TildRPairNormalizerTest(unittest.TestCase):
    def test_all(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/feature/movement/LiSi.config"
        structure = DStructure.from_file(
                        file_format="pwmat",
                        file_path=atom_config_path)

        scaling_matrix = [3, 3, 3]
        reformat_mark = True
        coords_are_cartesian = True

        center_atomic_number = 3    # Li
        nbr_atomic_number = 14      # Si
        rcut = 6.5
        rcut_smooth = 6.0
        
        ### Step 0. 计算某个结构的 TildeR
        struct_nbr = StructureNeighborsDescriptor.create(
                        'v1',
                        structure=structure,
                        rcut=rcut,
                        scaling_matrix=scaling_matrix,
                        reformat_mark=reformat_mark,
                        coords_are_cartesian=coords_are_cartesian)
        dpse_tildeR_pair_Li = DpseTildeRPairDescriptor.create(
                        'v1',
                        structure_neighbors=struct_nbr,
                        center_atomic_number=center_atomic_number,
                        nbr_atomic_number=3,
                        rcut=rcut,
                        rcut_smooth=rcut_smooth)
        #print(dpse_tildeR_pair.dp_feature_pair_tildeR)
        tildeRs_array_Li = dpse_tildeR_pair_Li.get_tildeR(max_num_nbrs=100)
        
        dpse_tildeR_pair_Si = DpseTildeRPairDescriptor.create(
                        'v1',
                        structure_neighbors=struct_nbr,
                        center_atomic_number=center_atomic_number,
                        nbr_atomic_number=nbr_atomic_number,
                        rcut=rcut,
                        rcut_smooth=rcut_smooth)
        #print(dpse_tildeR_pair.dp_feature_pair_tildeR)
        tildeRs_array_Si = dpse_tildeR_pair_Si.get_tildeR(max_num_nbrs=80)
        # (48, 100, 4) + (48, 80, 4) = (48, 180, 4)
        tildeRs_array = np.concatenate([tildeRs_array_Li, tildeRs_array_Si], axis=1)

    
        ### Step 1. 
        print()
        print("Step 1. ")
        normalizer = TildeRPairNormalizer(tildeRs_array=tildeRs_array)
        davg_unit, dstd_unit = normalizer.calc_stats()
        print("Step 1.1. The davg of environment matrix is :", end='\t')
        print(davg_unit)
        print("Step 1.2. The dstd of environment matrix is :", end='\t')
        print(dstd_unit)
        
        ### Step 2. 
        new_atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/feature/movement/LiSi_32.config"
        new_structure = DStructure.from_file(
                        file_format="pwmat",
                        file_path=new_atom_config_path)
        new_struct_nbr = StructureNeighborsDescriptor.create(
                        'v1',
                        structure=new_structure,
                        rcut=rcut,
                        scaling_matrix=scaling_matrix,
                        reformat_mark=reformat_mark,
                        coords_are_cartesian=coords_are_cartesian)
        new_dpse_tildeR_pair = DpseTildeRPairDescriptor.create(
                        'v1',
                        structure_neighbors=new_struct_nbr,
                        center_atomic_number=center_atomic_number,
                        nbr_atomic_number=nbr_atomic_number,
                        rcut=rcut,
                        rcut_smooth=rcut_smooth)
        #print(dpse_tildeR_pair.dp_feature_pair_tildeR)
        new_tildeRs_array = new_dpse_tildeR_pair.dp_feature_pair_tildeR
        
        print()
        print("Step 2. Using a new environment matrix, after normalize...")
        print("Step 2.1. The max value of environment is : ", end="\t")
        print(np.max(normalizer.normalize(tildeRs_array=new_tildeRs_array)))
        print("Step 2.2. The min value of environment is : ", end="\t")
        print(np.min(normalizer.normalize(tildeRs_array=new_tildeRs_array)))



if __name__ == "__main__":
    unittest.main()