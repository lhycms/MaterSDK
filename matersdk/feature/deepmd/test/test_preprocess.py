import unittest
from matersdk.io.publicLayer.structure import DStructure
from matersdk.io.publicLayer.neigh import StructureNeighborsDescriptor
from matersdk.feature.deepmd.se_pair import DpseTildeRPairDescriptor

# python3 -m matersdk.feature.deepmd.test.test_preprocess
from ..preprocess import TildeRPairNormalizer


class TildRNormalizerTest(unittest.TestCase):
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
        dpse_tildeR_pair = DpseTildeRPairDescriptor.create(
                        'v1',
                        structure_neighbors=struct_nbr,
                        center_atomic_number=center_atomic_number,
                        nbr_atomic_number=nbr_atomic_number,
                        rcut=rcut,
                        rcut_smooth=rcut_smooth)
        #print(dpse_tildeR_pair.dp_feature_pair_tildeR.shape)
        tildeRs_array = dpse_tildeR_pair.dp_feature_pair_tildeR
    
    
        ### Step 1. 
        normalizer = TildeRPairNormalizer(tildeRs_array=tildeRs_array)
        davg_unit, dstd_unit = normalizer.calc_stats()
        print(davg_unit)
        print(dstd_unit)



if __name__ == "__main__":
    unittest.main()