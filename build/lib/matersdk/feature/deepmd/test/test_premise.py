import unittest

# python3 -m matersdk.feature.deepmd.test.test_premise
from ....io.publicLayer.structure import DStructure
from ....io.publicLayer.neigh import StructureNeighborsDescriptor
from ..premise import DpFeaturePairPremise



class DpFeatureTest(unittest.TestCase):
    def test_all(self):
        ### Step 0.1. 
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
        #atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/feature/movement/LiSi.config"
        scaling_matrix = [5, 5, 1]
        #scaling_matrix = [3, 3, 3]
        reformat_mark = True
        n_neighbors = 200    # 需要设置得大一些
        algorithm = "ball_tree"
        coords_are_cartesian = True
        
        structure = DStructure.from_file(
                        file_format="pwmat", 
                        file_path=atom_config_path)
        neighbors = StructureNeighborsDescriptor.create(
                        "v2",
                        structure=structure,
                        scaling_matrix=scaling_matrix,
                        reformat_mark=reformat_mark,
                        n_neighbors=n_neighbors,
                        #algorithm=algorithm,
                        coords_are_cartesian=coords_are_cartesian)
        #neighbors = StructureNeighborsV2(
        #                structure=structure,
        #                scaling_matrix=scaling_matrix,
        #                reformat_mark=reformat_mark,
        #                n_neighbors=n_neighbors,
        #                coords_are_cartesian=coords_are_cartesian)
        dp_feature = DpFeaturePairPremise(structure_neighbors=neighbors)        
        
        
        ### Step 1. 抽取一对 "中心原子-近邻原子" 的 DpFeaturePairPremise
        print()
        print("Step 1. extract_feature:")
        center_atomic_number = 42
        nbr_atomic_number = 42
        rcut = 3.2
        max_num_nbrs = 10   # 需要设置的大一些
        
        dp_feature_pair_an, dp_feature_pair_d, dp_feature_pair_rc = \
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
        print("1.3. Coords -- dp_feature_pair_rc:")
        print(dp_feature_pair_rc)
        
        
        ### Step 2. The embedding of `DpFeaturePairPremise`
        #print()
        #print("Step 2. extract feature pair embedding:")
        #center_atomic_number = 42
        #nbr_atomic_number = 42
        #rcut = 3.2
        #max_num_nbrs = 10   # 需要设置的大一些
        # 
        #dp_feature_pair_embedding = \
        #            dp_feature.extract_feature_pair_embedding(
        #                            center_atomic_number=center_atomic_number,
        #                            nbr_atomic_number=nbr_atomic_number,
        #                            rcut=rcut,
        #                            max_num_nbrs=max_num_nbrs)
        #print(dp_feature_pair_embedding)
        

if __name__ == "__main__":
    unittest.main()