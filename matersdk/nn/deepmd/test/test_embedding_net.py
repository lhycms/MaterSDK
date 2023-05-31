import unittest
import torch
import numpy as np
from matersdk.io.publicLayer.structure import DStructure
from matersdk.io.publicLayer.neigh import StructureNeighborsDescriptor
from matersdk.feature.deepmd.se_pair import DpseTildeRPairDescriptor

# python3 -m matersdk.nn.deepmd.test.test_embedding_net
from ..embedding_net import DpEmbeddingNet

import warnings
warnings.filterwarnings("ignore")


class DpEmbeddingNetTest(unittest.TestCase):
    def test_all(self):
        ### Custom parameters
        embedding_sizes = [64, 32, 16]
        
        ### Step 1. 生成 DeepPot-SE 的特征
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/feature/movement/LiSi.config"
        structure = DStructure.from_file(file_path=atom_config_path, file_format="pwmat")
        scaling_matrix = [3, 3, 3]
        center_atomic_number = 3
        nbr_atomic_number = 14
        rcut = 6.5
        rcut_smooth = 6.0
        max_num_nbrs = 36
    
        struct_nbr = StructureNeighborsDescriptor.create(
                        "v1",
                        structure=structure,
                        rcut=rcut,
                        scaling_matrix=scaling_matrix,
                        reformat_mark=True,
                        coords_are_cartesian=True)
        
        tilde_r = DpseTildeRPairDescriptor.create(
                        'v1',
                        structure_neighbors=struct_nbr,
                        center_atomic_number=center_atomic_number,
                        nbr_atomic_number=nbr_atomic_number,
                        rcut=rcut,
                        rcut_smooth=rcut_smooth).get_tildeR(max_num_nbrs=max_num_nbrs)
        print(tilde_r.shape)

        ### Step 2. 将 `tilde_R` 重复一次，假装有两个结构
        # shape = (2, 48, 200, 4)
        tilde_R = np.repeat(tilde_r[np.newaxis, :, :, :], 2, axis=0)
        tilde_R = torch.from_numpy(tilde_R)
        tilde_R = tilde_R.to(torch.float32)

        ### 2. 从 embedding_net 获取 g
        dp_embedding = DpEmbeddingNet(
                            input_dims=[1, 2, 3],
                            embedding_sizes=embedding_sizes,
                            M2=4,
                            #activation_fn
                            )
        
        dp_embedding.check()
        
        # shape = (2, 48, 16, 4). M1=16; M2=4
        D = dp_embedding(tilde_R)
        print(D.size())


if __name__ == "__main__":
    unittest.main()