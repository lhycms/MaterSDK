import numpy as np

from ...io.publicLayer.neigh import StructureNeighbors
from ...io.publicLayer.neigh import DpFeaturePairPremise


class DeepmdSeR(object):
    '''
    Description
    -----------
        1. In `DeepPot-SE`, features is constructed as: 
            $D^{i} = (g^{i1})^T \widetilde{R}^i (\widetilde{R}^i)^T g^{i2}$
        2. In this class, we will calculate:
            $\widetilde{R}^i$
        3. The format of $\widetilde{R}^i$ is specific, you can check it in:
            - Zhang L, Han J, Wang H, et al. End-to-end symmetry preserving inter-atomic potential energy model for finite and extended systems[J]. Advances in neural information processing systems, 2018, 31.
    '''
    def __init__(
                self,
                structure_neighbors:int,
                center_atomic_number:int,
                nbr_atomic_number:int,
                rcut:float,
                rcut_smooth:float,
                max_num_nbrs:int):
        '''
        Parameters
        ----------
            1. 
        '''
        self.dp_feature_pair_an, self.dp_feature_pair_d, self.dp_feature_pair_rc = \
                    self._get_premise(
                            structure_neighbors=structure_neighbors,
                            center_atomic_number=center_atomic_number,
                            nbr_atomic_number=nbr_atomic_number,
                            rcut=rcut,
                            max_num_nbrs=max_num_nbrs
                    )
        self.dp_feature_pair_s = self._get_s(rcut=rcut, rcut_smooth=rcut_smooth)
        
    
    def _get_premise(
                    self,
                    structure_neighbors:StructureNeighbors,
                    center_atomic_number:int,
                    nbr_atomic_number:int,
                    rcut:float,
                    max_num_nbrs:int):
        '''
        Description
        -----------
            1. 得到计算 $\widetilde{R}^i$ 所需的基本信息:
                1) 近邻原子的元素种类
                2) 近邻原子距中心原子的距离
                3) 近邻原子距中心原子的相对坐标
            2. 本类针对“中心原子-近邻原子”，即固定了中心原子、近邻原子的元素种类，因此我们将其称为`Pair`
        
        Parameters
        ----------
            1. dp_feature_pair_premise: DpFeaturePairPremise
                - 
            2. center_atomic_number: int 
                - 中心原子的原子序数
            3. nbr_atomic_number: int
                - 近邻原子的原子序数 
            4. rcut: float
                - dp 的截断半径
            5. max_num_nbrs: int 
                - 最大近邻原子数，仅针对此`Pair`
                
        Return 
        ------
            1. dp_feature_pair_an: np.ndarray, shape = (num_center, max_num_nbrs)
                - 近邻原子的元素种类
            2. dp_feature_pair_d: np.ndarray, shape = (num_center, max_num_nbrs)
                - 近邻原子距中心原子的欧式距离 (在笛卡尔坐标系下)
            3. dp_feature_pair_rc: np.ndarray, shape = (num_center, max_num_nbrs, 3)
                - 近邻原子距中心原子的相对坐标 (在笛卡尔坐标系下)
        '''
        dp_feature_pair_premise = DpFeaturePairPremise(structure_neighbors=structure_neighbors)
        dp_feature_pair_an, dp_feature_pair_d, dp_feature_pair_rc = \
                        dp_feature_pair_premise.extract_feature_pair(
                            center_atomic_number=center_atomic_number,
                            nbr_atomic_number=nbr_atomic_number,
                            rcut=rcut,
                            max_num_nbrs=max_num_nbrs
                        )
        return dp_feature_pair_an, dp_feature_pair_d, dp_feature_pair_rc
    
    
    def _get_s(self, rcut:float, rcut_smooth:float):
        '''
        Description
        -----------
            1. 由 `self.dp_feature_pair_d` 构建 `dp_feature_pair_s`
                - `s` 是一个分段函数
                - `s` 的具体形式见 Zhang L, Han J, Wang H, et al. End-to-end symmetry preserving inter-atomic potential energy model for finite and extended systems[J]. Advances in neural information processing systems, 2018, 31.
            2. 距离大于 `rcut` 的已经在 `DpFeaturePairPremise` 中被设置为 0 了
        
        Parameters
        ----------
            1. rcut: float
                - 
            2. rcut_smooth: float
                - 
            
        Note
        ----
            1. 距离大于 `rcut` 的已经在 `DpFeaturePairPremise` 中被设置为 0 了
        '''
        ### Step 1. 距离大于 `rcut` 的已经在 `DpFeaturePairPremise` 中被设置为 0 了
        
        ### Step 2. 获取 `dp_feature_pair_d_reciprocal` -- $\frac{1}{rji}$
        # (num_center, max_num_nbrs)
        dp_feature_pair_d_reciprocal = np.where(
                    self.dp_feature_pair_d==0,
                    0,
                    np.reciprocal(self.dp_feature_pair_d)
        )
        
        ### Step 3. 把`self.dp_feature_pair_d`全部转换为 rcut_smooth < r < rcut 时的形式
        # (num_center, max_num_nbrs)
        dp_feature_pair_d_scaled = np.where(
                    self.dp_feature_pair_d==0,
                    0,
                    dp_feature_pair_d_reciprocal * (1/2) * (np.cos(np.pi*(self.dp_feature_pair_d-rcut_smooth)/(rcut-rcut_smooth)) + 1)
        )
        
        ### Step 4. 根据 Step2. 和 Step3. 的结果筛选 
        dp_feature_pair_s = np.where(
                                (self.dp_feature_pair_d>rcut_smooth) & (self.dp_feature_pair_d<rcut),
                                dp_feature_pair_d_scaled,
                                dp_feature_pair_d_reciprocal)
        
        return dp_feature_pair_s