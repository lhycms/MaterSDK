import numpy as np
from abc import ABC, abstractclassmethod

from ...io.publicLayer.neigh import StructureNeighborsV1
from .premise import DpFeaturePairPremiseDescriptor



class DpseTildeRPairBase(ABC):
    @abstractclassmethod
    def _get_premise(self):
        pass
    
    @abstractclassmethod
    def _get_s(self):
        pass
    
    @abstractclassmethod
    def _get_tildeR(self):
        pass
    


class DpseTildeRPairDescriptor(object):
    registry = {}
    
    @classmethod
    def register(cls, name:str):
        def wrapper(subclass:DpseTildeRPairBase):
            cls.registry[name] = subclass
        return wrapper

    @classmethod
    def create(cls, name:str, *args, **kwargs):
        subclass = cls.registry[name]
        if subclass is None:
            raise ValueError(f"No DpseTildeRPair registered with name '{name}'")
        return subclass(*args, **kwargs)
    
    
    
@DpseTildeRPairDescriptor.register("v1")
class DpseTildeRPairV1(DpseTildeRPairBase):
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
                structure_neighbors:StructureNeighborsV1,
                center_atomic_number:int,
                nbr_atomic_number:int,
                rcut:float,
                rcut_smooth:float):
        '''
        Parameters
        ----------
            1. structure_neighbors: StructureNeighborsV1
            2. center_atomic_number: int
            3. nbr_atomic_number: int
            4. rcut: float
            5. rcut_smooth: float
        '''
        self.dp_feature_pair_an, self.dp_feature_pair_d, self.dp_feature_pair_rc = \
                    self._get_premise(
                            structure_neighbors=structure_neighbors,
                            center_atomic_number=center_atomic_number,
                            nbr_atomic_number=nbr_atomic_number)
        self.dp_feature_pair_tildeR = self._get_tildeR(rcut=rcut, rcut_smooth=rcut_smooth)
        
    
    def _get_premise(
                    self,
                    structure_neighbors:StructureNeighborsV1,
                    center_atomic_number:int,
                    nbr_atomic_number:int):
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
            1. structure_neighbors: DpFeaturePairPremise
                - 
            2. center_atomic_number: int 
                - 中心原子的原子序数
            3. nbr_atomic_number: int
                - 近邻原子的原子序数 
                
        Return 
        ------
            1. dp_feature_pair_an: np.ndarray, shape = (num_center, max_num_nbrs_real_element)
                - 近邻原子的元素种类
            2. dp_feature_pair_d: np.ndarray, shape = (num_center, max_num_nbrs_real_element)
                - 近邻原子距中心原子的欧式距离 (在笛卡尔坐标系下)
            3. dp_feature_pair_rc: np.ndarray, shape = (num_center, max_num_nbrs_real_element, 3)
                - 近邻原子距中心原子的相对坐标 (在笛卡尔坐标系下)
        '''
        dp_feature_pair_premise = DpFeaturePairPremiseDescriptor.create(
                                    "v1",
                                    structure_neighbors=structure_neighbors)
        dp_feature_pair_an, dp_feature_pair_d, dp_feature_pair_rc = \
                dp_feature_pair_premise.extract_feature_pair(
                    center_atomic_number=center_atomic_number,
                    nbr_atomic_number=nbr_atomic_number)
                
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
                - 近邻原子距中心原子的距离超过 `rcut` 时，不予考虑
            2. rcut_smooth: float 
                - 近邻原子距中心原子的距离超过 `rcut_smooth` 时，计算对应的分段函数形式
                
        Return
        ------
            1. dp_feature_pair_s: np.ndarray, shape=(num_center, max_num_nbrs_real_element)
                - s 是根据 `近邻原子距中心原子的距离` 计算得出的，是一个分段函数形式
            
        '''        
        ### Step 1. 获取 `dp_feature_pair_d_reciprocal` -- $\frac{1}{rji}$
        # (num_center, max_num_nbrs_real_element)
        dp_feature_pair_d_reciprocal = np.reciprocal(self.dp_feature_pair_d)
        
        ### Step 2. 把`self.dp_feature_pair_d`全部转换为 rcut_smooth < r < rcut 时的形式
        # (num_center, max_num_nbrs_real_element)
        dp_feature_pair_d_scaled = dp_feature_pair_d_reciprocal * (1/2) * (np.cos(np.pi*(self.dp_feature_pair_d-rcut_smooth)/(rcut-rcut_smooth)) + 1)
        
        ### Step 3. 根据 Step2. 和 Step3. 的结果筛选 
        # (num_center, max_num_nbrs_real_element)
        dp_feature_pair_s = np.where(
                                (self.dp_feature_pair_d>rcut_smooth) & (self.dp_feature_pair_d<rcut),
                                dp_feature_pair_d_scaled,
                                dp_feature_pair_d_reciprocal)

        return dp_feature_pair_s
    
    
    def _get_tildeR(self, rcut:float, rcut_smooth:float):
        '''
        Description
        -----------
            1. Get $\widetilde{R}$ in Zhang L, Han J, Wang H, et al. End-to-end symmetry preserving inter-atomic potential energy model for finite and extended systems[J]. Advances in neural information processing systems, 2018, 31.
            
        Parameters
        ----------
            1. rcut: float
                - 近邻原子距中心原子的距离超过 `rcut` 时，不予考虑
            2. rcut_smooth: float 
                - 近邻原子距中心原子的距离超过 `rcut_smooth` 时，计算对应的分段函数形式      
            
        Return
        ------
            1. dp_feature_pair_tildeR: np.ndarray, shape=(num_center, max_num_nbrs_real_element, 4)
                - $\widetilde{R}$ in Zhang L, Han J, Wang H, et al. End-to-end symmetry preserving inter-atomic potential energy model for finite and extended systems[J]. Advances in neural information processing systems, 2018, 31.
                
        Note
        ----
            1. Haven't use zero-padding.
        '''
        ### Step 1. 调用 `self._get_s()` 得到 `dp_feature_pair_s`
        # shape = (num_center, max_num_nbrs_real_element)
        dp_feature_pair_s = self._get_s(rcut=rcut, rcut_smooth=rcut_smooth)
        # shape = (num_center, max_num_nbrs_real_element, 1)
        dp_feature_pair_s = np.repeat(
                                dp_feature_pair_s[:, :, np.newaxis],
                                1,
                                axis=2)
        
        ### Step 2. 利用 `self.dp_feature_pair_rc` 计算 $\widetilde{R}$ 的后三列
        ### Step 2.1.
        # shape = (num_center, max_num_nbrs_real_element)
        dp_feature_pair_d_reciprocal = np.reciprocal(self.dp_feature_pair_d)

        # shape = (num_center, max_num_nbrs_real_element, 1)
        dp_feature_pair_d_reciprocal = np.repeat(
                                dp_feature_pair_d_reciprocal[:, :, np.newaxis],
                                1,
                                axis=2)
        
        ### Step 2.2.
        # shape = (num_center, max_num_nbrs_real_element, 3)
        tildeR_last3 = dp_feature_pair_s * dp_feature_pair_d_reciprocal * self.dp_feature_pair_rc
        
        ### Step 3. 合并 `dp_feature_pair` 和 `tildeR_last3`
        # shape = (num_center, max_num_nbrs_real_element, 4)
        dp_feature_pair_tildeR = np.concatenate(
                                        [dp_feature_pair_s, tildeR_last3],
                                        axis=2)
        
        return dp_feature_pair_tildeR
    
    
    def get_tildeR(self, max_num_nbrs:int):
        '''
        Description
        -----------
            1. 按照 `max_num_nbrs`，对 `self.dp_feature_pair_tildeR` 进行 zero-padding
        '''
        tilde_R = np.zeros((
                    self.dp_feature_pair_tildeR.shape[0],
                    max_num_nbrs,
                    4)
        )
        tilde_R[:, :self.dp_feature_pair_tildeR.shape[1], :] = \
                    self.dp_feature_pair_tildeR
        
        return tilde_R