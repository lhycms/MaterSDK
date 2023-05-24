from abc import ABC, abstractclassmethod
import numpy as np
from typing import Union

from ...io.publicLayer.neigh import (
                        StructureNeighborsV1,
                        StructureNeighborsV2,
                        StructureNeighborsV3)


class DpFeaturePairPremiseDescriptor(object):
    '''
    Description
    -----------
        1. Map str to Derived class of `StructureNeighborBase`.
    
    Usage
    -----
        1. Demo 1.
            ```python
            dppd = DpFeaturePairPremiseDescription("v1")
            ```
    
    ----
        1. 'v1': `DpFeaturePairPremise`
            - for `StructureNeighborsV1`, `StructureNeighborsV2`
        2. 'v2': 
            - for `StructureNeighborsV3`
    '''
    registry = {}

    @classmethod
    def register(cls, name:str):
        def wrapper(subclass:DpFeaturePairPremiseBase):
            cls.registry[name] = subclass
        return wrapper
    
    
    @classmethod
    def create(cls, name:str, *args, **kwargs):
        subclass = cls.registry[name]
        if subclass is None:
            raise ValueError(f"No DpFeaturePairPremise registered with name '{name}'")
        return subclass(*args, **kwargs)
    


class DpFeaturePairPremiseBase(ABC):
    @abstractclassmethod
    def extract_feature_pair(self):
        pass


@DpFeaturePairPremiseDescriptor.register("v1")
class DpFeaturePairPremiseV1(DpFeaturePairPremiseBase):
    def __init__(
            self,
            structure_neighbors:Union[StructureNeighborsV1, StructureNeighborsV2],
            ):
        ### `StructureNeighbors` contains informations:
        ###     1. `key_nbr_atomic_numbers` shape = (12, 60)
        ###     2. `key_nbr_distances`      shape = (12, 60)
        ###     3. `key_nbr_coords`         shape = (12, 60, 3)
        ###         12 : num_atoms
        ###         60 : n_neighbors
        ###         3  : x, y, z 
        self.structure_neighbors = structure_neighbors
        
        self.n_neighbors = structure_neighbors.key_nbr_atomic_numbers.shape[1]
        
        ### Note:
        # 函数 `self.extract_feature_pair_embedding()` 实现了PRL版DEEPMD，目前的DEEPMD是 Smooth edition
        #del self.extract_feature_pair_embedding
    
    
    def extract_feature_pair(
                self, 
                center_atomic_number:int, 
                nbr_atomic_number:int,
                rcut:float,
                max_num_nbrs:int):
        '''
        Description
        -----------
            1. 根据以下几个条件筛选所需要的信息：
                1.1. `中心原子的元素种类`
                1.2. `近邻原子的元素种类`
                1.3. `截断半径`
                1.4. `排除中心原子自身`
                1.5. `最大近邻原子数` (对于一种近邻元素来说)

        Parameters
        ----------
            1. center_atomic_number: int
                - 中心原子的原子序数
            2. nbr_atomic_number: int
                - 近邻原子的原子序数
            3. rcut: float
                - 局域描述符的截断半径
            4. max_num_nbrs: int
                - `选中的近邻原子`的最大近邻数 (决定了 `zero_padding` 的尺寸)
                - 一定要大一些，否则会报错
        
        Return 
        ------
            1. dp_feature_pair_an: shape = (num of atoms of center specie, max_num_nbrs)
                - 原子序数
                - 
            2. dp_feature_pair_d:  shape = (num of atoms of center specie, max_num_nbrs)
                - 近邻原子距中心原子的距离 (单位：埃)
                - 
            3. dp_feature_pair_rc:  shape = (num of atoms of center specie, max_num_nbrs)
                - 近邻原子的笛卡尔坐标 (不具有旋转对称性)
                - 
        
        Attribute
        ---------
            1. `num_centers`: int
                - 中心原子的数目
            2. `max_num_nbrs_num`: int
                - 近邻原子的最大数目（真实的，without zero-padding）
        
        Note
        ----
            1. `selected_knan`: 以某原子为中心，近邻原子的元素种类
                - 第一列是中心原子自身
            2. `selected_knd`: 以某原子为中心，近邻原子距其距离
                - 第一列是中心原子自身
            3. `selected_knrc`: 以某原子为中心，近邻原子距其`相对坐标`
                - 第一列是中心原子自身
        '''
        ### Step 1. 根据中心原子的原子序数，设置筛选条件 -- `mask_center`
        ### Step 1.1. 根据中心原子的原子序数的筛选条件
        mask_center = (self.structure_neighbors.key_nbr_atomic_numbers[:, 0] == center_atomic_number )        
        mask_center = np.array(mask_center).reshape(-1, 1)
        mask_center = np.repeat(mask_center, self.n_neighbors, axis=1)
        ### Step 1.2. 根据近邻原子的原子序数的筛选条件
        mask_nbr = (self.structure_neighbors.key_nbr_atomic_numbers == nbr_atomic_number)
        ### Step 1.3. 根据截断半径的筛选条件
        mask_distance = (self.structure_neighbors.key_nbr_distances <= rcut)    
        ### Step 1.4. 获得上述筛选条件的 &
        # shape = (num_centers, n_neighbors)
        mask_total = mask_center & mask_nbr & mask_distance
        ### Step 1.5. 排除中心原子自身
        mask_total[:, 0] = False      
        '''
        此步过后，
        mask_total
        ----------
            [[False False False False False False False False False False False False False False False False]
            [False False False False False False False False False False False False False False False False]
            [False False False False False False False False False False False False False False False False]
            [False False False False False False False False False False False False False False False False]
            [False False False False False False False False False False False False False False False False]
            [False False False False False False False False False False False False False False False False]
            [False False False False False False False False False False False False False False False False]
            [False False False False False False False False False False False False False False False False]
            [False False False False False False False  True  True  True  True  True True False False False]
            [False False False False False False False  True  True  True  True  True True False False False]
            [False False False False False False False  True  True  True  True  True True False False False]
            [False False False False False False False  True  True  True  True  True True False False False]]
        '''


        ### Step 2. 找到 `实际的近邻原子的最大数目` -- `max_num_nbrs_real`
        max_num_nbrs_real = np.max(np.sum(mask_total==True, axis=1), axis=0) 
        if (max_num_nbrs_real > max_num_nbrs):
            raise ValueError("Please check the atomic number you input!!! Maybe the system does not contain these kinds of element.")
        setattr(self, "max_num_nbrs_real", max_num_nbrs_real)
        setattr(self, "num_centers", np.count_nonzero(np.count_nonzero(mask_total==True, axis=1)) )
        
        ### Step 3. 根据输入的 `max_num_nbrs` 和 `num_centers` 规定输出的 np.ndarray 大小
        # shape = (4, 10)   4: 四个Mo原子 (作为中心原子)； 10: `max_num_nbr`: 根据输入设置的最大近邻原子数
        dp_feature_pair_an = np.zeros((self.num_centers, max_num_nbrs))
        # shape = (4, 10)
        dp_feature_pair_d = np.zeros((self.num_centers, max_num_nbrs))
        # shape = (4, 10, 3)
        dp_feature_pair_rc = np.zeros((self.num_centers, max_num_nbrs, 3))

        
        ### Step 4. 根据筛选条件 `mask_total` 填充 `dp_feature_pair_an`, `dp_feature_pair_d`, `dp_feature_pair_relative_c`
        ### Step 4.1. 不满足条件都化为 0 
        ### Note: 注意 `selected_knan`, `selected_knd`, `selected_knrc` 的第一列 (axis=1 为列) 均为中心原子本身
        selected_knan = np.where(
                            mask_total,
                            self.structure_neighbors.key_nbr_atomic_numbers,
                            0
        )
        selected_knd = np.where(
                            mask_total,
                            self.structure_neighbors.key_nbr_distances,
                            0
        )
        # shape = (num_center, n_neighbors, 3)
        center_coords = self.structure_neighbors.key_nbr_coords[:, 0, :]
        center_coords = np.repeat(
                            center_coords[:, np.newaxis, :],
                            1,
                            axis=1                   
        )
        knrc = self.structure_neighbors.key_nbr_coords - center_coords
        mask_total4c = np.repeat(mask_total[:, :, np.newaxis], 3, axis=2)
        selected_knrc = np.where(mask_total4c, knrc, 0)

        ### 4.2. 删除全为 0 的entry
        # shape = (12,)
        tmp_rm_zeros = ~np.all(selected_knan == 0, axis=1)
        '''
        selected_knan
        -------------
            [[ 0.  0.  0.  0.  0.  0.  0. 42. 42. 42. 42. 42. 42.  0.  0.  0.]
            [ 0.  0.  0.  0.  0.  0.  0. 42. 42. 42. 42. 42. 42.  0.  0.  0.]
            [ 0.  0.  0.  0.  0.  0.  0. 42. 42. 42. 42. 42. 42.  0.  0.  0.]
            [ 0.  0.  0.  0.  0.  0.  0. 42. 42. 42. 42. 42. 42.  0.  0.  0.]]
        '''
        selected_knan = selected_knan[tmp_rm_zeros]
        selected_knd = selected_knd[tmp_rm_zeros]
        # shape = (num_centers, n_neighbors, 3)
        selected_knrc = selected_knrc[tmp_rm_zeros]
        
        ### Step 5. 从下面这种稀疏矩阵中筛选有效信息
        for tmp_i, tmp_center_idx in enumerate(range(self.num_centers)):
            # tmp_num_nbrs: tmp_center_idx原子的 `近邻原子的原子序数`
            tmp_mask = (selected_knan[tmp_center_idx, :] != 0)
            tmp_num_nbrs = np.sum(tmp_mask, axis=0)
            dp_feature_pair_an[tmp_i, :tmp_num_nbrs] = selected_knan[tmp_center_idx][tmp_mask]
            dp_feature_pair_d[tmp_i, :tmp_num_nbrs] = selected_knd[tmp_center_idx][tmp_mask]
            dp_feature_pair_rc[tmp_i, :tmp_num_nbrs, :] = selected_knrc[tmp_center_idx][tmp_mask].reshape(-1, 3)

        return dp_feature_pair_an, dp_feature_pair_d, dp_feature_pair_rc
    
    

@DpFeaturePairPremiseDescriptor.register("v2")
class DpFeaturePairPremiseV2(DpFeaturePairPremiseBase):
    def __init__(
                self,
                structure_neighbors:StructureNeighborsV3,
                ):
        self.structure_neighbors = structure_neighbors
        self.max_num_nbrs_real = structure_neighbors.key_nbr_atomic_numbers.shape[1]
    
    
    def extract_feature_pair(
                self,
                center_atomic_number:int,
                nbr_atomic_number:int):
        '''
        Description
        -----------
            1. 根据以下几个条件筛选所需要的信息：
                1.1. `中心原子的元素种类`
                1.2. `近邻原子的元素种类`
                1.3. `排除中心原子`
                1.4. `最大近邻原子数` (对于一种近邻元素来说)，决定了 zero_padding 的尺寸
        
        Parameters
        ----------
            1. center_atomic_number: int
                - 中心原子的原子序数
            2. nbr_atomic_number: int
                - 近邻原子的原子序数
        '''
        ### Step 1. 根据中心原子的原子序数，设置筛选条件 -- `mask_center`
        ### Step 1.1. 根据中心原子的原子序数的筛选条件
        mask_center = (self.structure_neighbors.key_nbr_atomic_numbers[:, 0] == center_atomic_number)
        mask_center = np.array(mask_center).reshape(-1, 1)
        mask_center = np.repeat(mask_center, self.max_num_nbrs_real, axis=1)
        ### Step 1.2. 根据近邻原子的原子序数的筛选条件
        mask_nbr = (self.structure_neighbors.key_nbr_atomic_numbers == nbr_atomic_number)
        ### Step 1.3. 获取上述筛选条件的 &
        mask_total = mask_center & mask_nbr
        ### Step 1.4. 排除中心原子自身
        # shape = (num_centers, max_num_nbrs_real)
        '''
        [[False False False False False False False False False False False False False]
        [False False False False False False False False False False False False False]
        [False False False False False False False False False False False False False]
        [False False False False False False False False False False False False False]
        [False False False False False False False False False False False False False]
        [False False False False False False False False False False False False False]
        [False False False False False False False False False False False False False]
        [False False False False False False False False False False False False False]
        [False False False False False False False  True  True  True  True  True True]
        [False False False False False False False  True  True  True  True  True True]
        [False False False False False False False  True  True  True  True  True True]
        [False False False False False False False  True  True  True  True  True True]]
        '''
        mask_total[:, 0] = False
        #print(mask_total)
        
        ### Step 2.
        setattr( self, "num_centers", np.count_nonzero(np.count_nonzero(mask_total, axis=1)) )
        setattr( self, "max_num_nbrs_real_element", np.max(np.count_nonzero(mask_total, axis=1)))
        
        ### Step 3. 初始化返回的数组 (这些数组都不包括中心原子自身)
        #       1. `dp_feature_pair_an`:
        #       2. `dp_feature_pair_d`:
        #       3. `dp_feature_pair_rc`:
        dp_feature_pair_an = np.zeros((self.num_centers, self.max_num_nbrs_real_element))
        dp_feature_pair_d = np.zeros((self.num_centers, self.max_num_nbrs_real_element))
        dp_feature_pair_rc = np.zeros((self.num_centers, self.max_num_nbrs_real_element, 3))

        ### Step 4. 计算相对坐标
        center_coords = self.structure_neighbors.key_nbr_coords[:, 0, :]
        center_coords = np.repeat(
                            center_coords[:, np.newaxis, :],
                            1,
                            axis=1)
        # shape = (num_centers, max_num_nbrs_max+1, 3) = (12, 13, 3)
        relative_coords = self.structure_neighbors.key_nbr_coords - center_coords
        
        ### Step 5. 根据筛选条件 `mask_total` 填充 Step 3 的三个数组
        real_center_i = 0
        for tmp_i in range(mask_total.shape[0]):
            tmp_num_nbrs = np.count_nonzero(mask_total[tmp_i, :])
            if tmp_num_nbrs == 0:
                '''
                0; 0; 0; 0; 0; 0; 0; 0; 6; 6; 6; 6
                '''
                continue
            dp_feature_pair_an[real_center_i, :tmp_num_nbrs] = self.structure_neighbors.key_nbr_atomic_numbers[tmp_i][mask_total[tmp_i, :]]
            dp_feature_pair_d[real_center_i, :tmp_num_nbrs] = self.structure_neighbors.key_nbr_distances[tmp_i][mask_total[tmp_i, :]]
            dp_feature_pair_rc[real_center_i, :tmp_num_nbrs, :] = relative_coords[tmp_i][mask_total[tmp_i, :]]
            real_center_i += 1
        
        return dp_feature_pair_an, dp_feature_pair_d, dp_feature_pair_rc