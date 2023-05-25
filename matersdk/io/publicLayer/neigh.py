import numpy as np
from sklearn.neighbors import NearestNeighbors
from typing import List, Union
from abc import ABC, abstractmethod

from .structure import DStructure


class StructureNeighborUtils(object):
    @staticmethod
    def get_max_num_nbrs_real(
                    structure:DStructure,
                    scaling_matrix:List[int],
                    rcut:float,
                    coords_are_cartesian:bool=True):
        ### Step 0. 获取 primitive_cell 中原子在 supercell 中的 index
        key_idxs = structure.get_key_idxs(scaling_matrix=scaling_matrix)
        supercell = structure.make_supercell_(
                                scaling_matrix=scaling_matrix,
                                reformat_mark=True)
        
        ### Step 1. 获取supercell所有原子的(笛卡尔)坐标 -- `supercell_coords`
        if coords_are_cartesian:
            supercell_coords = supercell.cart_coords
        else:
            supercell_coords = supercell.frac_coords
        
        
        ### Step 2. 计算在截断半径内的最大原子数
        max_num_nbrs = 0
        for tmp_i, tmp_center_idx in enumerate(key_idxs):
            ### Step 2.1. 计算该中心原子与近邻原子的距离
            # shape = (3,) -> (1, 3)
            tmp_center_coord = supercell_coords[tmp_center_idx].reshape(1, 3)
            # shape = (num_supercell, 3)
            tmp_relative_coords = supercell_coords - tmp_center_coord
            # shape = (num_supercell,)
            tmp_distances = np.linalg.norm(tmp_relative_coords, axis=1)
            
            ### Step 2.2. 判断哪些近邻原子在截断半径内
            tmp_mark_rcut = np.where(tmp_distances<rcut, True, False)
            
            ### Step 2.3. 计算在截断半径内的原子数目
            tmp_num_nbrs = np.count_nonzero(tmp_mark_rcut)
            
            ### Step 2.4. 判断
            if tmp_num_nbrs > max_num_nbrs:
                max_num_nbrs = tmp_num_nbrs
        
        return max_num_nbrs - 1
                       


class StructureNeighborsDescriptor(object):
    '''
    Description
    -----------
        1. Map str to Drived class of `StructureNeighborBase`.
    
    Usage
    -----
        1. Demo 1
            ```python
            snd_v1 = StructureNeighborsDescriptor.create("v1")
            ```
    
    -----
        1. 'v1': `StructureNeighborsV1`
    '''
    registry = {}
    
    @classmethod
    def register(cls, name:str):
        def wrapper(subclass:StructureNeighborsBase):
            cls.registry[name] = subclass
        return wrapper
    
    @classmethod
    def create(cls, name:str, *args, **kwargs):
        subclass = cls.registry[name]
        if subclass is None:
            raise ValueError(f"No StructureNeighbors registered with name '{name}'")
        return subclass(*args, **kwargs)



class StructureNeighborsBase(ABC):
    @abstractmethod
    def _get_key_neighs_info(self):
        pass
    
    
    def _get_max_num_nbrs(self):
        '''
        v1
        '''
        pass
    


@StructureNeighborsDescriptor.register("v1")
class StructureNeighborsV1(StructureNeighborsBase):
    '''
    Description
    -----------
        1. Work as `StructureNeighborsV3`, but set `rcut` not `n_neighbors`.
        2. Save images(frames) in different folders, and their neighbors' size are different.
    '''
    def __init__(
                self,
                structure:DStructure,
                scaling_matrix:List[int]=[3,3,3],
                rcut:float=6.5,
                reformat_mark:bool=True,
                coords_are_cartesian:bool=True,
                max_nbrs_num:Union[bool, int]=False):
        '''
        Parameters
        ----------
            1. structure: DStructure
                - 结构
            2. scaling_matrix: List[int]
                - 扩胞系数
            3. rcut: float:
                - 截断半径
            4. reformat_mark: bool
                - 是否按照原子序数排序。
                - 这个参数一定要设置为 `True`
            5. coords_are_cartesian: bool
                - 是否使用笛卡尔坐标
            6. max_nbrs_num: 
                - False: 若不设置 `max_nbrs_num`, 则 `max_nbrs_num = max_nbrs_num_real + 1`
                - int: 超出 `max_nbrs_num_real` 的部分用 `zero padding` 填充。
        '''
        ### Step 1. 
        self.structure = structure
        self.supercell = self.structure.make_supercell_(
                                scaling_matrix=scaling_matrix,
                                reformat_mark=reformat_mark)
        self.max_nbrs_num:Union[False, int] = max_nbrs_num

        ### Step 2.
        self.key_nbr_atomic_numbers, self.key_nbr_distances, self.key_nbr_coords = \
                self._get_key_neighs_info(
                        scaling_matrix=scaling_matrix,
                        rcut=rcut,
                        coords_are_cartesian=coords_are_cartesian)
    
    
    def _get_key_neighs_info(
                self,
                scaling_matrix:List[int],
                rcut:float,
                coords_are_cartesian:bool):
        '''
        Description
        -----------
            1. 
        
        Parameters
        ----------
            1. scaling_matrix: List[int]
                - 
            2. rcut: float 截断半径
                - 
            3. coords_are_cartesian: bool
                - 
        
        Return
        ------
            1. nbr_atomic_numbers: np.ndarray, shape = (num_center, n_neighbors)
                - 
            2. nbr_distances: np.ndarray, shape = (num_center, n_neighbors)
                - 
            3. nbr_coords: np.ndarray, shape = (num_center, n_neighbors, 3)
                - 
        '''
        ### Step 0. 获取 primitive_cell 中的原子在 supercell 中的 index
        key_idxs = self.structure.get_key_idxs(scaling_matrix=scaling_matrix)
        
        ### Step 1. 获取 supercell 的各种信息 (sites的原子序数、坐标)，便于后边直接从其中抽取信息
        ###         1) supercell_atomic_numbers     2) supercell_coords
        ### Step 1.1. 获取 supercell 的各位点的原子序数 -- `supercell_atomic_numbers`
        supercell_atomic_numbers = np.array([tmp_site.specie.Z for tmp_site in self.supercell.sites])
        
        ### Step 1.2. 获取 supercell 的 (笛卡尔) 坐标 -- `supercell_coords`
        if coords_are_cartesian:
            supercell_coords = self.supercell.cart_coords
        else:
            supercell_coords = self.supercell.frac_coords
        
        ### Step 2. 初始化需要返回的三个 np.ndarray
        #   nbr_atomic_numbers: 近邻原子的元素种类 (原子序数)
        #   nbr_distances: 近邻原子距中心原子的距离
        #   nbr_coords: 近邻原子的坐标
        # shape = (num_center, n_neighbors)
        nbr_atomic_numbers = np.zeros((len(key_idxs), supercell_atomic_numbers.shape[0]))
        # shape = (num_center, n_neighbors)
        nbr_distances = np.zeros((len(key_idxs), supercell_atomic_numbers.shape[0]))
        # shape = (num_center, n_neighbors, 3)
        nbr_coords = np.zeros((len(key_idxs), supercell_atomic_numbers.shape[0], 3))
        
        
        ### Step 2.1. 每个 primitiv_cell 中的原子，循环一次
        max_num_nbrs = 0
        for tmp_i, tmp_center_idx in enumerate(key_idxs):
            '''
            Note
            ----
                1. `tmp_i`: 从 0 开始
                2. `tmp_center_idx`: primitive_cell 的原子在 supercell 中的 index
                3. `tmp_nbr_idxs`: 将 `supercell 中所有原子的索引`按照距中心原子距离的远近排序
            '''
            ### Step 2.1.1. 计算所有原子距该中心原子的距离
            # shape = (3,) -> (1, 3)
            tmp_center_coord = supercell_coords[tmp_center_idx].reshape(1, 3)
            # shape = (num_supercell, 3)
            tmp_relative_coords = supercell_coords - tmp_center_coord
            # shape = (num_supercell,)
            tmp_distances = np.linalg.norm(tmp_relative_coords, axis=1)
            
            ### Step 2.1.2. 将 `supercell 中所有原子的索引`按照距中心原子距离的远近排序（这个索引指的是在supercell中的索引）
            tmp_num_nbrs = np.count_nonzero(tmp_distances<=rcut)    # 该中心原子在截断半径内的近邻原子数 (包括自身)
            
            if tmp_num_nbrs > max_num_nbrs:
                max_num_nbrs = tmp_num_nbrs
            tmp_sorted_nbr_idxs = np.argsort(tmp_distances)[:tmp_num_nbrs]
            #print(tmp_sorted_nbr[1000:1050])
            
            ### Step 2.1.3.
            nbr_atomic_numbers[tmp_i, :tmp_num_nbrs] = supercell_atomic_numbers[tmp_sorted_nbr_idxs]
            nbr_distances[tmp_i, :tmp_num_nbrs] = tmp_distances[tmp_sorted_nbr_idxs]
            nbr_coords[tmp_i, :tmp_num_nbrs, :] = supercell_coords[tmp_sorted_nbr_idxs, :]
        
        ### Step 3. 
        if self.max_nbrs_num:
            nbr_atomic_numbers = nbr_atomic_numbers[:, :self.max_nbrs_num]
            nbr_distances = nbr_distances[:, :self.max_nbrs_num]
            nbr_coords = nbr_coords[:, :self.max_nbrs_num, :]
        else:
            nbr_atomic_numbers = nbr_atomic_numbers[:, :max_num_nbrs]
            nbr_distances = nbr_distances[:, :max_num_nbrs]
            nbr_coords = nbr_coords[:, :max_num_nbrs, :]
        
        return nbr_atomic_numbers, nbr_distances, nbr_coords