import numpy as np
from sklearn.neighbors import NearestNeighbors
from typing import List
from abc import ABC, abstractmethod

from .structure import DStructure


class StructureNeighborsBase(ABC):
    @abstractmethod
    def _get_key_neighs_info(self):
        pass
    
    def _get_nnbrs(self):
        '''
        Note
        ----
            1. Only implemented in `StructureNeighborsV1`
        '''
        pass
    

class StructureNeighborsV1(StructureNeighborsBase):
    def __init__(
                self,
                structure:DStructure,
                scaling_matrix:List[int]=[3,3,1],
                reformat_mark:bool=True,
                coords_are_cartesian:bool=True,
                n_neighbors:int=60,
                algorithm:str="ball_tree"):
        '''
        Description
        -----------
            1. 获取需要被分析的结构(`primitive_cell`)的近邻原子环境
            
        Parameters
        ----------
            1. structure: DStructure
                - 需要被分析的结构 （未扩包）
            2. scaling_matrix: List[int]
                - 此列表控制着 primitive_cell 的扩包倍数
                - 此处扩包是为了处理周期性边界条件
            3. reformat_mark: bool
                - 扩包后的 supercell 是否按照原子序数排序
                - 若按原子序数排序的话，则需要调用 `self.structure.get_bidx2aidx_supercell(scaling_matrix)` 获取 {`排序前site的index`: `排序后site的index`}
                - e.g. `bidx2aidx = {1: 0, 2: 1, 4: 2, 5: 3, 7: 4, 8: 5, 10: 6, 11: 7, ...}`
                - 上述 bidx2aidx 中，我们需要获取 `bidx2aidx[0], bidx2aidx[1], ..., bidx2aidx[11]` 对应的原子的neighbors
            4. coords_are_cartesian: bool
                - True: 使用笛卡尔坐标
                - False: 使用分数坐标 
            5. n_neighbors: int
                - `sklearn.neighbors.NearestNeighbors` 的参数
                - 近邻原子的数目
            6. algorihtm: str
                - `sklearn.neighbors.NearestNeighbors` 的参数
                - 采用的算法名
        '''
        self.structure = structure
        self.supercell = self.structure.make_supercell_(
                                    scaling_matrix=scaling_matrix,
                                    reformat_mark=reformat_mark)
        # key 代指 primitive_cell 中的原子
        self.key_nbr_atomic_numbers, self.key_nbr_distances, self.key_nbr_coords = \
                        self._get_key_neighs_info(
                                        scaling_matrix=scaling_matrix,
                                        n_neighbors=n_neighbors,
                                        algorithm=algorithm,
                                        coords_are_cartesian=coords_are_cartesian)

    def _get_key_neighs_info(
                        self,
                        scaling_matrix:List[int],
                        n_neighbors:int,
                        algorithm:str,
                        coords_are_cartesian:bool):
        '''
        Description
        -----------
            1. 在 supercell 中的 primitive_cell原子 的 `原子的近邻原子情况`
        
        Parameters
        ----------  
            1. scaling_matrix: List[int]
                - 扩包的倍数
                - 扩包是为了处理周期性边界条件
            2. n_neighbors: int
                - 近邻原子数目
            3. algorithm:
                - sklearn.neighbors.NearestNeighbors() 的参数
            4. coords_are_cartesian: bool
                - True: 笛卡尔坐标
                - False: 分数坐标
                
        Return
        ------
            1. nbr_species: np.ndarray
                - shape = (primitive cell中的原子数, n_neighbors)
            2. nbr_distances: np.ndarray
                - shape = (primitive cell中的原子数, n_neighbors)
            3. nbr_coords: np.ndarray
                - shape = (primitive cell中的原子数, n_neighbors, 3)
        '''
        ### Step 1. 获取 primitive_cell 原子在 supercell 中的坐标 -- `keys_coords`
        key_idxs = self.structure.get_key_idxs(scaling_matrix=scaling_matrix)
        if coords_are_cartesian:
            keys_coords = self.supercell.cart_coords[key_idxs]
        else:
            keys_coords = self.supercell.frac_coords[key_idxs]
        
        ### Step 2. 获取 `sklearn.neighbors.NearestNeighbors` object
        nnbrs = self._get_nnbrs(
                        n_neighbors=n_neighbors,
                        algorithm=algorithm,
                        coords_are_cartesian=coords_are_cartesian)
        
        ### Step 3. 得到近邻原子的元素种类(原子序数)、距中心原子的距离、坐标
        # nbr_atomic_numbers: 近邻原子的元素种类(原子序数)
        # nbr_distances: 近邻原子距中心原子的距离
        # nbr_coords: 近邻原子的坐标
        nbr_distances, nbr_idxs = nnbrs.kneighbors(keys_coords)
        ### Step 3.1. nbr_atomic_numbers --> (primitive_cell原子数, n_neighbors)
        nbr_atomic_numbers = np.array([
                            self.supercell.species[tmp_nbr_idx].Z \
                                    for tmp_nbr_idxs_tuple in nbr_idxs \
                                        for tmp_nbr_idx in tmp_nbr_idxs_tuple]).reshape(-1, n_neighbors)
        #print(nbr_species)
        
        ### Step 3.2. nbr_distances --> (primitive_cell原子数, n_neighbors)
        # 已经由 `sklearn.neighbors.NearestNeighbors.kneighbors()` 获得
        #print(nbr_distances.shape)
        
        ### Step 3.3. nbr_coords --> (primitive_cell原子数, n_neighbors, 3)
        if coords_are_cartesian:
            nbr_coords = np.array([
                    self.supercell.cart_coords[tmp_nbr_idx]
                        for tmp_nbr_idxs_tuple in nbr_idxs \
                            for tmp_nbr_idx in tmp_nbr_idxs_tuple]).flatten().reshape(-1, n_neighbors, 3)
        else:
            nbr_coords = np.array([
                    self.supercell.frac_coords[tmp_nbr_idx]
                        for tmp_nbr_idxs_tuple in nbr_idxs \
                            for tmp_nbr_idx in tmp_nbr_idxs_tuple]).flatten().reshape(-1, n_neighbors, 3)
            
        return nbr_atomic_numbers, nbr_distances, nbr_coords     


    def _get_nnbrs(
                self,
                n_neighbors:int,
                algorithm:str,
                coords_are_cartesian:bool):
        ### Step 1. 获取 sklearn.neighbors.NearestNeighbors
        nnbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm=algorithm)
        
        ### Step 2. 获取supercell 的原子坐标
        if coords_are_cartesian == True:
            coords = self.supercell.cart_coords
        else:
            coords = self.supercell.frac_coords
        
        ### Step 3. fit
        nnbrs.fit(coords)
                
        return nnbrs



class StructureNeighborsV2(StructureNeighborsBase):
    def __init__(
                self,
                structure:DStructure,
                scaling_matrix:List[int]=[3,3,1],
                reformat_mark:bool=True,
                coords_are_cartesian:bool=True,
                n_neighbors:int=60):
        self.structure = structure
        self.supercell = self.structure.make_supercell_(
                                scaling_matrix=scaling_matrix,
                                reformat_mark=reformat_mark)
        # key 代指 primitive_cell 中的原子
        self.key_nbr_atomic_numbers, self.key_nbr_distances, self.key_nbr_coords = \
                self._get_key_neighs_info(
                        scaling_matrix=scaling_matrix,
                        n_neighbors=n_neighbors,
                        coords_are_cartesian=coords_are_cartesian)
    
    
    def _get_key_neighs_info(
                self,
                scaling_matrix:List[int],
                n_neighbors:int,
                coords_are_cartesian:bool):
        '''
        Description
        -----------
            1. 
        
        Parameters
        ----------
            1. scaling_matrix: List[int]
                - 
            2. n_neighbors: int 
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
        
        ### Step 1. 获取 supercell 的各种信息，便于后面直接从其中抽取信息填写
        ###             1) nbr_atomic_numbers, nbr_coords
        ### Step 1.1. 获取 supercell 的各位点的原子序数 -- `supercell_species`
        supercell_species = np.array([tmp_site.specie.Z for tmp_site in self.supercell.sites])
        
        ### Step 1.2. 获取 supercell 的(笛卡尔)坐标 -- `supercell_coords`
        if coords_are_cartesian:
            supercell_coords = self.supercell.cart_coords
        else:
            supercell_coords = self.supercell.frac_coords

        ### Step 2. 初始化需要返回的三个 np.ndarray
        #   nbr_atomic_numbers: 近邻原子的元素种类(原子序数)
        #   nbr_distances: 近邻原子距中心原子的距离
        #   nbr_coords: 近邻原子的坐标
        # shape = (num_center, n_neighbors)
        nbr_atomic_numbers = np.zeros((len(key_idxs), n_neighbors))
        # shape = (num_center, n_neighbors)
        nbr_distances = np.zeros((len(key_idxs), n_neighbors))
        # shape = (num_center, n_neighbors, 3)
        nbr_coords = np.zeros((len(key_idxs), n_neighbors, 3))
        
        ### Step 2.1. 每个 primitive_cell 中的原子，循环一次
        for tmp_i, tmp_center_idx in enumerate(key_idxs):
            '''
            Note
            ----
                1. `tmp_i`: 从 0 开始
                2. `tmp_center_idx`: primitive_cell的原子在supercell中的index
                3. `tmp_nbr_idxs`: 距中心原子最近的 `n_neighbors` 个原子的索引 (这个索引指的是在supercell中的索引)
            '''
            ### Step 2.1.1. 计算所有原子距离该中心原子的距离
            # shape = (3,) -> (1,3)
            tmp_center_coord = supercell_coords[tmp_center_idx].reshape(1, 3)
            # shape = (num_supercell, 3)
            tmp_relative_coords = supercell_coords - tmp_center_coord
            # shape = (num_supercell,)
            distances = np.linalg.norm(tmp_relative_coords, axis=1)
            
            ### Step 2.1.2. 取出距中心原子最近的 `n_neighbors` 个原子的索引（这个索引指的是在supercell中的索引）
            ### Note: 包括了中心原子本身
            # shape = (n_neighbors,), `tmp_nbr_idxs` 是指在 supercell 中的索引
            tmp_nbr_idxs = np.argsort(distances)[:n_neighbors]
            
            ### Step 2.1.3. 根据 `Step 2.1.2` 的索引取出 atomic_number (近邻原子的原子序数，包括自身)
            # shape = (n_neighbors,)
            tmp_nbr_atomic_numbers = supercell_species[tmp_nbr_idxs]
            nbr_atomic_numbers[tmp_i, :] = tmp_nbr_atomic_numbers.reshape(1, -1)
            
            ### Step 2.1.4. 根据 `Step 2.1.2` 的索引取出 distance (近邻原子距中心原子的距离，包括自身)
            # shape = (n_neighbors,)
            tmp_nbr_distances = np.array(sorted(distances, key=lambda tmp_d: tmp_d)[:n_neighbors])
            nbr_distances[tmp_i, :] = tmp_nbr_distances.reshape(1, -1)
            
            ### Step 2.1.5. 根据 `Step 2.1.2` 的索引取出 nbr_coord (近邻原子的坐标，包括自身)
            # shape = (n_neighbors, 3)
            tmp_nbr_coord = supercell_coords[tmp_nbr_idxs, :]
            nbr_coords[tmp_i, :, :] = tmp_nbr_coord
            
        return nbr_atomic_numbers, nbr_distances, nbr_coords
