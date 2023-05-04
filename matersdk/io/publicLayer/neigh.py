import numpy as np
from sklearn.neighbors import NearestNeighbors
from typing import List
from abc import ABC, abstractmethod

from .structure import DStructure


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
            - sklearn.NearestNeighbors
        2. 'v2': `StructureNeighborsV2`
            - Custom KNN
        3. 'v3': `StructureNeighborsV3`
    '''
    registry = {}
    
    @classmethod
    def register(cls, name:str):
        def wrapper(subclass):
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
    
    def _get_nnbrs(self):
        '''
        Note
        ----
            1. Only implemented in `StructureNeighborsV1`
        '''
        pass
    
    #@abstractmethod
    def _judge_rationality(self):
        '''
        v1, v2
        '''
        pass
    
    def _get_max_num_nbrs(self):
        '''
        v3
        '''
        pass
    
    
    
@StructureNeighborsDescriptor.register("v1")
class StructureNeighborsV1(StructureNeighborsBase):
    def __init__(
                self,
                structure:DStructure,
                scaling_matrix:List[int]=[3,3,1],
                reformat_mark:bool=True,
                coords_are_cartesian:bool=True,
                n_neighbors:int=200,
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
            
        Note
        ----
            1. `n_neighbors` is the number of all kinds of atoms
                - So, it must be larger than `max_num_nbrs` in `DeepPot-SE`
            2. You must tune `scaling_matrix` larger when tuning `n_neighbors` larger
        '''
        ### Step 1.
        self.structure = structure
        self.supercell = self.structure.make_supercell_(
                                    scaling_matrix=scaling_matrix,
                                    reformat_mark=reformat_mark)
        
        ### Step 2.
        if self._judge_rationality(n_neighbors=n_neighbors):
            raise ValueError("Now the number of atoms in supercell is less than n_neighbors! You must tune scaling_matrix larger.")
        
        ### Step 3.
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
            1. nbr_atomic_numbers: np.ndarray
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
    
    def _judge_rationality(self, n_neighbors:int):
        '''
        Description
        ----------- 
            1. 判断经过 `scaling_matrix`
        '''
        return (self.supercell.num_sites < n_neighbors)



@StructureNeighborsDescriptor.register("v2")
class StructureNeighborsV2(StructureNeighborsBase):
    '''
    Description
    -----------
        1. Work as `StructureNeighborsV1`, but rewrite KNN myself.
    '''
    def __init__(
                self,
                structure:DStructure,
                scaling_matrix:List[int]=[3,3,1],
                reformat_mark:bool=True,
                coords_are_cartesian:bool=True,
                n_neighbors:int=200):
        '''
        Note
        ----
            1. `n_neighbors` is the number of all kinds of atoms
                - So, it must be larger than `max_num_nbrs` in `DeepPot-SE`
        '''
        ### Step 1. 
        self.structure = structure
        self.supercell = self.structure.make_supercell_(
                                scaling_matrix=scaling_matrix,
                                reformat_mark=reformat_mark)

        ### Step 2.
        if self._judge_rationality(n_neighbors=n_neighbors):
            raise ValueError("Now the number of atoms in supercell is less than n_neighbors! You must tune scaling_matrix larger.")

        ### Step 3.         
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
        
        ### Step 1. 获取 supercell 的各种信息（sites的原子序数、坐标），便于后面直接从其中抽取信息填写
        ###             1) supercell_atomic_numbers, supercell_coords
        ### Step 1.1. 获取 supercell 的各位点的原子序数 -- `supercell_atomic_numbers`
        supercell_atomic_numbers = np.array([tmp_site.specie.Z for tmp_site in self.supercell.sites])
        
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
            tmp_nbr_atomic_numbers = supercell_atomic_numbers[tmp_nbr_idxs]
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


    def _judge_rationality(self, n_neighbors:int):
        '''
        Description
        ----------- 
            1. 判断经过 `scaling_matrix`
        '''
        return (self.supercell.num_sites < n_neighbors)
    


@StructureNeighborsDescriptor.register("v3")
class StructureNeighborsV3(StructureNeighborsBase):
    '''
    Description
    -----------
        1. Work as `StructureNeighborsV3`, but set `rcut` not `n_neighbors`.
    '''
    def __init__(
                self,
                structure:DStructure,
                scaling_matrix:List[int]=[3,3,1],
                reformat_mark:bool=True,
                coords_are_cartesian:bool=True,
                rcut:float=6.5):
        '''
        Parameters
        ----------
            1. structure: DStructure
                - 结构
            2. scaling_matrix: List[int]
                - 扩胞系数
            3. reformat_mark: bool
                - 是否按照原子序数排序。
                - 这个参数一定要设置为 `True`
            4. coords_are_cartesian: bool
                - 是否使用笛卡尔坐标
            5. rcut: float:
                - 截断半径
        '''
        ### Step 1. 
        self.structure = structure
        self.supercell = self.structure.make_supercell_(
                                scaling_matrix=scaling_matrix,
                                reformat_mark=reformat_mark)

        ### Step 2.
        #self.max_num_nbrs = self.get_max_num_nbrs(
        #                            scaling_matrix=scaling_matrix,
        #                            rcut=rcut,
        #                            coords_are_cartesian=coords_are_cartesian)

        ### Step 3. 
        # key 代指 primitive 中的原子
    
    
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
            
            ### Step 2.1.3.
            nbr_atomic_numbers[tmp_i, :tmp_num_nbrs] = supercell_atomic_numbers[tmp_sorted_nbr_idxs]
            nbr_distances[tmp_i, :tmp_num_nbrs] = tmp_distances[tmp_sorted_nbr_idxs]
            nbr_coords[tmp_i, :tmp_num_nbrs, :] = supercell_coords[tmp_sorted_nbr_idxs, :]
        
        ### Step 3. 
        nbr_atomic_numbers = nbr_atomic_numbers[:, :max_num_nbrs]
        nbr_distances = nbr_distances[:, :max_num_nbrs]
        nbr_coords = nbr_coords[:, :max_num_nbrs, :]
        
        return nbr_atomic_numbers, nbr_distances, nbr_coords
            
    
    def get_max_num_nbrs(
                        self,
                        scaling_matrix:List[int],
                        rcut:float,
                        coords_are_cartesian:bool=True):
        '''
        Description
        -----------
            1. 得到在一定的截断半径内，各中心原子的最大近邻原子数
        
        Parameters
        ----------
            1. scaling_matrix: List[int]
                - 扩胞倍数
            2. rcut: float
                - 截断半径
            3. coords_are_cartesian: bool
                - 是否使用笛卡尔坐标
        
        Note
        ----
            1. 包含原子自身，因此近邻原子数其实应该在此基础上 `-1`
        '''
        ### Step 0. 获取 primitive_cell 中原子在 supercell 中的 index
        key_idxs = self.structure.get_key_idxs(scaling_matrix=scaling_matrix)
        
        ### Step 1. 获取 supercell 的各种信息(sites的原子序数、坐标)，便于后边直接提取
        #               1) supercell_atomic_numbers, supercell_coords
        ### Step 1.1. 获取 supercell 的各位点的原子序数 -- `supercell_atomic_numbers`
        #supercell_atomic_numbers = np.array([tmp_site.specie.Z for tmp_site in self.supercell.sites])
        
        ### Step 1.2. 获取 supercell 的(笛卡尔)坐标 -- `supercell_coords`
        if coords_are_cartesian:
            supercell_coords = self.supercell.cart_coords
        else:
            supercell_coords = self.supercell.frac_coords
        
        ### Step 2. 计算在截断半径内的最大原子数
        max_num_nbrs = 0
        for tmp_i, tmp_center_idx in enumerate(key_idxs):
            ### Step 2.1. 计算该中心原子与近邻原子的距离
            # shape = (3,) -> (1,3)
            tmp_center_coord = supercell_coords[tmp_center_idx].reshape(1, 3)
            # shape = (num_supercell, 3)
            tmp_relative_coords = supercell_coords - tmp_center_coord
            # shape = (num_supercell)
            distances = np.linalg.norm(tmp_relative_coords, axis=1)

            ### Step 2.2. 判断哪些近邻原子在截断半径内
            tmp_mark_rcut = np.where(distances<rcut, True, False)
            #print(tmp_mark_rcut)
            
            ### Step 2.3. 计算在截断半径内的近邻原子数目
            tmp_num_nbrs = np.count_nonzero(tmp_mark_rcut)
            
            ### Step 2.4. 判断
            if tmp_num_nbrs > max_num_nbrs:
                max_num_nbrs = tmp_num_nbrs
        
        return max_num_nbrs