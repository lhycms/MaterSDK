import numpy as np
import itertools
from sklearn.neighbors import NearestNeighbors
from typing import List

from .structure import DStructure


class StructureNeighbors(object):
    def __init__(
                self,
                structure:DStructure,
                scaling_matrix:List[int]=[3,3,1],
                reformat_mark:bool=True,
                coords_are_cartesian:bool=False,
                n_neighbors:int=50,
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
        key_idxs = self._get_key_idxs(scaling_matrix=scaling_matrix)
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
        

    
    def _get_key_idxs(
                    self,
                    scaling_matrix:List[int]
                    ):
        '''
        Description
        -----------
            1. 在 supercell 中找到 primitive_cell 所对应的原子
                1. 扩包后，primitive_cell 对应的原子处于最前几个 sites (得益于我们重写的 `DStructure.make_supercell_()`)
                2. 经历`原子序数重拍`，primitive_cell 所对应的原子不在处于最前几个sites
        
        Parameters
        ----------
            1. scaling_matrix: List[int]
                - 处理周期性边界条件时，指定的扩包倍数
        
        Return
        ------
            1. idxs_primitive_cell: List[int]
                - primitive_cell 中的原子在 supercell 中的index
                - 仍按照原子序数从小到大排序
        '''
        ### Step 1. 获取 bidx2aidx_supercell
        # bidx2aidx_supercell: supercell `重排前原子的index`: `重排后原子的index`
        bidx2aidx_supercell = self.structure.get_bidx2aidx_supercell(scaling_matrix=scaling_matrix)
        
        ### Step 2. 获取 primitive_cell 中的原子在 supercell中的index (经历`DStructure.make_supercell_()`)
        idxs_primitive_cell_ = []
        # num_atoms: primitive cell中的原子数
        num_atoms = self.structure.num_sites
        for idx in range(num_atoms):
            idxs_primitive_cell_.append(bidx2aidx_supercell[idx])
        
        ### Step 3. 按照index大小排序，本身index小意味着对应的原子序数小
        idxs_primitive_cell = sorted(idxs_primitive_cell_, key=lambda tmp_idx: tmp_idx)
        
        return idxs_primitive_cell        


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
    
    
    def extract_certain_nbr():
        '''
        Description
        -----------
            1. 仅
        '''
        pass
    
    
class AtomNeighbors(object):
    def __init__(self,
                species_array:np.ndarray,
                distances_array:np.ndarray,
                coords_array:np.ndarray):
        pass
    

class AdjacentMatrix(object):
    '''
    Description
    -----------
        1. 用于计算邻接矩阵(Adjacent Matrix)
    '''
    def __init__(self, 
                structure:DStructure,
                rcut:float,
                scaling_matrix:list
                ):
        '''
        Parameters
        ----------
            1. structure: DStructure
                - 结构
            2. rcut: float
                - 截断半径
            3. scaling_matrix: List[int]
                - 扩包倍数
                - Note: 只能是奇数
        '''
        self.structure = structure
        self.rcut = rcut
        self.scaling_matrix = scaling_matrix
        
        # 确保扩包倍数为奇数
        assert ( (self.scaling_matrix[0] > 0) and (self.scaling_matrix[0] / 2 != 0) )
        assert ( (self.scaling_matrix[1] > 0) and (self.scaling_matrix[1] / 2 != 0) )
        assert ( (self.scaling_matrix[2] > 0) and (self.scaling_matrix[2] / 2 != 0) )
    
    
    def get_adjacent_matrix(self):
        '''
        Description
        -----------
            1. 得到邻接矩阵 (adjacent matrix)
        
        
        '''
        ### Step 1. 得到原子总数和 primitive_cell 的向量基矢
        num_atoms = len(self.structure.species)
        basis_vectors_array = self.structure.lattice.matrix
        
        ### Step 2. 扩包后，周围会有很多primitive_cell，需要得到新primitive_cell 的 `分数原子坐标`
        ###     Note: 包含了本身的 primitive_cell
        neigh_primitive_cell_frac_coords = self._get_neigh_primtive_cell_frac_coords()
        
        ### Step 3. 初始化 adjacent_matrix
        adjacent_matrix = np.zeros((num_atoms, num_atoms))
        
        ### Step 4. 根据上述信息获取 adjacent matrix
        for idx_atom_1 in range(num_atoms):
            for idx_atom_2 in range(idx_atom_1+1, num_atoms):
                # shape = (3, 27)
                atom_frac_diff = neigh_primitive_cell_frac_coords[idx_atom_2] - neigh_primitive_cell_frac_coords[idx_atom_1]
                distance_ij = np.dot(basis_vectors_array.T, atom_frac_diff)
                if sum(np.linalg.norm(distance_ij, axis=0) <= self.rcut) > 0:
                    adjacent_matrix[idx_atom_1, idx_atom_2] = 1
                    adjacent_matrix[idx_atom_2, idx_atom_1] = 1
        print(adjacent_matrix)
        return adjacent_matrix
        
        
        
    
    def _get_neigh_primtive_cell_frac_coords(self):
        '''
        Description
        -----------
            1. 扩包后，周围会有很多primitive_cell，需要得到新primitive_cell 的 `分数原子坐标`
            2. Note: 包含本身的 primitive_cell
            
        Return
        ------
            1. neigh_primitive_cell_frac_coords: np.ndarray
                - shape = (12, 3, 27)
                    - `12`: 12 个原子
                    - `3` : x, y, z 三个坐标
                    - `27`: 3 * 3 * 3 (扩包倍数)
        '''
        ### Step 1. 获取扩胞前primitive_cell的 `原子总数` 和 `原子分数坐标`
        num_atoms = len(self.structure.species)
        frac_coords = np.array(self.structure.frac_coords)
        
        ### Step 2. 扩包时，各个primitive_cell相对于未扩包的primitive_cell的移动
        '''
        frac_shift_matrix 
        -----------------
            - e.g. self.scaling_matrix = [3, 3, 3]
                [[-1 -1  0]
                [-1  0  0]
                [-1  1  0]
                [ 0 -1  0]
                [ 0  0  0]
                [ 0  1  0]
                [ 1 -1  0]
                [ 1  0  0]
                [ 1  1  0]]
                ...
                ...
        '''
        shift_x = int( (self.scaling_matrix[0]-1) / 2 )
        shift_y = int( (self.scaling_matrix[1]-1) / 2 )
        shift_z = int( (self.scaling_matrix[2]-1) / 2 )
        # 转置后 shape = (3, 27)
        frac_shift_matrix_ = np.array(
                            list(itertools.product(
                                        list(range(-shift_x, shift_x+1)),
                                        list(range(-shift_y, shift_y+1)),
                                        list(range(-shift_z, shift_z+1))
                            ))
        ).T
        #print(frac_shift_matrix)
        
        ### Step 3. `neigh_primitive_cell_frac_coords`
        # shape = (12, 3, 27)
        neigh_primitive_cell_frac_coords_ = np.repeat(
                                    frac_coords[:, :, np.newaxis],
                                    self.scaling_matrix[0] * self.scaling_matrix[1] * self.scaling_matrix[2],
                                    axis=2)
        # shape = (1, 3, 27)
        frac_shift_matrix = np.repeat(
                                    frac_shift_matrix_[np.newaxis, :, :],
                                    1,
                                    axis=0)
        # shape = (12, 3, 27)
        neigh_primitive_cell_frac_coords = neigh_primitive_cell_frac_coords_ - frac_shift_matrix
        return neigh_primitive_cell_frac_coords
        

        
