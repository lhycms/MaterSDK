import numpy as np
from sklearn.neighbors import NearestNeighbors
from typing import List

from .structure import DStructure


class StructureNeighbors(object):
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



class StructureNeighborsV2(object):
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


    
class AtomNeighbors(object):
    def __init__(self,
                species_array:np.ndarray,
                distances_array:np.ndarray,
                coords_array:np.ndarray):
        pass
        


class DpFeaturePairPremise(object):
    def __init__(
                self,
                structure_neighbors:StructureNeighbors,
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
                1.4. `最大近邻原子数` (对于一种近邻元素来说)

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
            
        Step
        ----
            1. 按照 `中心原子的原子序数` 筛选
            2. 按照 `截断半径` 筛选
            3. 按照 `近邻原子的原子序数` 筛选
            4. 
        '''
        ### Step 1. 根据中心原子的原子序数，设置筛选条件 -- `mask_center`
        ###     e.g. 12原子的 MoS2: `mask_center = [False False False False False False False False  True  True  True  True]`
        mask_center = (self.structure_neighbors.key_nbr_atomic_numbers[:, 0] == center_atomic_number )
        '''
        mask_center.reshape(-1, 1)
        --------------------------
                [
                    [False]
                    [False]
                    [False]
                    [False]
                    [False]
                    [False]
                    [False]
                    [False]
                    [ True]
                    [ True]
                    [ True]
                    [ True]
                ]
        '''
        # shape = (12, 1)
        mask_center = np.array(mask_center).reshape(-1, 1)
        # shape = (12, 60)
        mask4key_nbr_atomic_number = np.repeat(
                                            mask_center,
                                            self.n_neighbors,
                                            axis=1)
    
        ### Step 2. 筛选出 `center_specie` 的：
        ###     1. selected_k_neighbors_atomic_numbers:
        ###     2. selected_k_neighbors_distances:
        ###     3. selected_k_neighbors_coords:
        ### Step 2.1. 获取所有 center_specie 的近邻元素种类
        # shape = (4, 60)
        selected_knan = self.structure_neighbors.key_nbr_atomic_numbers[mask4key_nbr_atomic_number].reshape(-1, self.n_neighbors)
        ### Step 2.2. 获取所有 center_specie 的近邻距离
        # shape = (4, 60)
        selected_knd = self.structure_neighbors.key_nbr_distances[mask4key_nbr_atomic_number].reshape(-1, self.n_neighbors)
        ### Step 2.3. 获取所有 center_specie 的近邻原子坐标
        # shape = (12, 60, 3)
        mask4knc = np.repeat(
                        mask4key_nbr_atomic_number[:, :, np.newaxis],
                        3,
                        axis=2)
        ### Note: 计算各近邻原子距中心原子的相对坐标 (在笛卡尔坐标系下)
        # shape = (12, 3)
        center_coords = self.structure_neighbors.key_nbr_coords[:, 0, :]
        # shape = (12, 1, 3)
        center_coords = np.repeat(center_coords[:, np.newaxis, :], 1, axis=1)
        # shape = (12, 60, 3)
        selected_knrc = self.structure_neighbors.key_nbr_coords[:, :, :] - center_coords
        # shape = (4, 60, 3)
        selected_knrc = selected_knrc[mask4knc].reshape(-1, self.n_neighbors, 3)
        
        ### Step 3. 按照 `距中心原子的距离` 筛选 (`rcut`: 局域描述符的截断半径)
        # shape = (4, 60)
        mask4key_nbr_distances = (selected_knd <= rcut)
        # shape = (4, 60)
        selected_knan = np.where(mask4key_nbr_distances, selected_knan, 0)
        # shape = (4, 60)
        selected_knd = np.where(mask4key_nbr_distances, selected_knd, 0)
        # shape = (4, 60, 3)
        mask4knc = np.repeat(
                    mask4key_nbr_distances[:, :, np.newaxis],
                    3, 
                    axis=2)
        # shape = (4, 60, 3)
        selected_knrc = np.where(mask4knc, selected_knrc, 0).reshape(-1, self.n_neighbors, 3)


        ### Step 4. 按照近邻原子的元素种类 (`nbr_atomic_number`) 筛选 ("Mo"-"S"), 此步则筛选近邻的S        
        ###     1. `selected_knan`: atomic_number
        ###     2. `selected_knd` : distances
        ###     3. `selected_knrc` :  relative coords (相对于center_atom的坐标，在笛卡尔坐标系下)
        # shape = (4, 60)
        mask4nbr_atomic_number = (selected_knan == nbr_atomic_number)
        # shape = (4, 60)
        '''
        selected_knan
        -------------
            [[42  0  0  0  0  0  0 42 42 42 42 42 42  0  0  0  0  0  0  0  0  0  0  0
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
            0  0  0  0  0  0  0  0  0  0  0  0]
            [42  0  0  0  0  0  0 42 42 42 42 42 42  0  0  0  0  0  0  0  0  0  0  0
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
            0  0  0  0  0  0  0  0  0  0  0  0]
            [42  0  0  0  0  0  0 42 42 42 42 42 42  0  0  0  0  0  0  0  0  0  0  0
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
            0  0  0  0  0  0  0  0  0  0  0  0]
            [42  0  0  0  0  0  0 42 42 42 42 42 42  0  0  0  0  0  0  0  0  0  0  0
            0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
            0  0  0  0  0  0  0  0  0  0  0  0]]
            
        Note
        ----
            1. selected_knan 的第一列是中心原子自身
        '''
        # shape = (4, 60)  --  (self.num_centers, self.neighbors)
        selected_knan = np.where(mask4nbr_atomic_number, selected_knan, 0)
        # shape = (4, 60)
        selected_knd = np.where(mask4nbr_atomic_number, selected_knd, 0)
        # shape = (4, 60, 3)
        mask4knc = np.repeat(
                        mask4nbr_atomic_number[:, :, np.newaxis],
                        3,
                        axis=2)
        selected_knrc = np.where(mask4knc, selected_knrc, 0)
        
        
        ### Step 5. 找到 `实际的近邻原子的最大数目` -- `max_num_nbrs_real`
        mask4knan = (selected_knan != 0)
        max_num_nbrs_real = np.max(np.sum(mask4knan, axis=1), axis=0) - 1
        setattr(self, "max_num_nbrs_real", max_num_nbrs_real)
        setattr(self, "num_centers", selected_knan.shape[0])
        
        if self.max_num_nbrs_real > max_num_nbrs:
            raise ValueError("max_num_nbrs_real > max_num_nbr. Please modify max_num_nbr larger!")

        ### Step 6. 根据输入的 `max_num_nbrs` 和 `num_centers` 规定输出的 np.ndarray 大小
        # shape = (4, 10)   4: 四个Mo原子 (作为中心原子)； 10: `max_num_nbr`: 根据输入设置的最大近邻原子数
        dp_feature_pair_an = np.zeros((self.num_centers, max_num_nbrs))
        # shape = (4, 10)
        dp_feature_pair_d = np.zeros((self.num_centers, max_num_nbrs))
        # shape = (4, 10, 3)
        dp_feature_pair_rc = np.zeros((self.num_centers, max_num_nbrs, 3))

        ### Step 7. 根据之前信息填充 `dp_feature_pair_an`, `dp_feature_pair_d`, `dp_feature_pair_relative_c`
        ### Note: 注意 `selected_knan`, `selected_knd`, `selected_knrc` 的第一列 (axis=1 为列) 均为中心原子本身
        # shape = (4, 60)
        mask4not_zero = (selected_knan != 0)
        # shape (4, 60)
        mask4not_zero[:, 0] = False
        
        for idx_center in range(mask4not_zero.shape[0]):    # 每个 center_atom 循环一轮
            # shape = (60,)
            tmp_mask = mask4not_zero[idx_center, :]
            # tmp_num_nbrs: int -- 这个中心原子(center_atomic_number)的近邻原子(nbr_atomic_number)
            tmp_num_nbrs = np.sum(tmp_mask, axis=0) 
            
            ### Step 7.1. Set `dp_feature_pair_atomic_number` -- tmp_selected_knan
            tmp_selected_knan = selected_knan[idx_center, :][tmp_mask]
            dp_feature_pair_an[idx_center, :tmp_num_nbrs] = tmp_selected_knan

            ### Step 7.2. Set `dp_feature_pair_distances` -- tmp_selected_knd
            tmp_selected_knd = selected_knd[idx_center, :][tmp_mask]
            dp_feature_pair_d[idx_center, :tmp_num_nbrs] = tmp_selected_knd
        
            ### Step 7.3. Set `dp_feature_pair_relative_coords` -- tmp_selected_knrc
            tmp_mask_c = np.repeat(tmp_mask[:, np.newaxis], 3, axis=1)
            tmp_selected_knrc = selected_knrc[idx_center, :][tmp_mask_c].reshape(-1, 3)
            dp_feature_pair_rc[idx_center, :tmp_num_nbrs] = tmp_selected_knrc
    
        return dp_feature_pair_an, dp_feature_pair_d, dp_feature_pair_rc
    
    
    
    
    
    
    
    def extract_feature_pair_embedding(
                        self,
                        center_atomic_number:int,
                        nbr_atomic_number:int,
                        rcut:float,
                        max_num_nbrs:int):
        '''
        Description
        -----------
            1. 计算 `deepmd feature pair embedding`
                - (1/Rij, xij/Rij^2, yij/Rij^2, zij/Rij^2)
            2. 
        
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
        
        '''
        ### Step 1. 调用 `self.extract_feature_pair()` 获取信息
        #       dp_feature_pair_an: (num_center, max_num_nbrs)
        #       dp_feature_pair_d : (num_center, max_num_nbrs)
        #       dp_feature_pair_rc : (num_center, max_num_nbrs, 3)
        dp_feature_pair_an, dp_feature_pair_d, dp_feature_pair_rc = \
                        self.extract_feature_pair(
                            center_atomic_number=center_atomic_number,
                            nbr_atomic_number=nbr_atomic_number,
                            rcut=rcut,
                            max_num_nbrs=max_num_nbrs)
        
        ### Step 2. 根据上述信息构造 deepmd embedding
        ### Step 2.1. Rij^2, shape = (4, 10)
        dp_feature_pair_Rij2_ = np.power(dp_feature_pair_d, 2)
        # shape = (4, 10, 1)
        dp_feature_pair_Rij2 = np.repeat(
                                    dp_feature_pair_Rij2_[:, :, np.newaxis],
                                    1,
                                    axis=2)
        
        ### Step 2.2. lasting three element of embedding (xij/Rij^2, yij/Rij^2, zij/Rij^2)
        # shape = (4, 10, 3)
        dp_feature_pair_xyz = dp_feature_pair_rc / dp_feature_pair_Rij2
        # shape = (4, 10, 3)
        dp_feature_pair_xyz = np.where(np.isnan(dp_feature_pair_xyz), 0, dp_feature_pair_xyz)
        
        ### Step 2.3. concatenate
        # shape = (4, 10)
        dp_feature_pair_Rij = dp_feature_pair_d
        # shape = (4, 10)
        dp_feature_pair_Rij_r = np.where(dp_feature_pair_d==0, 0, np.reciprocal(dp_feature_pair_d))
        # shape = (4, 10, 1)
        dp_feature_pair_Rij_r = np.repeat(
                                    dp_feature_pair_Rij_r[:, :, np.newaxis],
                                    1,
                                    axis=2)
        
        # shape = (4, 10, 4)    # (num_centers, max_num_nbrs, embedding_size)
        dp_feature_pair_embedding = np.concatenate([dp_feature_pair_Rij_r, dp_feature_pair_xyz], axis=2)
        
        return dp_feature_pair_embedding