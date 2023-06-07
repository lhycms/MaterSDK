import numpy as np
import h5py
from typing import Union, List, Dict

from ...data.deepmd.data_system import DpLabeledSystem
from ...io.publicLayer.structure import DStructure
from ...io.publicLayer.neigh import StructureNeighborsDescriptor
from ...feature.deepmd.se_pair import DpseTildeRPairDescriptor


class TildeRNormalizer(object):
    '''
    Description
    -----------
        1. 计算 $\tilde{R}$ + normalize 一条龙
        2. 可以设置中心原子为不同的元素
            - e.g. 计算下列 pair 的 `davg`, `dstd`
                1) Li-Li/Si 的 $\tilde{R}$
                2) Si-Li/Si 的$\tilde{R}$
        3. 在 `TildeRNormalizer` 内部会调用：
            1. `TildeRPairNormalizer`
                - 
            2. `NormalizerPremise`
                - 
    '''
    def __init__(
            self,
            rcut:float,
            rcut_smooth:float,
            center_atomic_numbers:Union[List[int], np.ndarray],
            nbr_atomic_numbers:Union[List[int], np.ndarray],
            max_num_nbrs:List[int],
            scaling_matrix:List[int],
            dp_labeled_system:Union[DpLabeledSystem, bool]=False,
            structure_indices:Union[List[int], bool]=False,
            davgs:Union[np.ndarray, bool]=False,
            dstds:Union[np.ndarray, bool]=False
            ):
        '''
        Parameters
        ----------
            1. rcut: float
                - DeepPot-SE 的截断半径
            2. rcut_smooth: float
                - DeepPot-SE 的 smooth cutoff
            3. center_atomic_numbers: List[int]
                - 中心原子的原子序数，按照从小到大排序
            4. nbr_atomic_numbers: List[int]
                - 近邻原子的原子序数，按照从小到大排序
            5. max_num_nbrs: List[int]
                - 中心原子近邻的nbr元素最多的原子数目
                - Note: 与 nbr_atomic_numbers 一一对应
            6. scaling_matrix: List[int]
                - 
            7. dp_labeled_system: Union[DpLabeledSystem, bool]
                - Way 1
            8. structure_indices: Union[List[int], bool]
                - Way 1
            9. davgs: Union[np.ndarray, bool]
                - Way 2
            10. dstds: Union[np.ndarray, bool]
                - Way 2
            
        Note
        ----    
            1. 有两种方法初始化 `TildeRNormalizer`
                1) dp_labeled_system + structure_indices
                2) davgs + dstds
        '''
        self.rcut = rcut
        self.rcut_smooth = rcut_smooth
        self.center_atomic_numbers = center_atomic_numbers
        self.nbr_atomic_numbers = nbr_atomic_numbers
        self.max_num_nbrs = max_num_nbrs
        self.scaling_matrix = scaling_matrix
        
        ### Step 1. Get the `davgs` and `dstds`
        if (dp_labeled_system is not False) and (structure_indices is not False):
            self.davgs, self.dstds = self.calc_stats(
                            dp_labeled_system=dp_labeled_system,
                            structure_indices=structure_indices
            )
        
        elif (davgs is not False) and (dstds is not False):
            self.davgs = davgs
            self.dstds = dstds
            
            assert (davgs.shape[0] == len(self.center_atomic_numbers))
            assert (davgs.shape[1] == 4)
            assert (dstds.shape[0] == len(self.center_atomic_numbers))
            assert (dstds.shape[1] == 4)
        
        else:
            raise ValueError("You must specify the davgs and dstds!")
    
    
    def __str__(self):
        return self.__repr__()
    
    
    def __repr__(self):
        print("{0:*^80s}".format(" TildeRNormalizer Summary "))
        
        print("\t * {0:26s}: {1:14f}".format("rcut", self.rcut))
        print("\t * {0:26s}: {1:14f}".format("rcut_smooth", self.rcut_smooth))
        print("\t * {0:26s}:\t".format("center_atomic_numbers:"), self.center_atomic_numbers)
        print("\t * {0:26s}:\t".format("nbr_atomic_numbers:"), self.nbr_atomic_numbers)
        print("\t * {0:26s}:\t".format("max_num_nbrs"), self.max_num_nbrs)
        print("\t * {0:26s}:\t".format("scaling_matrix"), self.scaling_matrix)
        print("\t * {0:26s}:\t".format("davgs"))
        print(self.davgs)
        print("\t * {0:26s}:\t".format("dstds"))
        print(self.dstds)        
        
        print("{0:*^80s}".format("**"))
        return ""
    
    
    def calc_stats(
                self, 
                dp_labeled_system:DpLabeledSystem,
                structure_indices:List[int]
                ):
        '''
        Description
        -----------
            1. 如果初始化的时候使用 `dp_labeled_system` 和 `structure_indices`，则会调用这个函数
        
        Parameters
        ----------
            1. dp_labeled_system: DpLabeledSystem
                - 
            2. structure_indices: List[int]
                - 
        
        Return
        ------
            1. davgs: np.ndarray
                - .shape = (num_types, 4)
            2. dstds: np.ndarray
                - .shape = (num_types, 4)
        '''
        ### Step 1. 
        #       e.g. {
        #            3: (48, 1800, 4)
        #           14: (24, 1800, 4)
        #}
        davgs_lst = []
        dstds_lst = []
        for tmp_center_an in self.center_atomic_numbers:
            ### Step 1.1. Calcuate `tildeRs_array`
            # e.g.
            #   Li .shape = (48, 1800, 4)
            #   Si .shape = (24, 1800, 4)
            tmp_tildeRs_array = NormalizerPremise.concat_tildeRs(
                    dp_labeled_system=dp_labeled_system,
                    structure_indices=structure_indices,
                    rcut=self.rcut,
                    rcut_smooth=self.rcut_smooth,
                    center_atomic_number=tmp_center_an,
                    nbr_atomic_numbers=self.nbr_atomic_numbers,
                    max_num_nbrs=self.max_num_nbrs,
                    scaling_matrix=self.scaling_matrix
            )
            
            ### Step 1.2. Calculate `davg`, `dstd`
            tmp_normalizer = TildeRPairNormalizer(tildeRs_array=tmp_tildeRs_array)
            # tmp_normalizer.davg.shape = (1, 4)
            # tmp_normalizer.dstd.shape = (1, 4)
            davgs_lst.append(tmp_normalizer.davg)
            dstds_lst.append(tmp_normalizer.dstd)
        
        ### Step 2.
        # shape = (num_types, 4)
        davgs = np.concatenate(davgs_lst, axis=0)
        # shape = (num_types, 4)
        dstds = np.concatenate(dstds_lst, axis=0)
        
        return davgs, dstds
    
    
    def normalize(self, structure:DStructure):
        '''
        Description
        -----------
            1.
        
        Parameters
        ----------
            1. structure: DStructure
                - 
        
        Return
        ------
            1. tildeR_dict: Dict[str, np.ndarray]
                - e.g. {
                        "3_3": np.ndarray,
                        "3_14": np.ndarray,
                        "14_3": np.ndarray,
                        "14_14": np.ndarray
                    }
        '''
        tildeR_dict:Dict[str, np.ndarray] = {}
        ### Step 1. 计算新结构的 `StructureNeighbors`
        struct_neigh = StructureNeighborsDescriptor.create(
                    'v1',
                    structure=structure,
                    rcut=self.rcut,
                    scaling_matrix=self.scaling_matrix,
                    reformat_mark=True,
                    coords_are_cartesian=True)
        
        ### Step 2. 计算 Normalized 的 $\tilde{R}$
        for tmp_idx_center_an, tmp_center_an in enumerate(self.center_atomic_numbers):
            for tmp_idx_nbr_an, tmp_nbr_an in enumerate(self.nbr_atomic_numbers):
                ### Step 2.1. 计算 Environment Matrix -- $\tilde{R}$ 
                tmp_tilde_R = DpseTildeRPairDescriptor.create(
                            'v1',
                            structure_neighbors=struct_neigh,
                            center_atomic_number=tmp_center_an,
                            nbr_atomic_number=tmp_nbr_an,
                            rcut=self.rcut,
                            rcut_smooth=self.rcut_smooth
                ).get_tildeR(max_num_nbrs=self.max_num_nbrs[tmp_idx_nbr_an])
                
                ### Step 2.2. Normalize Environment Matrix -- $\tilde{R}$
                tmp_davg = self.davgs[tmp_idx_center_an]
                tmp_dstd = self.dstds[tmp_idx_center_an]
                tmp_tilde_r_pair_normalizer = TildeRPairNormalizer(
                            davg=tmp_davg,
                            dstd=tmp_dstd
                )
                tmp_normalized_tildeR = tmp_tilde_r_pair_normalizer.normalize(tildeRs_array=tmp_tilde_R)
                tildeR_dict.update({"{0}_{1}".format(tmp_center_an, tmp_nbr_an): tmp_normalized_tildeR})
        
        return tildeR_dict


    def to(self, hdf5_file_path:str):
        h5_file = h5py.File(hdf5_file_path, 'w')
        
        h5_file.create_dataset("rcut", data=self.rcut)
        h5_file.create_dataset("rcut_smooth", data=self.rcut_smooth)
        h5_file.create_dataset("center_atomic_numbers", data=self.center_atomic_numbers)
        h5_file.create_dataset("nbr_atomic_numbers", data=self.nbr_atomic_numbers)
        h5_file.create_dataset("max_num_nbrs", data=self.max_num_nbrs)
        h5_file.create_dataset("scaling_matrix", data=np.array(self.scaling_matrix))
        h5_file.create_dataset("davgs", data=self.davgs)
        h5_file.create_dataset("dstds", data=self.dstds)
        
        h5_file.close()
        
        
    @classmethod
    def from_file(cls, hdf5_file_path:str):
        ### Step 1. Extract information from hdf5 file
        hdf5_file = h5py.File(hdf5_file_path, 'r')
        rcut = hdf5_file["rcut"][()]
        rcut_smooth = hdf5_file["rcut_smooth"][()]
        center_atomic_numbers = hdf5_file["center_atomic_numbers"][()]
        nbr_atomic_numbers = hdf5_file["nbr_atomic_numbers"][()]
        max_num_nbrs = hdf5_file["max_num_nbrs"][()]
        scaling_matrix = hdf5_file["scaling_matrix"][()]
        davgs = hdf5_file["davgs"][()]
        dstds = hdf5_file["dstds"][()]
        hdf5_file.close()
        
        ### Step 2. Initialize the `TildeRNormailzer`
        tilde_r_normalizer = cls(
                        rcut=rcut,
                        rcut_smooth=rcut_smooth,
                        center_atomic_numbers=center_atomic_numbers,
                        nbr_atomic_numbers=nbr_atomic_numbers,
                        max_num_nbrs=max_num_nbrs,
                        scaling_matrix=scaling_matrix,
                        davgs=davgs,
                        dstds=dstds
        )
        
        return tilde_r_normalizer
        
        
        


class TildeRPairNormalizer(object):
    '''
    Description
    -----------
        1. 中心原子确定
            - e.g. 计算下列 pair 的 `davg`, `dstd`
                1) Li-Li/Si 的 $\tilde{R}$
                2) Si-Li/Si 的 $\tilde{R}$
                
    '''
    def __init__(self,
                tildeRs_array:Union[np.ndarray, bool]=False,
                davg:Union[np.ndarray, bool]=False,
                dstd:Union[np.ndarray, bool]=False):
        '''
        Description
        -----------
            1. Calculate the `davg` and `dstd` of Environment matrix ($\widetilde{R}^i = (s, sx/r, sy/r, sz/r)$)
            2. Normalize the Environment matrix ($\widetilde{R}^i = (s, sx/r, sy/r, sz/r)$)

        Parameters
        ----------
            1. tildeRs_array: 
                - Note:
                    1. You should calculate tildeR of `Li-Li` and `Li-Si`, and concat them, then calculate `davg` and `dstd`
            2. davg: np.ndarray 
                - `davg.shape = (1, 4)`
            3. dstd: np.ndarray
                - `dstd.shape = (1, 4)`
        '''
        # shape = (1, 4)
        if (davg is not False) and (dstd is not False):
            self.davg = davg
            self.dstd = dstd
        elif (tildeRs_array is not False):
            # shape: (num_frames, num_centers, max_num_nbrs, 4)     e.g. (48, 26, 4)
            #   ->
            # shape: (num_frames * num_centers * max_num_nbrs, 4)   e.g. (1248, 4)
            self.davg, self.dstd = self.calc_stats( tildeRs_array.reshape(-1, 4) )
        else:
            raise ValueError("You should check davg and dstd!")
    
    
    def calc_stats(self, tildeRs_array:np.ndarray):
        '''
        Description
        -----------
            1. 计算 DeepPot-SE 中 TildeR 的平均值(`avg`)和方差(`std`)
            2. 如果初始化的时候使用 `dp_labeled_system` 和 `structure_indices`，则会调用这个函数
        
        Parameters
        ----------
            1. tildeRs_array: np.ndarray
                - .shape = (num_frames * num_centers * max_num_nbrs, 4)
        '''
        ### Step 1. 分别获取径向信息(`info_radius`)和角度信息(`info_angles`)
        # shape: (num_frames * num_centers * max_num_nbrs, 1)
        info_radius = tildeRs_array[:, 0].reshape(-1, 1)
        # shape: (num_frames * num_centers * max_num_nbrs, 3)
        info_angles = tildeRs_array[:, 1:].reshape(-1, 3)
        
        ### Step 2. 分别获取径向信息(`info_radius`)和角度信息(`info_angles`)的一些量:
        #     1) sum: 求和
        #     2) sum^2: 平方和 -- 先平方后求和
        #     3) total_num_pairs: num_centers * max_num_nbrs
        ### Step 2.1. sum: 求和
        sum_info_radius = np.sum(info_radius)
        sum_info_angles = np.sum(info_angles) / 3.0
        
        ### Step 2.2. sum^2: 平方和 -- 先平方后求和
        sum2_info_radius = np.sum(
                    np.multiply(info_radius, info_radius)
        )
        sum2_info_angles = np.sum(
                    np.multiply(info_angles, info_angles)
        ) / 3.0
        
        ### Step 2.3. total_num_pairs.shape: num_centers * max_num_nbrs
        total_num_pairs = info_radius.shape[0] # Error: np.count_nonzero(info_radius.flatten() != 0.)

        
        ### Step 3. 计算平均值 -- davg_unit
        davg_unit = [sum_info_radius / (total_num_pairs + 1e-15), 0, 0, 0]
        # shape = (1, 4)
        davg_unit = np.array(davg_unit).reshape(-1, 4)
        
        
        ### Step 4. 计算方差 -- dstd_unit
        dstd_unit = [
            self._calc_std(sum2_value=sum2_info_radius, sum_value=sum_info_radius, N=total_num_pairs),
            self._calc_std(sum2_value=sum2_info_angles, sum_value=sum_info_angles, N=total_num_pairs),
            self._calc_std(sum2_value=sum2_info_angles, sum_value=sum_info_angles, N=total_num_pairs),
            self._calc_std(sum2_value=sum2_info_angles, sum_value=sum_info_angles, N=total_num_pairs)
        ]
        # shape = (1, 4)
        dstd_unit = np.array(dstd_unit).reshape(-1, 4)
        
        return davg_unit, dstd_unit
        
    
    def _calc_std(self, sum2_value:float, sum_value:float, N:int):
        '''
        Description
        -----------
            1. 计算标准差
        
        Parameters
        ----------
            1. sum2_value: float
                - sum2_value = \sum_i^N{x_i^2}，先平方后求和
            2. sum_value: float
                - sum_value  = \sum_i^N{x_i}
        '''
        if (N == 0):
            return 1e-2
        std = np.sqrt(
                sum2_value / N - np.multiply(sum_value/N, sum_value/N)
        )
        if np.abs(std) < 1e-2:
            std = 1e-2
        return std
        
    
    def normalize(self, tildeRs_array:np.ndarray):
        '''
        Description
        -----------
            1. 
        
        Parameters
        ----------
            1. tildeR_array: np.ndarray
                - shape = (num_frames, num_centers, max_num_nbrs, 4)
        
        Note
        ----    
            1. You can input environment matrix for `single frame` or `many frames`
                - single frame: .shape = (num_centers, max_num_nbrs, 4) 
                - many frames : .shape = (num_frames, num_centers, max_num_nbrs, 4)
        '''
        if (len(tildeRs_array.shape) == 3):
            # single frame: tildeR.shape = (num_centers, max_num_nbrs, 4)
            davg = self.davg.reshape(1, 1, 4)
            dstd = self.dstd.reshape(1, 1, 4)
        elif (len(tildeRs_array.shape) == 4):
            # many frames: tildeR.shape = (num_frames, num_centers, max_num_nbrs, 4)
            davg = self.davg.reshape(1, 1, 1, 4)
            dstd = self.dstd.reshape(1, 1, 1, 4)
        
        ### Normalize the environment matrix
        result = (tildeRs_array - davg) / dstd
        
        return result
    


class NormalizerPremise(object):
    @staticmethod
    def concat_tildeRs(
                dp_labeled_system:DpLabeledSystem,
                structure_indices:List[int],
                rcut:float, 
                rcut_smooth:float,
                center_atomic_number:int,
                nbr_atomic_numbers:List[int],
                max_num_nbrs:List[int],
                scaling_matrix:List[int]):
        '''
        Description
        -----------
            1. 中心原子的元素种类是确定的！！！
                - e.g. 计算多个结构的 Li-Li, Li-Si 的 $\tilde{R}$ 并合并
                    - Li-Li: (48, 100, 4) -- (num_centers, max_num_nbrs, 4)
                    - Li-Si: (48, 80, 4)  -- (num_centers, max_num_nbrs, 4)
                    - 取10个结构并计算`avg`和`std`。最终的 concated_tildeRs.shape = (48, 1800, 4)
        
        Parameters
        ----------
            1. 
        '''
        ### Step 1. 得到 DStructure 的列表
        structures_lst = [
            dp_labeled_system.structures_lst[tmp_idx] for tmp_idx in structure_indices
        ]
        
        ### Step 2. 
        all_structures_tildeRs_lst = []    # 所有结构的 $\tilde{R}$
        for tmp_structure in structures_lst:
            tmp_struct_nbr = StructureNeighborsDescriptor.create(
                            'v1',
                            structure=tmp_structure,
                            rcut=rcut,
                            scaling_matrix=scaling_matrix,
                            reformat_mark=True,
                            coords_are_cartesian=True)
            
            tmp_all_nbrs_tildeRs_lst = []    # 同一结构，同一中心原子，不同近邻原子
            for tmp_idx_nbr_an, tmp_nbr_an in enumerate(nbr_atomic_numbers):
                # e.g. Li-Li : (48, 100, 4) -- (num_centers, max_num_nbrs ,4)
                # e.g. Li-Si : (48, 80,  4) -- (num_centers, max_num_nbrs, 4)
                tmp_nbr_tildeR = DpseTildeRPairDescriptor.create(
                            'v1',
                            structure_neighbors=tmp_struct_nbr,
                            center_atomic_number=center_atomic_number,
                            nbr_atomic_number=tmp_nbr_an,
                            rcut=rcut,
                            rcut_smooth=rcut_smooth).get_tildeR(
                                    max_num_nbrs=max_num_nbrs[tmp_idx_nbr_an])
                tmp_all_nbrs_tildeRs_lst.append(tmp_nbr_tildeR)
            # tmp_tildeR: Li-Li&Si
            # shape = (48, 180, 4)
            tmp_structure_tildeR = np.concatenate(tmp_all_nbrs_tildeRs_lst, axis=1)
        
            all_structures_tildeRs_lst.append(tmp_structure_tildeR)
            
        # tildeR_tot: all structures for Li-Li/Si
        # shape = (48, 1800, 4)
        tildeRs_array = np.concatenate(all_structures_tildeRs_lst, axis=1)
        
        return tildeRs_array