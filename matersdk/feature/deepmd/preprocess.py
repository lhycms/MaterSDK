import numpy as np
from typing import Union, List, Dict

from ...data.deepmd.data_system import DpLabeledSystem
from ...io.publicLayer.neigh import StructureNeighborsDescriptor
from ...feature.deepmd.se_pair import DpseTildeRPairDescriptor


class TildeRPairNormalizer(object):
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
        if (davg != False) and (dstd != False):
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
            2. 
        
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
            1. tildeRs_array: np.ndarray
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
                max_num_nbrs_dict:Dict[int, int],
                scaling_matrix:List[int]):
        '''
        Description
        -----------
            1. ...
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
            for tmp_nbr_an in nbr_atomic_numbers:
                # e.g. Li-Li : (48, 100, 4) -- (num_centers, max_num_nbrs ,4)
                # e.g. Li-Si : (48, 80,  4) -- (num_centers, max_num_nbrs, 4)
                tmp_nbr_tildeR = DpseTildeRPairDescriptor.create(
                            'v1',
                            structure_neighbors=tmp_struct_nbr,
                            center_atomic_number=center_atomic_number,
                            nbr_atomic_number=tmp_nbr_an,
                            rcut=rcut,
                            rcut_smooth=rcut_smooth).get_tildeR(
                                    max_num_nbrs=max_num_nbrs_dict[tmp_nbr_an])
                tmp_all_nbrs_tildeRs_lst.append(tmp_nbr_tildeR)
            # tmp_tildeR: Li-Li&Si
            # shape = (48, 180, 4)
            tmp_structure_tildeR = np.concatenate(tmp_all_nbrs_tildeRs_lst, axis=1)
        
            all_structures_tildeRs_lst.append(tmp_structure_tildeR)
            
        # tildeR_tot: all structures for Li-Li/Si
        # shape = (48, 1800, 4)
        tildeRs_array = np.concatenate(all_structures_tildeRs_lst, axis=1)
        
        return tildeRs_array