import numpy as np
from typing import List

from ...io.pwmat.output.movement import Movement
from ...io.publicLayer.neigh import (
                                StructureNeighborsDescriptor,
                                StructureNeighborsV1
                                )
from ...io.pwmat.utils.parameters import specie2atomic_number



class AvgBond(object):
    def __init__(self,
                movement_path:str,
                element_1:str,
                element_2:str,
                cutoff:float):
        '''
        Parameters
        ----------
            1. movement_path: str
                - MOVEMENT 的路径
            2. element_1: str
                - 成键的第一种元素
            3. element_2:
                - 成键的第二种元素
            4. cutoff:
                - 成键的标准
        '''
        self.frames_lst = Movement(movement_path=movement_path).get_all_frame_structures_info()[0]
        self.atomic_number_1 = specie2atomic_number[element_1]
        self.atomic_number_2 = specie2atomic_number[element_2]
        self.cutoff = cutoff
        self.num_sites = self.frames_lst[0].num_sites
        
    
    def get_frames_avg_bond(
                    self,
                    scaling_matrix:List[int]=[3, 3, 3],
                    reformat_mark:bool=True,
                    coords_are_cartesian:bool=True):
        '''
        Description
        -----------
            1. 计算 Movement 中所有 frame 的平均键长
        
        Parameters
        ----------
            1. scaling_matrix: List[int]
                - bulk: [3, 3, 3]
                - slab: [3, 3, 1]
            2. reformat_mark: bool
                - 是否按原子序数从小到大排序
            3. coords_are_cartesian: bool
                - 是否按照笛卡尔坐标计算
        '''
        avg_bond_lengths_lst = []
        for tmp_idx, tmp_struct in enumerate(self.frames_lst):
            struct_neigh = StructureNeighborsDescriptor.create(
                                'v1',
                                structure=tmp_struct,
                                rcut=self.cutoff,
                                scaling_matrix=scaling_matrix,
                                reformat_mark=reformat_mark,
                                coords_are_cartesian=coords_are_cartesian)
            tmp_avg_bond_length = self._get_avg_bond_length(
                                        struct_neigh=struct_neigh)
            avg_bond_lengths_lst.append(tmp_avg_bond_length)

        avg_bond_lengths_array = np.array(avg_bond_lengths_lst)
        
        return avg_bond_lengths_array
    
    
    def _get_avg_bond_length(
                    self,
                    struct_neigh:StructureNeighborsV1):
        '''
        Description
        -----------
            1. 根据 `StructureNeighborVX` 的信息计算平均键长 (e.g. 所有 `Ge-Te` 对的长度)

        Parameters
        ----------
            1. struct_neigh: StructureNeighborsV#
                - 近邻原子信息
        '''
        key_nbr_ans = struct_neigh.key_nbr_atomic_numbers
        key_nbr_distances = struct_neigh.key_nbr_distances
        
        ### Step 1. 按照中心原子种类筛选
        filter_center = np.where(
                        key_nbr_ans[:, 0] == self.atomic_number_1,
                        True,
                        False)
        filter_center = np.repeat(
                        filter_center[:, np.newaxis],
                        key_nbr_ans.shape[1],
                        axis=1)
        
        ### Step 2. 按照近邻原子种类筛选
        filter_nbr = np.where(
                        key_nbr_ans == self.atomic_number_2,
                        True,
                        False)
     
        ### Step 3. 取 Step_1 和 Step_2 的and
        filter_tot = filter_center & filter_nbr
        ### Note: 中心原子全部取 False !!!
        filter_tot[:, 0] = False

        ### Step 4. 计算有效的键长，并存储为 np.array 格式
        effective_bonds_array = np.where(
                            filter_tot,
                            key_nbr_distances,
                            0)
        
        ### Step 5. 计算键长之和与键的数目
        num_bonds = np.count_nonzero(effective_bonds_array)
        sumlength_bonds = np.sum(effective_bonds_array)
        
        return sumlength_bonds / num_bonds



class PairBond(object):
    def __init__(
                self,
                movement_path:str,
                atom1_idx:int,
                atom2_idx:int,
                ):
        '''
        Description
        -----------
            1. 
        
        Parameters
        ----------
            1. movement_path: str
            2. atom1_idx: int 
            3. atom2_idx: int
        '''
        self.frames_lst = Movement(movement_path=movement_path).get_all_frame_structures_info()[0]
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx
    
    
    def get_frames_pair_bond(self):
        pair_bond_lengths_lst = []
        for tmp_struct in self.frames_lst:
            tmp_coords = tmp_struct.cart_coords
            tmp_atom1_coord = tmp_coords[self.atom1_idx]
            tmp_atom2_coord = tmp_coords[self.atom2_idx]
            tmp_pair_bond = self._get_pair_bond_length(tmp_atom1_coord, tmp_atom2_coord)
            pair_bond_lengths_lst.append(tmp_pair_bond)
        
        return pair_bond_lengths_lst
            
    
    
    def _get_pair_bond_length(
                            self,
                            atom1_coord:np.ndarray,
                            atom2_coord:np.ndarray
                            ):
        '''
        Description
        -----------
            1. 计算两个坐标之间的距离
        '''
        return np.linalg.norm(atom2_coord - atom1_coord)