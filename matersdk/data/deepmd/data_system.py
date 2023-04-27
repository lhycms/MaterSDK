import numpy as np
from typing import List
from ...io.publicLayer.traj import Trajectory
from ...io.publicLayer.structure import DStructure


class DeepmdDataSystem(object):
    '''
    Description
    -----------
        1. `DeepDataSystem` 中包含很多frame的DStructure对象
        2. 这个类用于....?
        
    Attributions
    ------------
        1. self.num_structures: int 
        2. self.structures_lst: List[Dstructure]
        3. self.total_energys_lst: List[float]
        4. self.kinetic_energys_lst: List[float]
        5. self.potential_energys_lst: List[float]
        6. self.virial_tensors_lst: List[np.ndarray]
    '''
    def __init__(
                self,
                structures_lst: List[DStructure],
                total_energys_array: np.ndarray,
                kinetic_energys_array: np.ndarray,
                potential_energys_array: np.ndarray,
                virial_tensors_array: np.ndarray
                ):
        self.num_structures = len(structures_lst)
        self.structures_lst = structures_lst
        self.total_energys_array = total_energys_array
        self.kinetic_energys_array = kinetic_energys_array
        self.potential_energys_array = potential_energys_array
        self.virial_tensors_array = virial_tensors_array
    
    
    @staticmethod
    def from_trajectory(
                trajectory_object:Trajectory,
                ):
        '''
        Description
        -----------
            1. trajectory_object: Trajectory
                - 轨迹对象
                - e.g. `matersdk.io.pwmat.output.movement.Movement`
        '''
        ### Step 0. `DeepmdDataSystem` 所含的构型数目
        num_structures = len(trajectory_object.get_chunksize())
         
        ### Step 1. 设置 `structures_lst`
        structures_lst = []
        for tmp_idx in range(num_structures):
            tmp_structure = trajectory_object.get_frame_structure(idx_frame=tmp_idx)
            structures_lst.append(tmp_structure)
        
        ### Step 2. 设置 `total_energys_lst`
        total_energys_lst = []
        kinetic_energys_lst = []
        potential_energys_lst = []
        for tmp_idx in range(num_structures):
            tmp_total_energy, tmp_kinetic_energy, tmp_potential_energy = \
                                            trajectory_object.get_frame_energy(idx_frame=tmp_idx)
            total_energys_lst.append(tmp_total_energy)
            kinetic_energys_lst.append(tmp_kinetic_energy)
            potential_energys_lst.append(tmp_potential_energy)
        
        ### Step 3. 设置 `virial_tensor`
        virial_tensors_lst = []
        for tmp_idx in range(num_structures):
            tmp_virial_tensor = trajectory_object.get_frame_virial(idx_frame=tmp_idx)
            virial_tensors_lst.append(tmp_virial_tensor)
            
        ### Step 4. 初始化
        dp_data_system = DeepmdDataSystem(
                        structures_lst=structures_lst,
                        total_energys_array=np.array(total_energys_lst),
                        kinetic_energys_array=np.array(kinetic_energys_lst),
                        potential_energys_array=np.array(potential_energys_lst),
                        virial_tensors_array=np.array(virial_tensors_lst)
        )
        
        return dp_data_system