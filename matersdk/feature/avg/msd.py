import os
import numpy as np
from typing import List
import multiprocessing as mp

from ...io.publicLayer.structure import DStructure
from ...io.publicLayer.traj import Trajectory
from ...io.pwmat.output.movement import Movement



class Msd(object):
    def __init__(self, trajectory:Trajectory):
        self.structures_lst = trajectory.get_all_frame_structures()
    
    
    def calc_msd(self):
        '''
        Description
        -----------
            1. Parallel the `MsdParallelFunction.calc_msd_s(structure_1, structure_2)`
        '''
        parameters_lst:List[List] = []
        for tmp_frame_idx in range(1, len(self.structures_lst)):
            parameters_lst.append([
                        self.structures_lst[tmp_frame_idx-1],
                        self.structures_lst[tmp_frame_idx]]
            )
        
        msd_values_lst:List[int] = []
        with mp.Pool(processes=os.cpu_count()-2) as pool:
            msd_values_lst = pool.starmap(
                                MsdParallelFunction.calc_msd_s,
                                parameters_lst
            )
        return msd_values_lst
            
        
class MsdParallelFunction(object):
    @staticmethod
    def calc_msd_s(
                structure_1:DStructure,
                structure_2:DStructure):
        '''
        Description
        -----------
            1. Calculate `Mean squared displacement` -- `MSD`
                - MSD = 1/n sum_{i=1}^n [x_i(t) - x_origin]^2
        '''
        relative_coords = structure_2.cart_coords - structure_1.cart_coords
        num_atoms = relative_coords.shape[0]
        
        return np.sum(np.power(relative_coords, 2)) / num_atoms
