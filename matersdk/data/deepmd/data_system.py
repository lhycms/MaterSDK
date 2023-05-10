import os
import shutil
import numpy as np
from typing import List
import multiprocessing as mp

from ...io.publicLayer.traj import Trajectory
from ...io.publicLayer.structure import DStructure
from ...io.publicLayer.neigh import StructureNeighborsDescriptor


class DeepmdDataSystem(object):
    '''
    Description
    -----------
        1. `DeepDataSystem` 中包含很多frame的DStructure对象 (所有frame的原子数相同!!!)
        2. 这个类用于....?
        
    Attributions
    ------------
        1. self.num_structures: int 
        2. self.structures_lst: List[Dstructure]
        3. self.total_energys_lst: List[float]
        4. self.kinetic_energys_lst: List[float]
        5. self.potential_energys_lst: List[float]
        6. self.virial_tensors_lst: List[np.ndarray]
        7. self.atomic_numbers_lst: List[int]
        8. self.num_atoms
    '''
    def __init__(
                self,
                structures_lst: List[DStructure],
                total_energys_array: np.ndarray,
                potential_energys_array: np.ndarray,
                kinetic_energys_array: np.ndarray,
                virial_tensors_array: np.ndarray,
                rcut:float):
        self.rcut = rcut
        
        self.num_structures = len(structures_lst)
        self.structures_lst = structures_lst
        self.total_energys_array = total_energys_array
        self.potential_energys_array = potential_energys_array
        self.kinetic_energys_array = kinetic_energys_array
        self.virial_tensors_array = virial_tensors_array
        
        self.atomic_numbers_lst = self._get_atomic_numbers()    # 不重复
        self.num_atoms = self._get_num_atoms()
    
    
    def _get_num_atoms(self):
        '''
        Description
        -----------
            1. 得到每个system中原子的数目
        '''
        num_atoms = self.structures_lst[0].num_sites
        return num_atoms
    
    def _get_atomic_numbers(self):
        '''
        Description
        -----------
            1. 得到体系内的原子序数
                - e.g. [3, 14]
        '''
        atomic_numbers_lst = [tmp_specie.Z for tmp_specie in self.structures_lst[0].species]
        atomic_numbers_lst = list(dict.fromkeys(atomic_numbers_lst))
        return atomic_numbers_lst
    
    
    @staticmethod
    def from_trajectory_s(
                trajectory_object:Trajectory,
                rcut:float,
                ):
        '''
        Description
        -----------
            1. 串行提取信息
        
        Parameters
        ----------
            1. trajectory_object: Trajectory
                - 轨迹对象
                - e.g. `matersdk.io.pwmat.output.movement.Movement`
        '''         
        ### Step 1. 得到 Movement 中所有构型的相关信息
        (structures_lst, 
        total_energys_array,
        potential_energys_array,
        kinetic_energys_array,
        virial_tensors_array) = \
                trajectory_object.get_all_frame_structures_info()
            
            
        ### Step 2. 初始化
        dp_data_system = DeepmdDataSystem(
                        structures_lst=structures_lst,
                        total_energys_array=total_energys_array,
                        potential_energys_array=potential_energys_array,
                        kinetic_energys_array=kinetic_energys_array,
                        virial_tensors_array=virial_tensors_array,
                        rcut=rcut
        )
        
        return dp_data_system
    
    
    def save(
            self,
            dir_path:str,
            scaling_matrix:List[int]=[3,3,3]):
        '''
        Description
        -----------
            1. 将初始化的信息存入对应image文件夹下 (以npy形式)。
            2. 存储近邻原子信息 -- `nbr_info.npy`
            
        Parameters
        ----------
            1. dir_path: str
                - 输出的文件夹路径
                - 各个frame的存储路径为 `<dir_path>/IMAGE_000`
        '''
        ### Step 1. 存储:
        #               1. box.npy
        #               2. coord.npy
        #               3. total_energy.npy
        #               4. kinetic_energy.npy
        #               5. potential_energy.npy
        #               6. atomic_energy.npy
        #               7. atomic_force.npy
        #               8. virial.npy
        num_bits = len(str(self.num_structures))
        for tmp_idx in range(self.num_structures):
            tmp_image_dir_path = os.path.join(dir_path, f"%0{num_bits}d" % tmp_idx)
            ### Step 1.1. 存在-->删除
            if os.path.exists(tmp_image_dir_path):
                shutil.rmtree(tmp_image_dir_path)
            # 创建
            os.mkdir(tmp_image_dir_path)
            
            ### Step 1.2. 保存相关信息
            self.structures_lst[tmp_idx].to(
                                    output_file_format="pwmat",
                                    output_file_path=os.path.join(tmp_image_dir_path, "atom.config"))
            # 1. box.npy
            np.save(
                    file=os.path.join(tmp_image_dir_path, "box.npy"),
                    arr=self.structures_lst[tmp_idx].lattice.matrix
            )
            # 2. coord.npy
            np.save(
                    file=os.path.join(tmp_image_dir_path, "coord.npy"),
                    arr=self.structures_lst[tmp_idx].cart_coords
            )
            # 3\4\5. energy.npy
            np.save(
                    file=os.path.join(tmp_image_dir_path, "energy.npy"), 
                    arr=np.array([
                            self.total_energys_array[tmp_idx],
                            self.kinetic_energys_array[tmp_idx],
                            self.potential_energys_array[tmp_idx]
                    ])
            )
            # 6. atomic_energy
            try:
                np.save(
                        file=os.path.join(tmp_image_dir_path, "atomic_energy.npy"),
                        arr=self.structures_lst[tmp_idx].get_atomic_energy()
                )
            except AttributeError:
                pass
            # 7. atomic_force
            np.save(
                    file=os.path.join(tmp_image_dir_path, "atomic_force.npy"),
                    arr=self.structures_lst[tmp_idx].get_atomic_force()
            )
            # 8. virial
            np.save(
                    file=os.path.join(tmp_image_dir_path, "virial.npy"),
                    arr=self.virial_tensors_array[tmp_idx]
            )
            # 9. atomic_number
            np.save(
                    file=os.path.join(tmp_image_dir_path, "atomic_number.npy"),
                    arr=np.array(self.atomic_numbers_lst)
            )
            # 10. num_atoms
            #num_atoms_dict = dict.fromkeys(self.atomic_numbers_lst, 0)
            #for tmp_atomic_number in self.atomic_numbers_lst:
            #    num_atoms_dict[tmp_atomic_number] += 1
            #np.save(
            #        file=os.path.join(tmp_image_dir_path, "num_atoms.npy"),
            #        arr=np.array(list(num_atoms_dict.values()))
            #)
        
        # 11. nbr_info.npy: 存储 nbr_info 的部分多进程并行
        parameters_lst = [(
                        os.path.join(dir_path, f"%0{num_bits}d" % tmp_idx),
                        self.structures_lst[tmp_idx],
                        scaling_matrix,
                        self.rcut) for tmp_idx in range(self.num_structures)]
        
        with mp.Pool(os.cpu_count()-2) as pool:
            pool.starmap(ParallelFunction.save_struct_nbr, parameters_lst)
            
            

class ParallelFunction(object):
    '''
    Description
    -----------
        1. 一些需要多进程并行的函数
    '''    
    @staticmethod
    def save_struct_nbr(
                    tmp_image_dir_path:int,
                    structure:DStructure,
                    scaling_matrix:List[int],
                    rcut:float):
        '''
        Description
        -----------
            1. 存储 `struct_nbr` 的信息，在 `DeepmdDataSystem.save()` 中调用

        Parameters
        ----------
            1. tmp_image_dir_path: str
                - Image文件夹的路径
            2. structure: DStructure
                - 结构
            3. scaling_matrix: List[int]
                - 扩包倍数
        '''
        # 10. nbr_info.npy
        struct_nbr = StructureNeighborsDescriptor.create(
                    'v3',
                    structure=structure,
                    scaling_matrix=scaling_matrix,
                    reformat_mark=True,
                    coords_are_cartesian=True,
                    rcut=rcut,
        )

        np.save(
                file=os.path.join(tmp_image_dir_path, "nbrs_atomic_numbers.npy"),
                arr=struct_nbr.key_nbr_atomic_numbers
        )
        np.save(
                file=os.path.join(tmp_image_dir_path, "nbrs_distances.npy"),
                arr=struct_nbr.key_nbr_distances
        )
        np.save(
                file=os.path.join(tmp_image_dir_path, "nbrs_coords.npy"),
                arr=struct_nbr.key_nbr_coords
        )