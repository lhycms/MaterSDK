import linecache
import numpy as np
from abc import ABC, abstractmethod

from .lineLocator import LineLocator
from ...publicLayer.atom import Atom
from .parameters import atomic_number2specie



class AtomConfigExtractorBase(ABC):
    @abstractmethod
    def get_num_atoms(self):
        pass
    
    @abstractmethod
    def get_basis_vectors_lst(self):
        pass
    
    @abstractmethod
    def get_virial_tensor(self):
        pass
    
    @abstractmethod
    def get_atomic_numbers_lst(self):
        pass
    
    @abstractmethod
    def get_coords_lst(self):
        pass
    
    @abstractmethod
    def get_atomic_forces_lst(self):
        pass
    
    @abstractmethod
    def get_atomic_velocitys_lst(self):
        pass
    
    @abstractmethod
    def get_atomic_energys_lst(self):
        pass
    
    @abstractmethod
    def get_magnetic_moments(self):
        pass
    
    #@abstractmethod
    def get_atoms_lst(self):
        pass
    
    
    
class AtomConfigExtractor(AtomConfigExtractorBase):
    '''
    Description
    -----------
        1. 提取 atom.config 文件中的信息：
            1. num_atoms: int
                体系内的原子总数
            2. basis_vectors_array: list
                体系的基矢 (二维list)
            3. species_array: list of str
                各个 site 的原子种类
            4. coords_array: np.array
                各个 site 的分数坐标
    
        2. 当提取的是 MOVEMENT 中某一帧信息时，还可以提取到：
            1. `Force (eV/Angstrom)`
            2. `Velocity (bohr/fs)`
            3. `Atomic-Energy`
            4. ...
    '''
    def __init__(self,
                atom_config_path: str):
        '''
        Parameters
        ----------
            1. atom_config_path: str
                The absolute path of `atom.config` file
        '''
        self.atom_config_path = atom_config_path
        self.num_atoms = self.get_num_atoms()
        self.basis_vectors_array = self.get_basis_vectors_lst()
        self.species_array = [atomic_number2specie[atomic_number] for atomic_number in self.get_atomic_numbers_lst() ]
        self.coords_array = self.get_coords_lst()
        self.magnetic_moments = self.get_magnetic_moments()


    def get_num_atoms(self):
        '''
        Description
        -----------
            1. 得到体系的总原子数目
        '''
        # atom.config 文件第一行综合描述这个体系
        first_row = linecache.getline(self.atom_config_path, 1)
        num_atoms = int( first_row.split()[0] )
        return num_atoms


    def get_basis_vectors_lst(self):
        '''
        Description
        -----------
            1. 得到材料的基矢

        Return
        ------
            1. basis_vectors_lst: list, 是一个二维列表
            e.g. [ [14.376087, 0.0, 0.0],
                    [0.0, 12.357182, 0.0], 
                    [0.0, 0.0, 10.594184]
                ]
            
        '''
        basis_vectors_lst = []

        ### Step 1. 得到所有原子的原子序数、坐标
        content = "LATTICE"    # 此处需要大写
        idx_row = LineLocator.locate_all_lines(
                                    file_path=self.atom_config_path,
                                    content=content)[0]
        with open(self.atom_config_path, 'r') as f:
            atom_config_content = f.readlines()

        ### Step 2. 获取基矢向量
        for row_idx in [idx_row+1, idx_row+2, idx_row+3]:
            row_content = linecache.getline(self.atom_config_path, row_idx).split()[:3]
            
            single_direction_vector = [float(value) for value in row_content]
            basis_vectors_lst.append(single_direction_vector)
                
        return np.array(basis_vectors_lst)
    
    
    def get_virial_tensor(self):
        '''
        Description
        -----------
            1. 得到材料的维里张量 (virial tensor)

        Return
        ------
            1. virial_tensor: np.array, 是一个二维 np.ndarray            
        '''
        virial_tensor = []

        ### Step 1. 得到所有原子的原子序数、坐标
        content = "LATTICE"    # 此处需要大写
        idx_row = LineLocator.locate_all_lines(
                                    file_path=self.atom_config_path,
                                    content=content)[0]
        with open(self.atom_config_path, 'r') as f:
            atom_config_content = f.readlines()

        ### Step 2. 获取基矢向量
        for row_idx in [idx_row+1, idx_row+2, idx_row+3]:
            row_content = linecache.getline(self.atom_config_path, row_idx).split()[5:]
            single_direction_vector = [float(value) for value in row_content]
            virial_tensor.append(single_direction_vector)
                
        return np.array(virial_tensor)
        

    
    def get_atomic_numbers_lst(self):
        '''
        Description
        -----------
            1. 得到体系内所有的原子序数 (各个原子的 atomic_numbers)

        Note
        ----
            1. 重复
        '''
        ### Part I. 得到所有原子的原子序数、坐标
        ###         所有原子的原子序数：atomic_numbers_lst
        ###         所有原子的坐标：coordinations_lst
        content = "POSITION"    # 此处需要大写
        idx_row = LineLocator.locate_all_lines(
                                    file_path=self.atom_config_path,
                                    content=content)[0]
        with open(self.atom_config_path, 'r') as f:
            atom_config_content = f.readlines()
        
        # 1. 得到所有原子的原子序数（注意将读取的 str 转换为 int）
        atomic_numbers_content = atom_config_content[idx_row:idx_row + self.num_atoms]
        atomic_numbers_lst = [int(row.split()[0]) for row in atomic_numbers_content]

        return np.array(atomic_numbers_lst)


    def get_coords_lst(self):
        '''
        Description
        -----------
            1. 得到体系内所有的坐标
        '''
        coords_lst = []
        content = "POSITION"    # 此处需要大写
        idx_row = LineLocator.locate_all_lines(
                                    file_path=self.atom_config_path,
                                    content=content)[0]
        with open(self.atom_config_path, 'r') as f:
            atom_config_content = f.readlines()
        
        # 1. 得到所有原子的原子序数（注意将读取的 str 转换为 int）
        """
        row_content
        -----------
            '29         0.377262291145329         0.128590184800933         0.257759805813488     1  1  1'
        """
        for row_content in atom_config_content[idx_row:idx_row + self.num_atoms]:
            row_content_lst = row_content.split()
            coord_tmp = [float(value) for value in row_content_lst[1:4]]
            coords_lst.append(np.array(coord_tmp))
        
        return np.array(coords_lst)
    
    
    def get_atomic_forces_lst(self):
        '''
        Description
        -----------
            1. 得到体系内所有原子的受力
        '''
        try:    # atom.config 中有关于原子受力的信息
            ### Step 1. 得到 atom.config 文件中所有信息，以列表的形式组合
            forces_lst = []
            content = "Force".upper()
            idx_row = LineLocator.locate_all_lines(
                                    file_path=self.atom_config_path,
                                    content=content)[0]
            with open(self.atom_config_path, 'r') as f:
                atom_config_content = f.readlines()
            
            ### Step 2. 将力的信息，组织成
            for row_content in atom_config_content[idx_row:idx_row + self.num_atoms]:
                row_content_lst = row_content.split()
                force_tmp = [float(value) for value in row_content_lst[1:4]]
                forces_lst.append(np.array(force_tmp))
        
            return np.array(forces_lst)
        except: # atom.config 中没有关于原子受力的信息
            return np.zeros((self.num_atoms, 3))


    def get_atomic_velocitys_lst(self):
        '''
        Description
        -----------
            1. 得到体系内所有原子的速度
        '''
        try:    # atom.config 中有关于原子速度的信息
            ### Step 1. 得到 atom.config 文件中所有信息，以列表的形式组合
            velocitys_lst = []
            content = "Velocity (bohr/fs)".upper()
            idx_row = LineLocator.locate_all_lines(
                                    file_path=self.atom_config_path,
                                    content=content)[0]
            with open(self.atom_config_path, 'r') as f:
                atom_config_content = f.readlines()
            
            ### Step 2. 将速度的信息，组织成 np.ndarray 形式
            for row_content in atom_config_content[idx_row:idx_row + self.num_atoms]:
                row_content_lst = row_content.split()
                force_tmp = [float(value) for value in row_content_lst[1:4]]
                velocitys_lst.append(np.array(force_tmp))
        
            return np.array(velocitys_lst)
        except: # atom.config 中没有关于原子速度的信息
            return np.zeros((self.num_atoms, 3))


    def get_atomic_energys_lst(self):
        '''
        Description
        -----------
            1. 得到体系内所有原子的能量
        '''
        try:    # atom.config 中有关于原子能量的信息
            ### Step 1. 得到 atom.config 文件中所有信息，以列表的形式组合
            energys_lst = []
            content = "Atomic-Energy, ".upper()
            idx_row = LineLocator.locate_all_lines(
                                    file_path=self.atom_config_path,
                                    content=content)[0]
            with open(self.atom_config_path, 'r') as f:
                atom_config_content = f.readlines()
            
            ### Step 2. 将能量的信息，组织成 np.ndarray 形式
            for row_content in atom_config_content[idx_row:idx_row + self.num_atoms]:
                row_content_lst = row_content.split()
                energy_tmp = [float(value) for value in row_content_lst[1:4]]
                energys_lst.append(np.array(energy_tmp))
        
            return np.array(energys_lst)
        except: # atom.config 中没有关于原子能量的信息
            return np.zeros((self.num_atoms, 3))
    

    def get_magnetic_moments(self):
        '''
        Description
        -----------
            1. 得到所有原子的磁矩，顺序与 `坐标` 的顺序一致
        '''
        content = "MAGNETIC"

        magnetic_moments_lst = []
        
        try:    # 处理异常：若 atom.config 中不包含原子的磁矩信息
            idx_row = LineLocator.locate_all_lines(
                                    file_path=self.atom_config_path,
                                    content=content)[-1]

            with open(self.atom_config_path, "r") as f:
                atom_config_content = f.readlines()
            
            magnetic_moments_content = atom_config_content[idx_row: idx_row+self.num_atoms]
            # MAGNETIC  
            # 3 0.0 # 原子序数 磁矩
            # ...
            magnetic_moments_lst = [float(tmp_magnetic_moment.split()[-1]) for tmp_magnetic_moment in magnetic_moments_content]
        except Exception as e:
            #print(e)
            magnetic_moments_lst = [0 for _ in range(self.num_atoms)]
        
        return magnetic_moments_lst


    
    def get_atoms_lst(self):
        '''
        Description
        -----------
            1. 得到材料体系内所有 `ht_battery.cores.Atom` 对象
        '''
        atoms_lst = []
        
        ### Part I. 得到所有原子的原子序数、坐标
        ###         所有原子的原子序数：atomic_numbers_lst
        ###         所有原子的坐标：coordinations_lst
        content = "POSITION"
        idx_row = LineLocator.locate_all_lines(
                                    file_path=self.atom_config_path,
                                    content=content)[-1]
        with open(self.atom_config_path, 'r') as f:
            atom_config_content = f.readlines()
        
        # 1. 得到所有原子的原子序数（注意将读取的 str 转换为 int）
        atomic_numbers_content = atom_config_content[idx_row:idx_row + self.num_atoms]
        atomic_numbers_lst = [int(row.split()[0]) for row in atomic_numbers_content]

        # 2. 得到所有原子坐标的信息 (注意将读取的 str 转换为 int)
        coordinations_content = atom_config_content[idx_row:idx_row + self.num_atoms]
        coordinations_lst = [row.split()[1:4] for row in coordinations_content]
        for i in range(len(coordinations_lst)): # 转换
            for j in range(3):
                coordinations_lst[i][j] = float(coordinations_lst[i][j])

        ### Part II. 得到所有原子的磁矩 (注意将读取的 str 转换为 int)
        ###         所有原子的磁矩：magnetic_moments_lst
        try:
            content = "MAGNETIC"
            idx_row = LineLocator.locate_all_lines(
                                        file_path=self.atom_config_path,
                                        content=content)[-1]
            with open(self.atom_config_path, 'r') as f:
                atom_config_content = f.readlines()
            magnetic_moments_content = atom_config_content[idx_row: idx_row+self.num_atoms]
            magnetic_moments_lst = [float( row.split()[1] ) for row in magnetic_moments_content]
            #print("Information of magnetic moment exist.")
        except:
            #print("The program will extract magnetic moment from template atom.config...")
            magnetic_moments_lst = [0 for _ in range(self.num_atoms)]
        

        ### Part III
        for idx_atom in range(self.num_atoms):
            tmp_atom = Atom(
                        atomic_number=atomic_numbers_lst[idx_atom],
                        coordination=coordinations_lst[idx_atom],
                        magnetic_moment=magnetic_moments_lst[idx_atom]
                        )
            atoms_lst.append(tmp_atom)

        return atoms_lst
    
    
class AtomConfigStrExtractor(AtomConfigExtractorBase):
    '''
    Description
    -----------
        1. 从 str 中提取 atom.config 的信息
    '''
    def __init__(self, atom_config_str:str):
        self.atom_config_str = atom_config_str
        
    
    def get_num_atoms(self):
        pass
    
    
    def get_basis_vectors_lst(self):
        return super().get_basis_vectors_lst()
    
    
    def get_virial_tensor(self):
        return super().get_virial_tensor()
    
    
    def get_atomic_numbers_lst(self):
        return super().get_atomic_numbers_lst()
    
    
    def get_coords_lst(self):
        return super().get_coords_lst()
    
    
    def get_atomic_forces_lst(self):
        return super().get_atomic_forces_lst()
    
    
    def get_atomic_velocitys_lst(self):
        return super().get_atomic_velocitys_lst()
    
    
    def get_atomic_energys_lst(self):
        return super().get_atomic_energys_lst()
    
    
    def get_magnetic_moments(self):
        return super().get_magnetic_moments()
        