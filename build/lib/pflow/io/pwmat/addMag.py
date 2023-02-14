import numpy as np

from ..publicLayer.structure import DStructure
from ..publicLayer.atom import Atom


class AtomConfigMagTem(object):
    '''
    
    '''
    def __init__(self,
                atom_config_path_template:str,
                atom_config_path_sqs:str):
        '''
        Parameters
        ----------
            1. atom_config_path_template: str
                模板 atom.config 的路径（包含各个原子的磁矩）
            2. atom_config_path_sqs: str
                sqs 产生的 atom.config 的路径 (需要填写各个原子的磁矩)
        '''
        self.atom_config_path_template = atom_config_path_template
        self.atom_config_path_sqs = atom_config_path_sqs
        self.structure_template = DStructure.from_file(
                                        file_path=atom_config_path_template,
                                        file_format="pwmat",
                                        coords_are_cartesian=False
                                        )
        self.structure_new = DStructure.from_file(
                                        file_path=atom_config_path_sqs,
                                        file_format="pwmat",
                                        coords_are_cartesian=False
                                        )
    

    def assign_magnetic_moment(self, atomic_number_zero_lst:list=None):
        '''
        Description
        -----------
            1. 根据原子最近原则，根据 atom.config 模板的磁矩给新 atom.config 赋值

        Parameters
        ----------
            1. atomic_number_zero_lst: list
                哪些元素的磁矩需要设置为 0
        
        Return
        ------
            1. atoms_add_magnetic_moment: list of `Atom`
                按照模板填写磁矩
            2. atoms_add_magnetic_moment_zero: list of `Atom`
                将没有磁矩的原子的磁矩归零
        '''
        ### Part I.根据原子最近原则，根据 atom.config 模板的磁矩给新 atom.config 赋值
        atoms_add_magnetic_moment = []
        for atom_1 in self.structure_new.atoms_lst:
            distance2_min = 1000
            for atom_2 in self.structure_template.atoms_lst:
                distance2 = AtomConfigMagTem.calculate_distance2_bet_atoms(
                                                                atom_1=atom_1,
                                                                atom_2=atom_2)
                #print(atom_1.atomic_number, atom_2.atomic_number, distance2)
                if (distance2 < distance2_min):
                    distance2_min = distance2
                    atom_1.magnetic_moment = atom_2.magnetic_moment
            #print('*' * 60)
            atoms_add_magnetic_moment.append(atom_1)
        
        self.structure_new.Atoms_lst = atoms_add_magnetic_moment

        ### Part II. 设置某些元素的磁矩为零
        atoms_add_magnetic_moment_zero = []
        for tmp_Atom in self.structure_new.Atoms_lst:
            if tmp_Atom.atomic_number in atomic_number_zero_lst:
                tmp_Atom.magnetic_moment = 0
            atoms_add_magnetic_moment_zero.append(tmp_Atom)
        self.structure_new.Atoms_lst = atoms_add_magnetic_moment_zero
        
        return self.structure_new.Atoms_lst
    

    def append_mag_to_atomconfig(self):
        '''
        Description
        -----------
            1. 把用模板赋值后的 `ht_battery.cores.structure.Structure` 对象的原子磁矩添加到 `atom.config` 文件中 

        Note
        ---- 
            1. 在调用本函数前，请先调用 `self.assign_magnetic_moment()`
            2. 本函数调用了 `AtomConfigMagTem.remove_atom_config_mag()`
            3. 若 `self.atom_config_path_template` 文件本来就有磁矩信息，本函数会自动将其删除后再按照模板添加磁矩信息
        '''
        AtomConfigMagTem.remove_atom_config_mag(self.atom_config_path_sqs)
        with open(self.atom_config_path_sqs, "a") as f:
            f.write("\nMAGNETIC\n")
            for tmp_atom in self.structure_new.Atoms_lst:
                f.write("{0} {1}\n".format(tmp_atom.atomic_number, tmp_atom.magnetic_moment))
        
    
    @staticmethod
    def remove_atom_config_mag(atom_config_path: str):
        '''
        Description
        -----------
            1. 清除 `atom.config` 文件中的磁矩信息
        
        Parameters
        ----------
            1. atom_config_path: str
                需要清理磁矩信息的 `atom.config` 的路径
        '''
        structure = DStructure.from_file(file_path=atom_config_path,
                                        file_format="pwmat",
                                        coords_are_cartesian=False,
                                        )
        # atom.config 前6行为无效信息
        count_row = structure.num_sites + 6
        with open(atom_config_path, "r") as f:
            lines_content = f.readlines()[:count_row]
        #print(lines_content)
        with open(atom_config_path, "w") as f:
            f.writelines(lines_content[:-1])
            if lines_content[-1][:-1] == '\n':
                f.writelines([ lines_content[-1][:-1] ])
            else:
                f.writelines( [lines_content[-1][:]] )


    @staticmethod
    def calculate_distance2_bet_atoms(atom_1: Atom, atom_2: Atom):
        '''
        Description
        -----------
            1. 计算两原子间的距离
        '''
        coordination_1 = np.array( atom_1.coordination )
        coordination_2 = np.array( atom_2.coordination )
        diff_coordination = coordination_2 - coordination_1
        distance2 = np.sum( np.power(diff_coordination, 2) , axis=0)
        return distance2