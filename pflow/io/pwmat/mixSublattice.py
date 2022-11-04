'''
Author       : Liu Hanyu
Email        : hyliu2016@buaa.edu.cn
Date         : 2022-09-05 11:44:08
LastEditTime : 2022-11-03 15:09:59
FilePath     : /pflow/pflow/io/pwmat/mixSublattice.py
Description  : e.g. 将 Li子晶格(挖了空位) 和 MnPO4子晶格(做了掺杂) 混合
'''
from ..publicLayer.structure import DStructure



class MixSublattice(object):
    '''
    '''
    def __init__(self,
                atom_config,
                atom_config_sublattice,
                atom_config_output
                ):
        '''
        Description
        -----------
            1. 

        Steps
        -----
            1. 在 atom_config_sublattice 中挑选对应元素的行
            2. 在 atom_config 中删除对应元素的行
            3. 将第一步挑选的行，添加到 atom_config 文件中

        Parameters
        ----------
            1. atom_config: str
                以这个 atom_config 为主体，进行操作 (包含 MnPO4 子晶格的文件)
            2. atom_config_sublattice: str  
                从 atom_config_sublattice 中提取 sublattice 信息 (包含 Li 子晶格的文件)
        '''
        # 1. 包含 MnPO4 子晶格的文件
        self.atom_config = atom_config
        # 2. 包含 Li/Na 子晶格的文件
        self.atom_config_sublattice = atom_config_sublattice
        self.atom_config_output = atom_config_output
    

    def _pick_Atoms_accordingto_atomic_number(self,
                                    atomic_number:int):
        '''
        Step 1. 

        Description
        -----------
            1. 将 `atom_config_sublattice` 中的 element (Li) 元素的行提取出来
        
        Parameters
        ----------
            1. atomic_number: int
                需要挑选的原子

        Note
        ----
            1. 被 `write_atom_config` 调用
        '''
        picked_Atoms_row_lst = []

        with open(self.atom_config_sublattice, 'r') as f:
            contents = f.readlines()

        # 得到 POSITION 部分的行数
        num_atoms = Structure.from_file(self.atom_config_sublattice).num_atoms

        # 遍历 POSITION 部分，防止包含 MAGNETIC 时报错
        for row_content in contents[6: num_atoms+6]:
            atomic_number_row = row_content.split()[0]
            if ( int(atomic_number_row) == atomic_number) :
                picked_Atoms_row_lst.append(row_content)
        
        return picked_Atoms_row_lst

    
    def _delete_Atoms_accordingto_atomic_number(self,
                                    atomic_number:int):
        '''
        Step 2. 处理提供 MnPO4 子晶格的 atom.config


        Note
        ----
            1. 被 `self.write_atom_config()` 调用
        '''
        picked_Atoms_row_lst = []

        with open(self.atom_config, 'r') as f:
            contents = f.readlines()

        # Part I. 前6行，信息保持一致
        for row in contents[:6]:
            picked_Atoms_row_lst.append(row)
        
        # Part II. 删除指定元素对应的行 -- 遍历 POSITION 部分
        # 得到 POSITION 部分的行数，防止包含 MAGNETIC 时报错
        num_atoms = Structure.from_file(self.atom_config).num_atoms
        for row_content in contents[6: num_atoms+6]:
            atomic_number_row = row_content.split()[0]
            if ( int(atomic_number_row) != atomic_number ):
                picked_Atoms_row_lst.append(row_content)

        return picked_Atoms_row_lst

    
    def write_atom_config(self,
                        atomic_number):
        '''
        Step 3. 
        '''
        picked_Atoms_row_lst_1 = self._delete_Atoms_accordingto_atomic_number(
                                                        atomic_number=atomic_number)
        picked_Atoms_row_lst_2 = self._pick_Atoms_accordingto_atomic_number(
                                                        atomic_number=atomic_number)

        picked_Atoms_row_lst_1 += picked_Atoms_row_lst_2
        # print(picked_Atoms_row_lst_1)
        with open(self.atom_config_output, 'w') as f:
            f.writelines( ["{0}\n".format(len(picked_Atoms_row_lst_1) - 6)] )
            f.writelines(picked_Atoms_row_lst_1[1:-1])
            if picked_Atoms_row_lst_1[-1][:-1] == '\n':
                f.writelines([ picked_Atoms_row_lst_1[-1][:-1] ])
            else:
                f.writelines( [picked_Atoms_row_lst_1[-1][:]] )