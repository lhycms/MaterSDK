import os 
import re
import numpy as np

from ..utils.lineLocator import LineLocator
from ...publicLayer.atom import Atom
from ..utils.parameters import atomic_number2specie
from ...publicLayer.structure import DStructure


class Movement(object):
    def __init__(self,
                movement_path:str
                ):
        '''
        Parameters
        ----------
            1. movement_path: str
                - MOVEMENT 文件的路径
            2. 
            
        Note
        ----
            1. MOVEMENT 第一步与其他步的 chunksize 不同。
        '''
        self.movement_path = movement_path
        self.chunksizes_lst = self.get_chunksize()
        
        
    def get_chunksize(self,
                    ):
        '''
        Description
        -----------
            1. 由于MOVEMENT 文件太大，因此将 MOVEMENT 的每一帧 (frame) 对应的内容，定义为一个chunk
            文件处理时 `一个chunk一个chunk处理`
            
        Return
        ------
            1. chunk_size: int 
                - 每一帧的行数
        
        Note
        ----
            1. chunksize 包括 `-------` 这一行
        '''
        chunksizes_lst = []
        content = "--------------------------------------"
        row_idxs = LineLocator.locate_all_lines(
                                    file_path=self.movement_path,
                                    content=content
                                    )

        chunksizes_lst.append(row_idxs[0])
        for idx in range(1, len(row_idxs)):
            chunksizes_lst.append(row_idxs[idx] - row_idxs[idx-1])
        
        return chunksizes_lst
    
    
    def _get_frame_str(self, idx_frame:int):
        '''
        Description
        -----------
            1. 得到某一帧的str
            2. Note: 帧数从 0 开始
        
        Parameter
        ---------
            1. idx_frame: int
                - 得到代表某一帧的 str
            
        Return
        ------
            1. str_frame: str
                - 某一帧的 str
        '''        
        ### Step 1. 得到对应的行数并组成str形式
        tmp_idx_frame = 0
        with open(self.movement_path, 'r') as f_mvt:
            while True: # 循环读取文件，直至文件末尾
                ### Step 1.1. 帧数从 0 开始计数；每次开始新的一帧都要用空的str
                tmp_chunk = ""
                tmp_chunksize = self.chunksizes_lst[tmp_idx_frame]
                
                ### Step 1.2. 收集这个 frame 对应的行数，组成 tmp_chunk
                for tmp_row_idx in range(tmp_chunksize):
                    tmp_row = f_mvt.readline()  # Read the next row from the file
                    if not tmp_row: 
                        break   # End of file reached
                    tmp_chunk += tmp_row
                    
                ### Step 1.3. 当处理到对应的帧数时，返回对应的 tmp_chunk
                if (tmp_idx_frame == idx_frame):
                    break
                
                ### Step 1.4. 将读取下一帧的 chunk
                tmp_idx_frame += 1 
                        
        return tmp_chunk
        
    
    def get_frame_structure(self, idx_frame:int):
        '''
        Description
        -----------
            1. 将某一帧的结构抽取出来，构建成 DStrucure 对象
        
        Parameters
        ----------
            1. idx_frame: int
                - 第几帧（从 0 开始计数）
        
        Return
        ------
            1. structure: DStructure
                - 
        '''
        str_frame = self._get_frame_str(idx_frame=idx_frame)
        structure = None
        
        with CreateAndRemove() as context:
            with open(context.tmp_struct_file, 'w') as f:
                f.write(str_frame)
                f.seek(0)   # restart
                structure = DStructure.from_file(
                            file_path=context.tmp_struct_file,
                            file_format="pwmat",
                            coords_are_cartesian=False
                            )
        return structure
    
    
    def get_frame_energy(self, idx_frame:int):
        '''
        Description
        -----------
            1. 获取某一帧的总能、势能、动能
                72 atoms,Iteration (fs) =   -0.1000000000E+01, Etot,Ep,Ek (eV) =   -0.1188642969E+05  -0.1188642969E+05   0.0000000000E+00, SCF =    16
                Lattice vector (Angstrom)
                0.8759519000E+01    0.0000000000E+00    0.0000000000E+00     stress (eV):  0.124196E+02  0.479262E+01  0.245741E+01
                0.2209000000E+00    0.7513335000E+01    0.0000000000E+00     stress (eV):  0.479308E+01  0.961132E+01  0.225365E+01
                0.4093050000E+00    0.2651660000E+00    0.1828974400E+02     stress (eV):  0.245631E+01  0.225430E+01 -0.198978E+01
        
        Paramters
        ---------
            1. idx_frame: int
                - 第几帧 (从第 0 帧开始计数)
        
        Return
        ------
            1. energy_tot: float
            2. energy_p: float
            3. energy_k: float
        '''
        ### 1. 获取 `idx_frame` 对应的 chunk
        frame_str = self._get_frame_str(idx_frame=idx_frame)
        
        ### 2. 获取
        pattern = re.compile(r"Etot,Ep,Ek (eV) =.*?,")
        content = pattern.search(frame_str)
        return content



class CreateAndRemove(object):
    '''
    Descripion
    ----------
        1. 上下文背景
        2. Enter: 创建文件
        3. Exit: 删除文件
    '''
    def __init__(self):
        current_path = os.getcwd()
        self.tmp_struct_file = os.path.join(current_path, "tmp_structure_file")
    
    
    def __enter__(self):
        # 创建 `tmp_structure_file` 文件
        if os.path.isfile(self.tmp_struct_file):
            os.remove(self.tmp_struct_file)
        os.mknod(self.tmp_struct_file)
        
        return self
        
        
    def __exit__(self, exc_type, exc_value, traceback):
        os.remove(self.tmp_struct_file)