import linecache
import numpy as np
import os 

from ..utils.lineLocator import LineLocator
from ...publicLayer.atom import Atom
from ..utils.parameters import atomic_number2specie
from ...publicLayer.structure import DStructure


class MovementExtractor(object):
    def __init__(self,
                movement_path:str
                ):
        '''
        Parameters
        ----------
            1. movement_path: str
                - MOVEMENT 文件的路径
            2. 
        '''
        self.movement_path = movement_path
        self.chunksize = self.get_chunksize()
        
        
    
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
        content = "-------------------------------------------------"
        row_idxs = LineLocator.locate_all_lines(
                                    file_path=self.movement_path,
                                    content=content
                                    )
        chunksize = row_idxs[0]
        return chunksize
    
    
    def _get_frame_str(self, idx_frame:int):
        '''
        Description
        -----------
            1. 得到某一帧的str
            2. Note: 帧数从 1 开始
        
        Parameter
        ---------
            1. idx_frame: int
                - 得到代表某一帧的 str
            
        Return
        ------
            1. str_frame: str
                - 某一帧的 str
        '''
        assert (idx_frame > 0)
        
        
        ### Step 1. 得到对应的行数并组成str形式
        tmp_idx_frame = 0
        with open(self.movement_path, 'r') as f_mvt:
            while True: # 循环读取文件，直至文件末尾
                ### Step 1.1. 帧数从 1 开始计数；每次开始新的一帧都要用空的str
                tmp_idx_frame += 1
                tmp_chunk = ""
                
                ### Step 1.2. 收集这个 frame 对应的行数，组成 tmp_chunk
                for tmp_row_idx in range(self.chunksize):
                    tmp_row = f_mvt.readline()  # Read the next row from the file
                    if not tmp_row: 
                        break   # End of file reached
                    tmp_chunk += tmp_row
                    
                ### Step 1.3. 当处理到对应的帧数时，返回对应的 tmp_chunk
                if (tmp_idx_frame == idx_frame):
                    break
                        
        return tmp_chunk
        
    
    def get_frame_structure(self, idx_frame:int):
        '''
        Description
        -----------
            1. 将某一帧的结构抽取出来，构建成 DStrucure 对象
        
        Parameter
        ---------
            1. idx_frame: int
                - 第几帧（从 1 开始计数）
        
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
        #os.remove(self.tmp_struct_file)
        pass