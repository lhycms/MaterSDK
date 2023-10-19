from abc import ABC, abstractmethod
from typing import List
import copy
import numpy as np

from .acextractor import ACstrExtractor
from .lineLocator import LineLocator

    
class MVExtractor(object):
    def __init__(
            self, 
            movement_path:str,
            virial:bool=False, 
            magmoms:bool=False, 
            eatoms:bool=False):
        self.mv_path = movement_path
        self.virial_mark = virial
        self.magmoms_mark = magmoms
        self.eatoms_mark = eatoms
        
        self.chunksizes:List[int] = self.get_chunksizes()
        self.chunkslices:List[int] = self.get_chunkslices()
    
    
    def get_chunksizes(self):
        '''
        Description
        -----------
            1. 由于MOVEMENT 文件太大，因此将 MOVEMENT 的每一帧 (frame) 对应的内容，定义为一个chunk
            文件处理时 `一个chunk一个chunk处理`
            
        Return
        ------
            1. chunk_size: int 
                - 每一帧的行数
                - e.g. [225, 382, 382, 382, 382, 382, ...]
        
        Note
        ----
            1. chunksize 包括 `-------` 这一行
        '''
        chunksizes_lst:List[int] = []
        content = "--------------------------------------"
        aim_row_idxs = LineLocator.locate_all_lines(
                        file_path=self.mv_path, 
                        content=content)

        chunksizes_lst.append(aim_row_idxs[0])
        for idx in range(1, len(aim_row_idxs)):
            chunksizes_lst.append(aim_row_idxs[idx] - aim_row_idxs[idx-1])
        return chunksizes_lst
    
    
    def get_chunkslices(self):
        '''
        Return
        ------
            1. chunkslices: np.ndarrary
                - e.g. [     0    225    607    989   1371, ... ]
        '''
        chunksizes = copy.deepcopy(self.chunksizes)
        chunksizes.insert(0, 0)
        chunksizes = np.cumsum(chunksizes)
        
        return chunksizes
    
    
    def get_frame_str(self, fidx:int):
        frame_str:str = ""
        with open(self.mv_path, "r") as mv:
            for idx_line, line in enumerate(mv):
                if idx_line in range(self.chunkslices[fidx], self.chunkslices[fidx+1]):
                    frame_str += line
                elif idx_line >= self.chunkslices[fidx+1]:
                    break
        return frame_str
        
    
    def get_frame_info(self, fidx:int):
        frame_str = self.get_frame_str(fidx=fidx)
        ace_str_extractor = ACstrExtractor(atom_config_str=frame_str)
        box = ace_str_extractor.get_basis_vectors()
        types = ace_str_extractor.get_types()
        coords = ace_str_extractor.get_coords()
        etot = ace_str_extractor.get_etot()
        fatoms = ace_str_extractor.get_fatoms()
        return box, types, coords, etot, fatoms
        
    
        
    def get_frames_info(self):
        pass
    

