import os 
from typing import List, Union, Optional

from .infoset import InfoSet
from ..io.pwmat.utils.mvextractor import MVExtractor
from ..io.pwmat.utils.parameters import specie2atomic_number, atomic_number2specie


class InfoSets(object):
    def __init__(self, dir_path:str, file_name:str, file_format:str):
        # Step 1. Search for all `file_name` in `dir_path`
        file_paths_lst:List[str] = []
        for root, dirs, files in os.walk(dir_path):
            for tmp_file in files:
                if tmp_file.upper() == file_name.upper():
                    file_paths_lst.append( os.path.join(root, tmp_file) )
        
        # Step 2. Get `infosets_atomic_numbers`
        infosets_atomic_numbers:List[int] = []
        infoset_lst:List[InfoSet] = []
        for tmp_file in file_paths_lst:
            infoset_lst.append(InfoSet(file_path=tmp_file, file_format=file_format))
        for tmp_infoset in infoset_lst:
            for tmp_element in tmp_infoset.formula_dict.keys():
                if not (specie2atomic_number[tmp_element] in infosets_atomic_numbers):
                    infosets_atomic_numbers.append(specie2atomic_number[tmp_element])
        
        # Step 3. Get `self.infoset_lst`
        self.infoset_lst:List[InfoSet] = []
        for tmp_file in file_paths_lst:
            self.infoset_lst.append(InfoSet(file_path=tmp_file, file_format=file_format, infosets_atomic_numbers=infosets_atomic_numbers))
            
        for tmp_idx1, tmp_infoset1 in enumerate( self.infoset_lst ):
            for tmp_idx2, tmp_infoset2 in enumerate( self.infoset_lst ):
                if (tmp_idx1 != tmp_idx2) and \
                        (tmp_infoset1.formula == tmp_infoset2.formula):
                    pass
                
        
        
    def get_num_frames(self):
        num_frames = 0
        for tmp_infoset in self.infoset_lst:
            num_frames += tmp_infoset.num_frames
        return num_frames
    
    
    @classmethod
    def from_dir(cls, dir_path:str, file_name:str, file_format:str):
        return cls(dir_path=dir_path, file_name=file_name, file_format=file_format)


    def to_dir(self, dir_path:str, part_size:Union[int, bool]=False):
        if part_size is False:
            part_size = self.get_num_frames()
        for tmp_infoset in self.infoset_lst:
            tmp_infoset.to_dir(dir_path=dir_path, part_size=part_size)
    