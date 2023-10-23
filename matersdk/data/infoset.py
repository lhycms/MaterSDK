from abc import ABC, abstractmethod
from typing import Dict, List, Union, Optional
import numpy as np
import os
import warnings

from ..io.pwmat.utils.mvextractor import MVExtractor
from ..io.pwmat.utils.parameters import atomic_number2specie


class InfoSet(object):
    def __init__(self, file_path:str, file_format:str):
        '''
        Description
        -----------
            1. 
        '''
        if (file_format.upper() == "PWMAT/MOVEMENT"):
            self.traj_extractor = MVExtractor(movement_path=file_path)
            self.virial_mark, self.eatoms_mark, self.magmoms_mark = self.traj_extractor._find_extra_properties();
        else:
            raise NameError("InfoSetError : No implementation for format of {0}".format(file_format))
    
        frames_info:List[np.ndarray] = self.traj_extractor.get_frames_info()
        self.box:np.ndarray = frames_info[0]
        self.types:np.ndarray = frames_info[1]
        self.coord:np.ndarray = frames_info[2]
        self.etot:np.ndarray = frames_info[3]
        self.fatoms:np.ndarray = frames_info[4]
        if self.virial_mark:
            self.virial:np.ndarray = frames_info[5]
        if self.eatoms_mark:
            if self.virial_mark:
                self.eatoms:np.ndarray = frames_info[6]
            else:
                self.eatoms:np.ndarray = frames_info[5]
        if self.magmoms_mark:
            if (self.virial_mark and self.eatoms_mark):
                self.magmoms:np.ndarray = frames_info[7]
            elif (self.virial_mark or self.eatoms_mark):
                self.magmoms:np.ndarray = frames_info[6]
            else:
                self.magmoms:np.ndarray = frames_info[5]
        
        self.num_frames:int = frames_info[0].shape[0]
        self.formula_dict:Dict[str, int] = {}
        for tmp_type, tmp_count in zip(np.unique(self.types[0], return_counts=True)[0], np.unique(self.types[0], return_counts=True)[1]):
            self.formula_dict.update({atomic_number2specie[tmp_type]: tmp_count})
        self.formula = ""
        for k, v in self.formula_dict.items():
            self.formula += "{0}{1}".format(k, v)
            

    @classmethod
    def from_file(cls, file_path:str, file_format:str):
        info_set = cls(file_path=file_path, file_format=file_format)
        return info_set

    
    def to_dir(self, dir_path:str, part_size:Union[int, bool]=False):
        if part_size is not False:
            num_parts = int(self.num_frames / part_size) + 1
        else:
            num_parts = 1
            part_size = self.num_frames
        
        if os.path.exists( os.path.join(dir_path, self.formula) ):
            warnings.warn("This dir exists: ".format(os.path.join(dir_path, self.formula)))
        os.mkdir(os.path.join(dir_path, self.formula))
        
        for part_rank in range(num_parts):
            start:int = part_rank * part_size
            if part_rank == (num_parts-1):
                end:int = self.num_frames
            else:
                end:int = start + part_size
            if (end == start):
                break
                        
            ### Step 1. 
            tmp_part_dir_path:str = os.path.join(dir_path, self.formula, "part.{0:0>3}".format(part_rank))
            if os.path.exists(tmp_part_dir_path):
                warnings.warn("This dir contains part.{0:0>3}".format(part_rank))
            os.mkdir(tmp_part_dir_path)
            np.save(file=os.path.join(tmp_part_dir_path, "box.npy"), arr=self.box[start:end, :])
            np.save(file=os.path.join(tmp_part_dir_path, "types.npy"), arr=self.types[start:end, :])
            np.save(file=os.path.join(tmp_part_dir_path, "coord.npy"), arr=self.coord[start:end, :])
            np.save(file=os.path.join(tmp_part_dir_path, "etot.npy"), arr=self.etot[start:end, :])
            np.save(file=os.path.join(tmp_part_dir_path, "fatom.npy"), arr=self.fatoms[start:end, :])
            if self.virial_mark:
                np.save(file=os.path.join(tmp_part_dir_path, "virial.npy"), arr=self.virial[start:end, :])
            if self.eatoms_mark:
                np.save(file=os.path.join(tmp_part_dir_path, "eatom.npy"), arr=self.eatoms[start:end, :])
            if self.magmoms_mark:
                np.save(file=os.path.join(tmp_part_dir_path, "magmoms"), arr=self.magmoms[start:end, :])
        