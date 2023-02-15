import linecache
import numpy as np


from ...publicLayer.structure import DStructure


class Inkpt(object):
    '''
    Description
    -----------
        1. This class is aimed to `IN.KPT`
    
    Note
    ----
        1. 注意单位换算：埃 <-> Bohr
            - IN.KPT 中用的是 `埃`
            - REPORT 中用的是 `Bohr`
        2. 目前仅支持 `iflag, a0 = 2, 0`
    '''
    def __init__(self, in_kpt_path:str):
        '''
        Parameters
        ----------
            1. in_kpt_path: str
                - IN.KPT 文件的路径
        '''
        self.in_kpt_path = in_kpt_path
        self.iflag = self._get_iflag()
        self.a0 = self._get_a0()
    
    
    def _get_iflag(self):
        '''
        Description
        -----------
            1. 得到 IN.KPT 中的 `iflag`
                - 用于判断 `K 点位置定义在 x, y, z 方向` 或 `K 点位置定义在 AL(3,3) 的倒格子中`
        
        Return
        ------
            1. a0: int
                - 
        '''
        line_2 = linecache.getline(self.in_kpt_path, 2)
        iflag = int( line_2.split()[0] )
        return iflag
    
    
    def _get_a0(self):
        '''
        Description
        -----------
            1. 得到 IN.KPT 中的 `a0`
                - 仅在 iflag = 1 时使用 (原子单位 Bohr)
        
        Return
        ------
            1. a0: int
                - 
        '''
        line_2 = linecache.getline(self.in_kpt_path, 2)
        a0 = int( line_2.split()[-1] )
        return a0
    
    
    def get_num_kpts(self):
        '''
        Description
        -----------
            1. 得到 kpoints 的个数
        
        Return
        ------
            1. num_kpts: int
        '''
        line_1 = linecache.getline(self.in_kpt_path, 1)
        num_kpts = int( line_1.split()[0] )
        return num_kpts
    
    
    def get_kpt_coords_frac(self):
        '''
        Description
        -----------
            1. 得到 IN.KPT 中所有 KPOINTS 的分数坐标
            
        Return
        ------
            1. kpt_coord_frac: np.ndarray
                - 所有 kpoints 的坐标
        '''
        ### Step 1. 读取 `仅包含坐标和权重的行`
        with open(self.in_kpt_path, "r") as f:
            lines_lst = f.readlines()
        lines_lst = lines_lst[2:]
        
        ### Step 2. 取出所有的 Kpoints 坐标
        coords_frac_lst = []
        for tmp_line in lines_lst:
            tmp_coord_frac = [float(tmp_value) \
                        for tmp_value in tmp_line.split()[:3]]
            coords_frac_lst.append(tmp_coord_frac)
        coords_frac_array = np.array(coords_frac_lst)
        
        return coords_frac_array
    

    def get_kpt_coords_A(
                    self, 
                    atom_config_path:str):
        '''
        Description
        -----------
            1. 得到所有 kpoints 的坐标 (单位：埃)，与 IN.KPT/OUT.KPT 对应
        '''
        ### Step 1. 得到倒易晶格（unit: 埃）
        structure = DStructure.from_file(
                            file_path=atom_config_path,
                            file_format="pwmat",
                            )
        reciprocal_lattice_array = np.array( structure.lattice.reciprocal_lattice.matrix )

        ### Step 2. 得到kpoints的分数坐标
        kpt_coords_frac = self.get_kpt_coords_frac()
        
        ### Step 3. kpoints的分数坐标 * 倒易晶格
        return np.dot(kpt_coords_frac, reciprocal_lattice_array)        
    
    
    def get_kpt_coords_Bohr(
                    self, 
                    atom_config_path:str):
        '''
        Description
        -----------
            1. 得到所有 kpoints 的坐标 (单位：Bohr)，与 REPORT 对应
        '''
        ### Step 1. 得到倒易晶格（unit: 埃）
        structure = DStructure.from_file(
                            file_path=atom_config_path,
                            file_format="pwmat",
                            )
        reciprocal_lattice_array = np.array( structure.lattice.reciprocal_lattice.matrix )

        ### Step 2. 得到kpoints的分数坐标
        kpt_coords_frac = self.get_kpt_coords_frac()
        
        ### Step 3. kpoints的分数坐标 * 倒易晶格
        return np.dot(kpt_coords_frac, reciprocal_lattice_array) * 0.529177249
    
    
    def get_kpt_weights(self):
        '''
        Description
        -----------
            1. 得到 kpoints 的权重
            
        Return
        ------
            1. weights_lst: List[float]
                - kpoints 的权重
        '''
        pass
    
    
    def get_hsp(self):
        '''
        Description
        -----------
            1. 得到 IN.KPT 中的高对称点和坐标
        
        Return
        ------
            1. hsp2coord_frac: Dict[str, list]
        '''
        pass
    
    
    def get_distance_from_gamma_A(self):
        '''
        Description
        -----------
            1. 
        '''
        pass
    
    
    def get_distance_from_gamma_Bohr(self):
        '''
        Description
        -----------
            1. 
        '''
        pass