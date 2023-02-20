import linecache
import numpy as np


class Report(object):
    COLUMN_PER_LINE = 5 # REPORT中每行出现几个 eigen energies
    
    def __init__(self, report_path:str):
        self.report_path = report_path
    
    
    def _search_aim(self, aim_content:str):
        '''
        Description
        -----------
            1. 查询REPORT文件中是否存在特定的内容(aim_content)，并确定所在的行数
        
        Parameters
        ----------
            1. aim_content: str
                - 特定的内容
        
        Return
        ------
            1. idxs_lst: List[int]
        '''
        with open(self.report_path, "r") as f:
            lines_lst = f.readlines()
        
        idxs_lst = []
        for tmp_idx, tmp_line in enumerate(lines_lst, 1):
            if aim_content in tmp_line:
                idxs_lst.append(tmp_idx)
        
        return idxs_lst
    
    
    def get_num_bands(self):
        '''
        Description
        -----------
            1. 得到能带数 (每个kpoint的本征态)
        
        Return
        ------
            1. num_bands: int
                - 能带的数目
        '''
        ### Step 1. 查询得到 `NUM_BAND` 所在的行
        aim_content = "NUM_BAND"
        idx_num_bands = self._search_aim(aim_content=aim_content)[0]
        #print(idx_num_bands)
        ### Step 2. 提取能带的数目
        specific_line = linecache.getline(self.report_path, idx_num_bands)
        num_bands = int( specific_line.split()[-1] )
        
        return num_bands
    
    
    def get_num_kpts(self):
        '''
        Description
        -----------
            1. 得到Kpoint的数目
        
        Return
        ------
            1. num_kpts: int
                - kpoints 的数目
        '''    
        ### Step 1. 查询得到 `NUM_KPT` 所在的行
        aim_content = "NUM_KPT"
        idx_num_kpts = self._search_aim(aim_content=aim_content)[0]
        
        ### Step 2. 提取kpoints的数目
        specific_line = linecache.getline(self.report_path, idx_num_kpts)
        num_kpts = int( specific_line.split()[-1] )
        
        return num_kpts
        
    
    def get_eigen_energies(self):
        '''
        Description
        -----------
            1. 
        
        Return
        ------
            1. spin2eigen_energies: Dict[str, np.ndarray]
                - np.ndarray 一维: kpoint
                - np.ndarray 二维: 某kpoint的本征能量
        '''
        ### Step 1. 
        ###     1. 初始化 `spin2eigen_energies`
        ###     2. 得到kpoints的数目 `num_kpts`
        ###     3. 得到bands的数目 `num_bands`
        ###     4. 得到 `idx_eigen_start_lst`
        spin2eigen_energies = {"up":[], "down":[]}
        num_kpts = self.get_num_kpts()
        num_bands = self.get_num_bands()
        aim_content_eigen = "eigen energies, in eV"
        idxs_eigen_start_lst = self._search_aim(
                                aim_content=aim_content_eigen)
        
        ### Step 2. 读取 REPORT 文件
        with open(self.report_path, "r") as f:
            lines_lst = f.readlines()
        
        ### Step 3. 得到每个kpoint的本征能量
        num_lines_for_band = int( np.ceil(num_bands / self.COLUMN_PER_LINE) )
        for tmp_idx, tmp_idx_eigen_start in enumerate(idxs_eigen_start_lst):
            tmp_eigen_energies_ = lines_lst[tmp_idx_eigen_start : tmp_idx_eigen_start+num_lines_for_band]
            tmp_eigen_energies = [float(eigen) for tmp_5_eigen in tmp_eigen_energies_ for eigen in tmp_5_eigen.split()]
            tmp_eigen_energies_array = np.array( tmp_eigen_energies )
            
            if tmp_idx < num_kpts:
                spin2eigen_energies["up"].append(tmp_eigen_energies_array)
            else:
                spin2eigen_energies["down"].append(tmp_eigen_energies_array)
        
        ### Step 4. 将 spin2igen_energies 的 values 变为 np.ndarray 形式
        spin2eigen_energies.update(
                        {"up": np.array( spin2eigen_energies["up"] )}
                        )
        spin2eigen_energies.update(
                        {"down": np.array( spin2eigen_energies["down"] )}
                        )
        
        ### Step 5. 当 ispin 打开时，自旋向上和向下的(kpoints, eigen_states)应该相等
        if spin2eigen_energies["down"].size != 0:
            assert (spin2eigen_energies["up"].shape != spin2eigen_energies["down"])
        
        
        return spin2eigen_energies