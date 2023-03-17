import pandas as pd
import numpy as np
from io import StringIO
from ..utils.lineLocator import LineLocator


class FatbandStructure(object):
    '''
    Description
    -----------
        1. 使用 `plot_fatbandstructure.x` 之后，会产生 `fatbandstruture_1.txt` 文件，
           此 class 用于解析这个文件的各类信息
    
    Parameters
    ----------
        1. fatbandstructure_txt_path: str
            - `fatbandstructure_1.txt` 文件的路径
        2. KPOINT 的距离单位转换为`埃`
    '''
    BOHR = 0.529177249
    
    def __init__(
                self,
                fatbandstructure_txt_path:str
                ):
        self.fatbandstructure_txt_path = fatbandstructure_txt_path
        self.num_bands = self._get_num_bands()
        self.num_kpoints = self._get_num_kpoints()
    
    
    def _get_num_bands(self):
        '''
        Description
        -----------
            1. 得到能带的数目
        '''
        idx_lines_lst = LineLocator.locate_all_lines(
                    file_path=self.fatbandstructure_txt_path,
                    content="BAND"
                    )
        return len(idx_lines_lst)

    
    def _get_num_kpoints(self):
        '''
        Description
        -----------
            1. 得到 kpoints 的数目
        '''
        idx_lines_lst = LineLocator.locate_all_lines(
                    file_path=self.fatbandstructure_txt_path,
                    content="BAND",
        )
        num_kpoints = idx_lines_lst[1] - idx_lines_lst[0] - 2
        return num_kpoints
    
    
    def _get_BAND_mark_idxs(self):
        '''
        Description
        -----------
            1. 找到 `BAND` 所在的行，并返回所在行的索引
        
        Note
        ----
            1. 索引是从 1 开始的，便于 `linecache.getline()` 调用
        '''
        idx_lines_lst = LineLocator.locate_all_lines(
            file_path=self.fatbandstructure_txt_path,
            content="BAND"
            )
        return idx_lines_lst


    def _preprocess(self):
        '''
        Description
        -----------
            1. 预处理，将 `fatbandstructure_1.txt` 读取成 pd.DataFrame，
            2. 例如有64条能带、29个kpoint的能带结构，会产生 64*29=1856 columns
        
        Note
        ----
            1. 跳过 `BAND 行` 和 `空行`
        '''
        ### Step 1. 跳过 `BAND 行` 和 `空行`
        ss = StringIO()
        with open(self.fatbandstructure_txt_path, 'r') as f:
            for line in f:
                if (line=='' or "BAND" in line):
                    continue
                else:
                    ss.write(line)
        ss.seek(0)   # "rewind" to the beginning of the StringIO object
        
        df = pd.read_csv(
                        ss,
                        delimiter='\s+'
        )
        df.loc[:, "KPOINT"] = df.loc[:, "KPOINT"] / self.BOHR
        ### Step 2. 每条能带组成新的 DataFrame
        ### type(dfs) = List[pd.DataFrame]
        dfs_lst = np.array_split(df, self.num_bands)
        assert ( len(dfs_lst) == self.num_bands )
        return dfs_lst
    
    
    def get_dfs(self):
        '''
        Description
        -----------
            1. 
            
        Return
        ------
            1. dfs_lst: List[pd.DataFrame]
                - 每个 DataFrame 代表一条能带
        '''
        dfs_lst = self._preprocess()
        return dfs_lst