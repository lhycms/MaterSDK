import pandas as pd
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
    '''
    def __init__(
                self,
                fatbandstructure_txt_path:str
                ):
        self.fatbandstructure_txt_path = fatbandstructure_txt_path
        self.num_bands = self._get_num_bands()
    
    
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
    
    
    def get_weights_orbitals(self):
        '''
        Description
        -----------
            1. 得到
        '''
        df_data = pd.read_csv(
                    self.fatbandstructure_txt_path,
                    delimiter='\s+',
        )
        
        print(df_data)