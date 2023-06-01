

class FFInfer(object):
    def __init__(self, ff_path:str):
        '''
        Parameters
        ----------
            1. ff_path: str
                - The path of ff file.
                
        '''
        pass
    
    
    def _get_ff_params(self, ff_path:str):
        '''
        Description
        -----------
            1. 从`.ff`文件中读取力场的参数
        
        Parameters
        ----------
            1. ff_path: str
                - `.ff` 文件的路径
        '''
        pass
    
    
    def _featurize(self,
                struct_file_path:str,
                struct_file_fmt:str):
        '''
        Description
        -----------
            1. 根据结构文件提取其 Deepmd 的feature
            
        Parameters
        ----------
            1. struct_file_path: str
                - 结构文件的路径
            2. struct_file_fmt: str
                - 结构文件的格式
        '''
        pass
    
    
    def _norm(self, davg, dstd):
        '''
        Description
        -----------
            1. 

        Parameters
        ----------
            1. davg: 
            2. dstd: 
        '''
        pass