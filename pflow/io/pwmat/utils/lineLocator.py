'''
Author       : Liu Hanyu
Email        : hyliu2016@buaa.edu.cn
Date         : 2022-10-31 15:20:26
LastEditTime : 2022-11-03 13:01:27
FilePath     : /pflow/pflow/io/utilitys/lineLocator.py
Description  : 
'''
class LineLocator(object):
    @staticmethod
    def locate_all_lines(file_path:str, content:str):
        '''
        Description
        -----------
            1. 定位某段文本所在的行 (返回所有行数)

        Parameters
        ----------
            1. file_path: str
                文件的绝对路径
            2. content: str
                需要定位的内容
        '''
        row_idxs_lst = []
        row_no = 0

        with open(file_path, "r") as f:
            for row_content in f:
                row_no += 1
                
                if content in row_content.upper():
                    row_idxs_lst.append(row_no)
        
        return row_idxs_lst