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
        
        Note
        ----
            1. content 必须为大写
        '''
        row_idxs_lst = []
        row_no = 0

        with open(file_path, "r") as f:
            for row_content in f:
                row_no += 1
                
                if content in row_content.upper():
                    row_idxs_lst.append(row_no)
        
        return row_idxs_lst