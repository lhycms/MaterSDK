import os 
import numpy as np

class Dostotalspin(object):
    def __init__(self,
                dos_totalspin_path:str):
        self.dos_totalspin_path = dos_totalspin_path
    

    def get_tdos(self):
        '''
        Description
        ----------- 
            1. 读取 DOS.totalspin, DOS.spinup, DOS.spindown
        
        Return
        ------
            1. tdos_array: np.ndarray
                    [[-65.075    0.       0.     ...   0.       0.       0.    ]
                    [-65.057    0.       0.     ...   0.       0.       0.    ]
                    [-65.039    0.       0.     ...   0.       0.       0.    ]
                    ...
                    [  6.0108   0.       0.     ...   0.       0.       0.    ]
                    [  6.0285   0.       0.     ...   0.       0.       0.    ]
                    [  6.0463   0.       0.     ...   0.       0.       0.    ]]
            
        '''
        ### Step 1. 删除 DOS.totalspin, DOS.spinup, DOS.spindown 行首的 "#"
        with open(self.dos_totalspin_path, "r") as f:
            lines_lst = f.readlines()
            # first_row: "#  Energy  Total  Mo-s  Mo-p  Mo-s  Mo-d  S-s  S-p"
            first_line = lines_lst[0]
            first_line_lst = first_line.split()
            try:
                first_line_lst.remove("#")
            except:
                pass
            new_first_line = "\t   ".join(first_line_lst)
            new_first_line += "\n"   # 'Energy\tTotal\tMo-s\tMo-p\tMo-s\tMo-d\tS-s\tS-p\n'
            
            ### Step 1.1. 删除原来的行首
            lines_lst.pop(0)        
            ### Step 1.2. 添加新的行首
            lines_lst.insert(0, new_first_line)
        
        ### Step 2. 写入新的 DOS.totalspin 到 `DOS.totalspin.bak`
        current_path = os.getcwd()
        filename = os.path.basename(self.dos_totalspin_path)
        bak_dos_totalspin_path = os.path.join(current_path, filename)
        with open(bak_dos_totalspin_path, "w") as f:
            f.writelines(lines_lst)
            
        ### Step 3. 读取 DOS.totalspin.bak 的数据
        ###     - row: 
        ###     - column: 
        tdos_array = np.loadtxt(
                        fname=bak_dos_totalspin_path,
                        #delimiter="\s+",
                        skiprows=1,     # 跳过第一行: 'Energy\tTotal\tMo-s\tMo-p\tMo-s\tMo-d\tS-s\tS-p\n'
                        )
        
        ### Step 4. 删除中间文件(temp file) -- DOS.totalspin.bak
        os.remove(bak_dos_totalspin_path)
        
        return tdos_array