import os
import numpy as np 
from typing import Dict, List
from pymatgen.symmetry.kpath import KPathSetyawanCurtarolo

from ...io.publicLayer.structure import DStructure


class KpathSampler(KPathSetyawanCurtarolo):
    HIGHK_file_path = "./HIGHK" # 输出的 HIGHK 路径
    gen_kpt_path = "./gen.kpt"  # 输出的 gen.kpt 路径，用于 split_kp.x

    def __init__(self,
                structure:DStructure,
                dimension:int,
                symprec:float=0.1,
                angle_tolerance:float=5,
                atol:float=1e-5,
                ):
        '''
        Parameters
        ----------
            1. structure: DStructure
            2. dimension: Optionnal[2|3]
                - 二维材料 / 三维材料
            3. symprec: float
                - 对称性误差
            4. angle_tolerance: float
                - 角度误差
            5. atol: float
                - 总误差
        '''
        super(KpathSampler, self).__init__(
                    structure=structure,
                    symprec=symprec,
                    angle_tolerance=angle_tolerance,
                    atol=atol,
                    )
        self.dimension = dimension

        ### Note: 必须先得到 self.kpoints，再得到 self.kpaths (利用 self.points 的坐标处理二维材料)
        # 字典：{"高对称点": "坐标", ...}
        self.kpoints = self._get_kpoints()
        # e.g. [['G', 'M', 'K', 'G', 'A', 'L', 'H', 'A'], ['L', 'M'], ['K', 'H']]
        self.kpaths = self._get_kpaths()

    
    @property
    def kpoints(self) -> Dict[str, np.ndarray]:
        return self._kpoints
    
    @kpoints.setter
    def kpoints(self, kpoints:Dict[str, np.ndarray]):
        self._kpoints = kpoints
        
    @property
    def kpaths(self) -> List[List[str]]:
        return self._kpath
    
    @kpaths.setter
    def kpaths(self, kpath:List[List[str]]):
        self._kpath = kpath
        
        
    @property
    def mark_vacuum_z(self) -> bool:
        return self._mark_vacuum_z
    
    @mark_vacuum_z.setter
    def mark_vacuum_z(self, mark: bool):
        self._mark_vacuum_z = mark
    
    
    def output_HIGHK_file(self):
        '''
        Description
        -----------
            1. 根据 self.kpoints(dict)，生成 HIGHK 文件，存储高对称点和高对称点的坐标
                High symmetry kpoints generated by Q-Flow, structure:Mo1 S2, fractional coordinates in reciprocal lattice
                    0.000000     0.000000     0.000000    G
                    0.000000     0.000000     0.500000    A
                    0.333333     0.333333     0.500000    H
                    0.333333     0.333333     0.000000    K
                    0.500000     0.000000     0.500000    L
                    0.500000     0.000000     0.000000    M        
        '''
        if os.path.exists(self.HIGHK_file_path):
            os.remove(self.HIGHK_file_path)
        
        with open(self.HIGHK_file_path, "a") as f:
            f.write("High symmetry kpoints generated by MaterSDK, structure:{0}, fractional coordinates in reciprocal lattice\n".format(self.structure.formula))

            for tmp_kpoint in self.kpoints.keys():
                tmp_coord = self.kpoints[tmp_kpoint]
                f.write("    {0:>12.6f}\t{1:>12.6f}\t{2:>12.6f}\t{3}\n".format(
                                            tmp_coord[0],
                                            tmp_coord[1],
                                            tmp_coord[2],
                                            tmp_kpoint,
                                            )
                )


    def output_gen_kpt(self, density:float=0.01):
        '''
        Description
        -----------
            1. 根据 self.kpaths 生成 gen.kpt 文件
                - e.g.
                -   self.kpaths = [['G', 'M', 'K', 'G', 'A', 'L', 'H', 'A'], ['L', 'M'], ['K', 'H']]
        '''
        if os.path.exists(self.gen_kpt_path):
            os.remove(self.gen_kpt_path)

        '''
        Note
        ----
            1. No `pymatgen.cores.Lattice.reciprocal_lattic_crystalloggraphic`: no factor of 2*pi
            2. The `density` user input has 2*pi, but doesn't deal with in code, 因为和 1. 约分了 (同除以 2*pi)
        '''
        reciprocal_lattice_crystalloggraphic = \
                self.structure.lattice.reciprocal_lattice_crystallographic.matrix   # no factor of 2*pi


        with open(self.gen_kpt_path, "a") as f:
            f.write("K-path generated by MaterSDK, structure:{0}, density:{1}, fractional coordinates in reciprocal lattice\n".format(self.structure.formula, density))
            #print("*"*40)
            #print(self.kpaths)
            #print("*"*40)
            for idx_kpath, tmp_kpath in enumerate(self.kpaths):
                for idx_point, tmp_point in enumerate(tmp_kpath):
                    if (idx_point == 0):
                        point_start = tmp_point
                        point_dest = tmp_point
                        position_dest = self.kpoints[tmp_point]
                    else:
                        point_start = point_dest
                        position_start = position_dest

                        point_dest = tmp_point
                        position_dest = self.kpoints[tmp_point]
                        
                        position_start_array = np.array(position_start)
                        position_dest_array = np.array(position_dest)
                        path_length = np.linalg.norm(
                                (position_dest_array - position_start_array).reshape(-1, 3) * \
                                        reciprocal_lattice_crystalloggraphic
                                )
                        num_bandpoint = np.ceil(path_length/density)
                    
                        f.write("{0}\n".format(num_bandpoint))
                        f.write("    {0:>12.6f}\t{1:>12.6f}\t{2:>12.6f}\t{3}\n".format(
                                    position_start[0],
                                    position_start[1],
                                    position_start[2],
                                    point_start
                                    )
                        )
                        f.write("    {0:>12.6f}\t{1:>12.6f}\t{2:>12.6f}\t{3}\n".format(
                                    position_dest[0],
                                    position_dest[1],
                                    position_dest[2],
                                    point_dest
                                    )
                        )

                    
    def _get_kpoints(self):
        '''
        Description
        ------------
            1. 获取 `Structure` 所有高对称点的名字和坐标
        
        Return 
        ------
            1. {高对称点(str): 坐标(np.array)}
            e.g.
                {'\\Gamma': array([0., 0., 0.]), 'A': array([0. , 0. , 0.5]), 'H': array([0.33333333, 0.33333333, 0.5       ]), 'K': array([0.33333333, 0.33333333, 0.        ]), 'L': array([0.5, 0. , 0.5]), 'M': array([0.5, 0. , 0. ])}
        
        Note
        ----
            1. 最后将 "\\Gamma" 置换成 "G"
        '''
        ### Step 1. 将 "\\Gamma" 置换成 "G"
        kpoint2coord = self.kpath["kpoints"]
        try:
            kpoint2coord.update({'G': kpoint2coord.pop('\\Gamma')})
        except:
            pass
        
        ### Step 2. 如果是二维材料，需要把 `z 坐标不为0的 Kpoint 删去`
        return_kpoint2coord = {}
        if (self.dimension == 2):   # 把 z 坐标不为0的 Kpoint 删去
            for tmp_kpt_name, tmp_kpt_coord in kpoint2coord.items():
                if tmp_kpt_coord[-1] == 0:
                    return_kpoint2coord.update({tmp_kpt_name: tmp_kpt_coord})
            return return_kpoint2coord
        
        elif (self.dimension == 3):
            return kpoint2coord
        
        else:
            raise Exception("The dimension of material must be 2 or 3!")


    def _get_kpaths(self):
        '''
        Description
        ------------
            1. 获取 `Structure` 所有对称路径

        Return 
        ------
            1. 多条高对称路径 (多个列表)
            e.g.
                [['G', 'M', 'K', 'G', 'A', 'L', 'H', 'A'], ['L', 'M'], ['K', 'H']]
        
        Note
        ----
            1. 最后将 "\\Gamma" 置换成 "G"
        '''
        ### Step 1. 将所有路径中的 "\\Gamma" 置换成 "G"
        new_kpaths_lst = []

        for tmp_kpath in self.kpath["path"]:
            if "\\Gamma" in tmp_kpath:
                tmp_kpath_new = ["G" if tmp_point=="\\Gamma" else tmp_point \
                                    for tmp_point in tmp_kpath]
                new_kpaths_lst.append(tmp_kpath_new)
            else:
                new_kpaths_lst.append(tmp_kpath)

        ### Step 2. 如果是二维材料，需要把 `z 坐标不为0的 Kpoint 删去`
        if (self.dimension == 2):
            return_kpaths_lst = []
            ### Step 2.1. 删除 self.kpoints 中不存在的点
            for tmp_kpath in new_kpaths_lst:
                tmp_new_kpath = [tmp_kpt for tmp_kpt in tmp_kpath if tmp_kpt in self.kpoints.keys()]
                return_kpaths_lst.append(tmp_new_kpath)
            tmp_new_kpath = [tmp_kpt for tmp_kpt in tmp_kpath if tmp_kpt in self.kpoints.keys()]
                        
            ### Step 2.2. 删除`长度<=1`的 kpath: List(str)
            return_kpaths_lst_ = [tmp_kpath for tmp_kpath in return_kpaths_lst if (len(tmp_kpath)>1)] 
            
            return return_kpaths_lst_

        elif (self.dimension == 3):
            return new_kpaths_lst
        
        else:
            raise Exception("The dimension of material must be 2 or 3!")