'''
Author       : Liu Hanyu
Email        : hyliu2016@buaa.edu.cn
Date         : 2022-10-31 16:24:03
LastEditTime : 2022-12-07 16:52:40
FilePath     : /pflow/pflow/io/publicLayer/structure.py
Description  : 
'''
import numpy as np
from pymatgen.core import Structure

from ..pwmat.utils.atomConfigExtractor import AtomConfigExtractor
from ..pwmat.utils.parameters import specie2atomic_number


class DStructure(Structure):
    '''
    Description
    -----------
        1. The Derivated class of `pymatgen.core.Structure`
        2. Overload the member method:
            1. from_file()
            2. to()
    
    New attributes
    --------------
        PWmat 的 atom.config 初始而成的 DStructure 多出以下两种 attributions
        1. self.atoms_lst
        2. self.atomic_numbers_lst
    '''
    @classmethod
    def from_file(cls,
                file_path:str,
                file_format:str,
                coords_are_cartesian:bool=False,
                ):
        '''
        Parameters
        ----------
            1. file_path: str
                结构文件路径，atom.config 或 POSCAR
            2. file_format: str
                1. "vasp"
                2. "pwmat"
                3. "cif"
                4. ...
            3. coords_are_cartesian: bool
                1. 坐标是否是笛卡尔形式，默认是分数形式
        
        Note
        ----
            1. Reads a structure from a file. For example, 
                anything ending in a `cif` is assumed to be a 
                Crystallographic Information Format file.
        '''
        if (file_format != "pwmat"):
            structure = Structure.from_file(filename=file_path)

        if (file_format == "pwmat"):
            atom_config_extractor = AtomConfigExtractor(
                                            atom_config_path=file_path,
                                            )

            structure = Structure(
                        lattice=atom_config_extractor.basis_vectors_array,
                        species=atom_config_extractor.species_array,
                        coords=atom_config_extractor.coords_array,
                        coords_are_cartesian=coords_are_cartesian,
                        site_properties={
                            "magmom": atom_config_extractor.magnetic_moments,
                            }
                        )

        structure.__class__ = cls
        return structure
        

    def to(self,
            output_file_path:str,
            output_file_format:str,
            include_magnetic_moments:bool=False,
            ):
        '''
        Desription
        ----------
            1. 将 Structure 对象输出成文件
        
        Parameters
        ----------
            1. output_file_path: str
                文件输出的路径
            2. 输出文件的格式
                1. "pwmat"
                2. "poscar" / "vasp"
                3. "cssr"
                4. "json"
                5. "xsf"
                6. "mcsqs"
                7. "prismatic"
                8. "yaml"
                9. "fleur-inpgen"
        '''
        if (output_file_format != "pwmat"):
            super(DStructure, self).to(
                                    fmt=output_file_format,
                                    filename=output_file_path)
        
        if (output_file_format == "pwmat"):
            with open(output_file_path, "w") as f:
                # 1. 
                f.write("  {0} atoms\n".format(self.num_sites))

                # 2. Lattice vector 信息
                f.write("Lattice vector (Angstrom)\n")
                f.write("   {0:<14E}    {1:<14E}    {2:<14E}\n".format(
                                            self.lattice.matrix[0, 0],
                                            self.lattice.matrix[0, 1],
                                            self.lattice.matrix[0, 2],
                                            )
                        )
                f.write("   {0:<14E}    {1:<14E}    {2:<14E}\n".format(
                                            self.lattice.matrix[1, 0],
                                            self.lattice.matrix[1, 1],
                                            self.lattice.matrix[1, 2],
                                            )
                        )
                f.write("   {0:<14E}    {1:<14E}    {2:<14E}\n".format(
                                            self.lattice.matrix[2, 0],
                                            self.lattice.matrix[2, 1],
                                            self.lattice.matrix[2, 2],
                                            )
                        )

                # 3. 
                f.write("Position (normalized), move_x, move_y, move_z\n")

                # 4. sites 的坐标信息
                for idx_site in range(self.num_sites):
                    f.write("  {0:<2d}         {1:<10f}         {2:<10f}         {3:<10f}     1  1  1\n".format(
                                    specie2atomic_number[str(self.species[idx_site])],
                                    self.frac_coords[idx_site, 0],
                                    self.frac_coords[idx_site, 1],
                                    self.frac_coords[idx_site, 2]
                                    )
                        )

                # 5. 向 atom.config 写入磁性信息
                if include_magnetic_moments:
                    f.write("Magnetic\n")
                    for idx_site in range(self.num_sites):
                        f.write("  {0:<3d} {1:<.2f}\n".format(
                                    specie2atomic_number[str(self.species[idx_site])],
                                    self.site_properties["magmom"][idx_site],
                                    )
                        )

    
    def judge_vacuum_exist(self):
        '''
        Description
        -----------
            1. 判断结构中是否存在 `真空层`
            2. structure.lattice.abc[-1] - (`原子最大z坐标 - 原子最小z坐标`) > 10
        
        Return
        ------
            1. vacuum_lst: list of bool
                - e.g. [True, True, False]: x方向有真空层，y方向有真空层，z方向没有真空层
        '''
        vacuum_lst = []

        for idx_direction in range(3):
            lattice_z_length = self.lattice.abc[idx_direction]
            
            coordination_z_lst = self.cart_coords[:, idx_direction]
            max_coordination_z = np.max(coordination_z_lst)
            min_coordination_z = np.min(coordination_z_lst)
            z_length = max_coordination_z - min_coordination_z
            
            if ( (lattice_z_length - z_length) > 10):
                vacuum_lst.append(True)
            
            else:
                vacuum_lst.append(False)
        
        return vacuum_lst
    
    
    def reformat_elements(self):
        '''
        Description
        -----------
            1. Reformat `DStructure` object in specified order of elements
                - 按照原子序数，从小到大排列
        '''
        self.sites.sort(key=lambda tmp_site: specie2atomic_number[str(tmp_site.specie)])

    
    def get_bidx2aidx_supercell(
                            self,
                            scaling_matrix:np.ndarray
                            ):
        '''
        Description
        -----------
            1. 
        
        Return
        ------
            1. bidx2aidx: Dict[int, int]
                - key: 扩胞前，primitive_cell 中原子对应的 index
                - value: 扩包后，primitive_cell 中原子对应的 index
                - e.g. 
                
        sorted_indexes    
        --------------
            [1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23, 25, 26, 28, 29, 31, 32, 34, 35, 37, 38, 40, 41, 43, 44, 46, 47, 49, 50, 52, 53, 55, 56, 58, 59, 61, 62, 64, 65, 67, 68, 70, 71, 73, 74, 76, 77, 79, 80, 82, 83, 85, 86, 88, 89, 91, 92, 94, 95, 97, 98, 100, 101, 103, 104, 106, 107, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105]\
        bidx2aidx
        ---------
            {1: 0, 2: 1, 4: 2, 5: 3, 7: 4, 8: 5, 10: 6, 11: 7, 13: 8, 14: 9, 16: 10, 17: 11, 19: 12, 20: 13, 22: 14, 23: 15, 25: 16, 26: 17, 28: 18, 29: 19, 31: 20, 32: 21, 34: 22, 35: 23, 37: 24, 38: 25, 40: 26, 41: 27, 43: 28, 44: 29, 46: 30, 47: 31, 49: 32, 50: 33, 52: 34, 53: 35, 55: 36, 56: 37, 58: 38, 59: 39, 61: 40, 62: 41, 64: 42, 65: 43, 67: 44, 68: 45, 70: 46, 71: 47, 73: 48, 74: 49, 76: 50, 77: 51, 79: 52, 80: 53, 82: 54, 83: 55, 85: 56, 86: 57, 88: 58, 89: 59, 91: 60, 92: 61, 94: 62, 95: 63, 97: 64, 98: 65, 100: 66, 101: 67, 103: 68, 104: 69, 106: 70, 107: 71, 0: 72, 3: 73, 6: 74, 9: 75, 12: 76, 15: 77, 18: 78, 21: 79, 24: 80, 27: 81, 30: 82, 33: 83, 36: 84, 39: 85, 42: 86, 45: 87, 48: 88, 51: 89, 54: 90, 57: 91, 60: 92, 63: 93, 66: 94, 69: 95, 72: 96, 75: 97, 78: 98, 81: 99, 84: 100, 87: 101, 90: 102, 93: 103, 96: 104, 99: 105, 102: 106, 105: 107}
        '''
        supercell = self.make_supercell_(
                            scaling_matrix=scaling_matrix,
                            reformat_mark=False)
        sorted_indexes = [
                        idx for idx, _ in \
                                sorted(
                                    enumerate(supercell.sites), 
                                    key=lambda tmp_entry: specie2atomic_number[str(tmp_entry[1].specie)]
                                    )
                        ]
        bidx2aidx = {sorted_indexes[i]: sorted_indexes.index(sorted_indexes[i]) \
                                        for i in range(len(sorted_indexes))}
        return bidx2aidx
    
    
    def remove_vacanies(self):
        '''
        Description
        -----------
            1. 删除结构中的空位
                - 空位的元素用 "X0+" 表示
        '''
        remove_indexes_lst = []
        for tmp_idx, site in enumerate(self.sites):
            if str(site.specie) == "X0+":
                remove_indexes_lst.append(tmp_idx)
    
        for tmp_idx in remove_indexes_lst:
            self.remove_sites(indices=remove_indexes_lst)
    
    
    def make_supercell_(self,
                    scaling_matrix: np.ndarray,
                    reformat_mark:bool=True):
        '''
        Description
        -----------
            1. 将自身扩包
        
        Parameters
        ----------
            1. scaling_matrix: np.array
            2. reformat_mark: bool
                - 是否按照原子顺序排序
        
        Note
        ----
            1. 扩胞倍数一定要是奇数，保证原胞在中心位置（便于构造特征等用途）
        
        Difference from pymatgen.core.Structure.make_supercell()
        --------------------------------------------------------
            1. 
        '''
        ### Step 1. 获取扩包后的晶格矢量(type=np.ndarray, shape=3*3)
        new_lattice_array = np.dot(
                                np.eye(3)*scaling_matrix, 
                                self.lattice.matrix)
        
        ### Step 2. 获取扩包后所有位点的坐标 (
        #               type=np.ndarray, 
        #               shape=(num_atoms * scaling_matrix[0]*scaling_matrix[1]*scaling_matrix[2], 3)
        #               )
        ### Step 2.1. 获取 `shift_matrix_frac`
        grid = np.meshgrid(
                        np.arange(scaling_matrix[0]),
                        np.arange(scaling_matrix[1]), 
                        np.arange(scaling_matrix[2]),
                        indexing="ij",
                        )
        ### grid[0]: shift_matrix_frac 所有点的 x 坐标
        ### grid[1]: shift_matrix_frac 所有点的 y 坐标
        ### grid[2]: shift_matrix_frac 所有点的 z 坐标
        '''
        shift_matrix_coeffs: np.ndarray
        -----------------
            [
                [-1. -1. -1.]
                [-1. -1.  0.]
                [-1. -1.  1.]
                [-1.  0. -1.]
                ...
            ]
        
        Note
        ----
            1. 需要让原胞仍处于supercell的中心（坐标系原点处）
        '''
        # Note: 需要让原胞仍处于supercell的中心（坐标系原点处）
        shift_matrix_coeffs = (
            np.vstack(
                    [grid[0].ravel() - (scaling_matrix[0]-1)/2,
                     grid[1].ravel() - (scaling_matrix[1]-1)/2,
                     grid[2].ravel() - (scaling_matrix[2]-1)/2])
            ).T
        # Note: 需要删除 np.array([0, 0, 0])
        # `np.any(a!=[0,0,0], axis=1)`: 全为 True，才返回 `True`
        mask = np.any(
                    shift_matrix_coeffs != np.array([0, 0, 0]),
                    axis=1
                    )

        shift_matrix_coeffs = shift_matrix_coeffs[mask] # (26, 3)
        # Note: 将 np.array([0, 0, 0]) 添加到第一个
        shift_matrix_coeffs = np.insert(shift_matrix_coeffs, 0, np.array([0,0,0]), axis=0) # (27, 3)
        
        ### Step 2.2. 获取平移后所有原胞的坐标信息
        '''
        tmp_shift_matrix_coeff
        ----------------------
            [1, 1, 0]
        
        tmp_shift_matrix_frac
        ---------------------
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 0]
            ]
            
        tmp_shift_matrix_cart = tmp_shift_matrix_frac * basis_vectors
        ---------------------
        
        tmp_shift_vector_cart = np.sum(tmp_shift_matrix_cart, axis=0)
        ---------------------
            1. 直接加到每个 site 的坐标上，就完成了 site 的平移
            2. 每个 tmp_shift_vector_cart 对应一个平移后的原胞
        
        tmp_pcell_coords_cart
        ---------------------
            1. 每个 tmp_shift_vector_cart 对应一个平移后的原胞
        '''
        pcell_coords_cart = []
        for tmp_shift_matrix_coeff in shift_matrix_coeffs:
            tmp_shift_matrix_frac = tmp_shift_matrix_coeff * np.eye(3)            
            tmp_shift_matrix_cart = np.dot(tmp_shift_matrix_frac, self.lattice.matrix)
            tmp_shift_vector_cart = np.sum(tmp_shift_matrix_cart, axis=0)      
            tmp_pcell_coords_cart = self.cart_coords + tmp_shift_vector_cart
            pcell_coords_cart.append(tmp_pcell_coords_cart)  # (27, num_atoms, 3)
        ### Step 2.3. `new_coords_cart`: supercell 所有位点的坐标
        ### new_coords_cart.shape = (
        #           scaling_matrix[0]*scaling_matrix[1]*scaling_matrix[2] * \
        #           num_atoms,
        #           3)
        new_coords_cart = np.array(pcell_coords_cart).reshape(-1, 3)
        #print(new_coords_cart.shape)    # (324, 3)
    
        ### Step 3. 获取扩包后所有位点的元素种类
        new_species = self.species * scaling_matrix[0] * scaling_matrix[1] * scaling_matrix[2]
        #print(len(new_species))
        #print(new_coords_cart.shape[0])
        assert len(new_species) == new_coords_cart.shape[0]
    
        ### Step 4. 用前三步得到的信息，初始化一个 DStructure 类
        supercell = DStructure(
                        lattice=new_lattice_array, # cartesian coordinates
                        species=new_species,
                        coords=new_coords_cart,  # cartesian coordinations
                        coords_are_cartesian=True
                        )

        ### Step 5. 是否按照原子序数，从小到大排序
        if reformat_mark:
            supercell.reformat_elements()
            
        return supercell