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