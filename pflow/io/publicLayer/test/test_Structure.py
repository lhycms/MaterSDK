'''
Author       : Liu Hanyu
Email        : hyliu2016@buaa.edu.cn
Date         : 2022-10-31 16:45:10
LastEditTime : 2022-12-07 16:48:21
FilePath     : /pflow/pflow/io/publicLayer/test/test_Structure.py
Description  : 
'''
import unittest
import numpy as np

# python3 -m pflow.io.publicLayer.test.test_Structure
from ..structure import DStructure



class StructureTest(unittest.TestCase):
    def test_structure_init(self):
        file_path = "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/atom_config/atom.config"
        file_format = "pwmat"
        coords_are_cartesian = False
        
        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )
        #print(structure.cart_coords)
        #print(np.max(structure.cart_coords[:, -1]))
        #print(np.min(structure.cart_coords[:, -1]))
        #print(structure.atoms_lst)
        #print(structure.atomic_numbers_lst)
    

    def test_structure_to_(self):
        file_format = "pwmat"
        # 1. 普通的测试
        file_path = "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/atom_config/atom.config"
        # 2. 测试 DStructure.site_properties["magnetic_properties"]
        file_path = "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/atom_config/addMag/output/atom.config"
        
        # 1. 
        output_file_path = \
                "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/atom_config/rnstr.in"
        # 2. 
        output_file_path = \
                "./atom.config"
        coords_are_cartesian = False

        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )
        
        print(structure.site_properties["magnetic_moments"])

        # 1.
        #structure.to(output_file_path=output_file_path,
        #            output_file_format="rndstr.in")
        
        # 2. 
        structure.to(output_file_path=output_file_path,
                    output_file_format="pwmat",
                    include_magnetic_moments=True,
                    )



    def test_judge_vacuum_exist(self):
        file_format = "vasp"
        file_path = "/Users/mac/Desktop/ReNbSSe/mc_2/0/POSCAR"
        coords_are_cartesian = False

        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )

        #print(structure.judge_vacuum_exist())
        



if __name__ == "__main__":
    unittest.main()