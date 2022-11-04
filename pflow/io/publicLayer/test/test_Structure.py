'''
Author       : Liu Hanyu
Email        : hyliu2016@buaa.edu.cn
Date         : 2022-10-31 16:45:10
LastEditTime : 2022-11-03 13:30:53
FilePath     : /pflow/pflow/io/publicLayer/test/test_Structure.py
Description  : 
'''
import unittest

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

        #print(structure)
        #print(structure.atoms_lst)
        #print(structure.atomic_numbers_lst)
    

    def test_structure_to_(self):
        file_format = "mcsqs"
        file_path = "/Users/mac/我的文件/Mycode/new/new2/pflow/rndstr.in"
        output_file_path = \
                "/Users/mac/我的文件/Mycode/new/new2/pflow/atom.config"
        coords_are_cartesian = False

        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )

        structure.to(output_file_path=output_file_path,
                    output_file_format="pwmat")
        



if __name__ == "__main__":
    unittest.main()