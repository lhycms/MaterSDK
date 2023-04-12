'''
Author       : Liu Hanyu
Email        : hyliu2016@buaa.edu.cn
Date         : 2022-10-31 16:45:10
LastEditTime : 2022-12-09 11:54:28
FilePath     : /pflow/pflow/io/publicLayer/test/test_Structure.py
Description  : 
'''
import unittest
import numpy as np

# python3 -m pflow.io.publicLayer.test.test_Structure
from ..structure import DStructure



class StructureTest(unittest.TestCase):
    def test_structure_init(self):
        file_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
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
        file_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        
        # 1. 
        output_file_path = \
                "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/test.config"
        # 2. 
        output_file_path = \
                "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/test.config"
        coords_are_cartesian = False

        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )

        #print(structure.site_properties["magmom"])

        # 1.
        #structure.to(output_file_path=output_file_path,
        #            output_file_format="rndstr.in")
        
        # 2. 
        structure.to(output_file_path=output_file_path,
                    output_file_format="pwmat",
                    include_magnetic_moments=True,
                   )


    def test_judge_vacuum_exist(self):
        file_format = "pwmat"
        file_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        coords_are_cartesian = False

        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )

        #print(structure.judge_vacuum_exist())
    
    
    def test_reformat_elements(self):
        file_format = "pwmat"
        file_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        coords_are_cartesian = False
        elements_lst = ["S", "Mo"]
        
        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )
        
        print("Step 1. Primimary structure:")
        print(structure)
        new_structure = structure.reformat_elements(elements_lst=elements_lst)
        print("\nStep 2. New structure:")
        print(new_structure)  
    
    
    def test_remove_vacanies(self):
        file_format = "pwmat"
        file_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        coords_are_cartesian = False
        
        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )
        
        #print("删除体系内所有空位前，体系内的原子数：{0}".format(
        #                                len(structure.atomic_numbers)
        #                                )
        #      )
        #structure.remove_vacanies()
        #print("删除体系内所有空位前，体系内的原子数：{0}".format(
        #                                len(structure.atomic_numbers)
        #                                )
        #      )
        #print("第 3 个 site 处的magmom: {0}".format(structure.sites[2].magmom)
        #      )
        #print(structure)
        #structure.make_supercell([2,2,2])
        #print(structure)
    
    
    def test_make_supercell(self):
        file_format = "pwmat"
        # 1. 普通的测试
        file_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        coords_are_cartesian = False
        
        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )
        scaling_matrix = np.array([3,3,3])
        
        supercell = structure.make_supercell_(
                                scaling_matrix=scaling_matrix,
                                )
        #print(supercell)
        #print(supercell)
        #supercell.to(output_file_path="./POSCAR",
        #            output_file_format="vasp",
        #            include_magnetic_moments=True,
        #           )
        
        
    def test_get_bidx2aidx_supercell(self):
        file_format = "pwmat"
        # 1. 普通的测试
        file_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/structure/atom.config"
        coords_are_cartesian = False
        scaling_matrix=np.array([3, 3, 1])
        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )
        bidx2aidx = structure.get_bidx2aidx_supercell(scaling_matrix=scaling_matrix)
        #print(bidx2aidx)
        


if __name__ == "__main__":
    unittest.main()