'''
Author       : Liu Hanyu
Email        : hyliu2016@buaa.edu.cn
Date         : 2022-12-08 17:58:14
LastEditTime : 2022-12-09 11:53:36
FilePath     : /pflow/pflow/calculation/kpath/test/test_kpathSampler.py
Description  : 
'''
import unittest
import warnings

# python3 -m pflow.calculation.kpath.test.test_kpathSampler
from ....io.publicLayer.structure import DStructure
from ..kpathSampler import KpathSampler


warnings.filterwarnings("ignore")


class KpathSamplerTest(unittest.TestCase):
    def test_get_kpoints(self):
        pass



    def test_get_kpath(self):
        ### Part I. Structure Config
        file_format = "pwmat"
        #file_path = "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/demo/atom.config"
        file_path = "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/atom_config/高对称点/MoS2_pri/atom.pwmat"
        coords_are_cartesian = False

        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )
        

        ### Part II. Symmetry Config
        symprec = 0.1
        angle_tolerance = 5
        atol = 1e-5

        density = 0.01


        ### Part III. file path config


        ### Part IV. test
        kpath_sampler = KpathSampler(
                            structure=structure,
                            symprec=symprec,
                            angle_tolerance=angle_tolerance,
                            atol=atol,
                            )
        
        
        #print(kpath_sampler.get_kpath())
        #print(kpath_sampler.get_kpoints())
        #kpath_sampler.output_HIGHK_file()
        #kpath_sampler.output_gen_kpt(density=density)



if __name__ == "__main__":
    unittest.main()