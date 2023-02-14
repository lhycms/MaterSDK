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
        ### Part I. Get DStructure from atom.config
        file_format = "pwmat"
        file_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/kpath/atom.config"
        coords_are_cartesian = False
        structure = DStructure.from_file(
                        file_path=file_path,
                        file_format=file_format,
                        coords_are_cartesian=coords_are_cartesian
                        )
        
        ### Part II. Setting torlenrance
        symprec = 0.1
        angle_tolerance = 5
        atol = 1e-5
        density = 0.01

        ### Part III. test
        kpath_sampler = KpathSampler(
                            structure=structure,
                            symprec=symprec,
                            angle_tolerance=angle_tolerance,
                            atol=atol,
                            )
        kpath_sampler.HIGHK_file_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/kpath/HIGHK"
        kpath_sampler.gen_kpt_path = "/data/home/liuhanyu/hyliu/code/pflow/demo/kpath/gen.kpt"
        
        print(kpath_sampler.kpoints)
        print(kpath_sampler.kpaths)
        kpath_sampler.output_HIGHK_file()
        kpath_sampler.output_gen_kpt(density=density)



if __name__ == "__main__":
    unittest.main()