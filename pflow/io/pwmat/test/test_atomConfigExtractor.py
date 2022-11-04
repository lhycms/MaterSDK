'''
Author       : Liu Hanyu
Email        : hyliu2016@buaa.edu.cn
Date         : 2022-10-31 15:48:42
LastEditTime : 2022-11-03 14:31:53
FilePath     : /pflow/pflow/io/pwmat/test/test_atomConfigExtractor.py
Description  : 
'''
import unittest

# python3 -m pflow.io.pwmat.test.test_atomConfigExtractor
from ..atomConfigExtractor import AtomConfigExtractor


class AtomConfigExtractorTest(unittest.TestCase):
    def test_extractor(self):
        atom_config_path = "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/atom_config/addMag/tem/atom.config"
        atom_config_extractor = AtomConfigExtractor(
                                    atom_config_path=atom_config_path
                                    )
        #print(atom_config_extractor.num_atoms)
        #print(atom_config_extractor.basis_vectors_array)
        #print(atom_config_extractor.species_array)
        #print(atom_config_extractor.coords_array.shape)
        print(atom_config_extractor.get_atoms_lst())


if __name__ == "__main__":
    unittest.main()