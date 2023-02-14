'''
Author       : Liu Hanyu
Email        : hyliu2016@buaa.edu.cn
Date         : 2022-11-03 13:52:41
LastEditTime : 2022-11-03 14:27:44
FilePath     : /pflow/pflow/io/pwmat/test/test_addMag.py
Description  : 
'''
import unittest

# python3 -m pflow.io.pwmat.test.test_addMag
from ..addMag import AtomConfigMagTem


class AtomConfigMagTemTest(unittest.TestCase):
    def test_add_mag(self):
        atom_config_path_template = "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/atom_config/addMag/tem/atom.config"
        atom_config_path_sqs = "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/atom_config/addMag/output/atom.config"
        atomic_number_zero_lst = [12]
        
        atom_config_mag_tem = AtomConfigMagTem(
                                atom_config_path_sqs=atom_config_path_sqs,
                                atom_config_path_template=atom_config_path_template
                                )
        
        tmp = atom_config_mag_tem.assign_magnetic_moment(
                                atomic_number_zero_lst=atomic_number_zero_lst
                                )
        #print(tmp)
        atom_config_mag_tem.append_mag_to_atomconfig()


if __name__ == "__main__":
    unittest.main()