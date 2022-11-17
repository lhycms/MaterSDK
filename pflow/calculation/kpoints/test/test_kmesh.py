'''
Author       : Liu Hanyu
Email        : hyliu2016@buaa.edu.cn
Date         : 2022-11-17 09:54:56
LastEditTime : 2022-11-17 12:07:25
FilePath     : /pflow/pflow/calculation/kpoints/test/test_kmesh.py
Description  : 
'''
import unittest

# python3 -m pflow.calculation.kpoints.test.test_kmesh
from ..kmesh import KMesh


class KMeshTest(unittest.TestCase):
    def test_get_lattice_info(self):
        file_format = "pwmat"
        file_path = "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/atom_config/atom.config"

        kmesh = KMesh(
            file_format=file_format,
            file_path=file_path
            )
        kmesh.get_lattice_info()
    


    def test_get_kmesh(self):
        file_format = "vasp"
        file_path = "/Users/mac/Desktop/ReNbSSe/mc_2/0/POSCAR"
        file_format = "pwmat"
        file_path = "/Users/mac/我的文件/Mycode/new/new2/pflow/test_data/atom_config/atom.config"
        density = 0.03

        kmesh = KMesh(
            file_format=file_format,
            file_path=file_path
        )
        print("KMesh when density = {0} (unit: 2pi/Angstrom)".format(density))
        print(kmesh.get_kmesh(density))


if __name__ == "__main__":
    unittest.main()