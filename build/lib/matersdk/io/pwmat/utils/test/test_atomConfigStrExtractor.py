import unittest

# python3 -m matersdk.io.pwmat.utils.test.test_atomConfigStrExtractor
from ...output.movement import Movement
from ..atomConfigExtractor import AtomConfigStrExtractor


class AtomConfigStrExtractorTest(unittest.TestCase):
    def test_all(self):
        #movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo1/PWdata/data1/MOVEMENT"
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
        movement = Movement(movement_path=movement_path)
        idx_frame = 2 # 帧数从 0 开始计数
        
        ### Step 0. 得到某 frame 的 string
        atom_config_string = movement._get_frame_str(idx_frame=idx_frame)
        #print(atom_config_string)
        
        ### Step 1. 得到体系内的原子数目
        atom_config_str_extractor = AtomConfigStrExtractor(
                                        atom_config_str=atom_config_string)
        print()
        print("Step 1. The number of atom in system:", end="\t")
        print(atom_config_str_extractor.get_num_atoms())
        
        ### Step 2. 得到体系的 basis vectors
        print()
        print("Step 2. The basis vectors of system:") 
        print(atom_config_str_extractor.get_basis_vectors_lst())
        
        
        ### Step 3. 得到体系的维里张量
        print()
        print("Step 3. The virial tensor of system:")
        print(atom_config_str_extractor.get_virial_tensor())
        
        
        ### Step 4. 得到体系的原子序数 (重复的)
        print()
        print("Step 4. The atomic number in system:")
        print(atom_config_str_extractor.get_atomic_numbers_lst())
        
        
        ### Step 5. 得到体系的坐标 (np.array 形式)
        print()
        print("Step 5. The coords of atoms in system:")
        print(atom_config_str_extractor.get_coords_lst())
        
        
        ### Step 6. 得到体系中各个原子的受力
        print()
        print("Step 6. The atomic forces of atoms in system:")
        print(atom_config_str_extractor.get_atomic_forces_lst())
        
        
        ### Step 7. 得到体系中各个原子的速度
        print()
        print("Step 7. The atomic velocitys of atoms in system:")
        print(atom_config_str_extractor.get_atomic_velocitys_lst())
        
        
        ### Stpe 8. 得到体中各个原子的能量
        print()
        print("Step 8. The atomic enegy of atoms in system:")
        print(atom_config_str_extractor.get_atomic_energys_lst())
        

        ### Stpe 9. 得到体中各个原子的磁矩
        print()
        print("Step 8. The magnetic moment of atoms in system:")
        print(atom_config_str_extractor.get_magnetic_moments())
        


if __name__ == "__main__":
    unittest.main()