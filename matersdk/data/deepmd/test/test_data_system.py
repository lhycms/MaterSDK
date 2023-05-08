import unittest

# python3 -m matersdk.data.deepmd.test.test_data_system
from ..data_system import DeepmdDataSystem
from ....io.pwmat.output.movement import Movement


class DeepmdDataSystemTest(unittest.TestCase):
    def test_from_trajectory(self):
        movement = Movement(movement_path="/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT")
        rcut = 6.5
        scaling_matrix = [3, 3, 3]
        
        dp_data_system = DeepmdDataSystem.from_trajectory_s(
                            trajectory_object=movement,
                            rcut=rcut)
        dir_path = "/data/home/liuhanyu/hyliu/code/test"
        
        ### Step 1.
        print()
        print("Step 1. DeepDataSystem 包含的结构数目:", end='\t')
        print(dp_data_system.num_structures)
        
        ### Step 2. 
        #print()
        #print("Step 2. DeepDataSystem.structures_lst:")
        #print(dp_data_system.structures_lst)
        
        ### Step 3. 
        #print()
        #print("Step 3. DeepDataSystem.total_energys_lst:")
        #print(dp_data_system.total_energys_array)
        
        
        ### Step 4.
        print()
        print("Step 4. Save information: running...")
        dp_data_system.save(
                        dir_path=dir_path,
                        scaling_matrix=scaling_matrix)
        
        
        ### Step 5.
        print()
        #print("Step 5. The atomic number in system:", end='\t')
        #print(dp_data_system.atomic_numbers_lst)
        
        
        ### Step 6.
        print()
        #print("Step 6. The number of atoms in system:", end='\t')
        #print(dp_data_system.num_atoms)


if __name__ == "__main__":
    unittest.main()