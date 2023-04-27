import unittest

# python3 -m matersdk.data.deepmd.test.test_data_system
from ..data_system import DeepmdDataSystem
from ....io.pwmat.output.movement import Movement


class DeepmdDataSystemTest(unittest.TestCase):
    def test_from_trajectory(self):
        movement = Movement(movement_path="/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT")
        dp_data_system = DeepmdDataSystem.from_trajectory(trajectory_object=movement)
        
        ### Step 1.
        print()
        print("Step 1. DeepDataSystem 包含的结构数目:", end='\t')
        print(dp_data_system.num_structures)
        
        ### Step 2. 
        print()
        print("Step 2. DeepDataSystem.structures_lst:")
        print(dp_data_system.structures_lst)
        
        ### Step 3. 
        print()
        print("Step 3. DeepDataSystem.total_energys_lst:")
        print(dp_data_system.total_energys_array)
        


if __name__ == "__main__":
    unittest.main()