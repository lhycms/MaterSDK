import unittest

# python3 -m matersdk.data.deepmd.test.test_data_system
from ..data_system import DeepmdDataSystem, SubDeepmdDataSystem
from ....io.pwmat.output.movement import Movement


class DeepmdDataSystemTest(unittest.TestCase):
    def all(self):
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
        assert (dp_data_system.num_structures == len(dp_data_system))
        print(dp_data_system.num_structures, " = ", len(dp_data_system))
        
        ### Step 2. 
        #print()
        #print("Step 2. DeepDataSystem.structures_lst:")
        #print(dp_data_system.structures_lst)
        
        ### Step 3. 
        #print()
        #print("Step 3. DeepDataSystem.total_energys_lst:")
        #print(dp_data_system.total_energys_array)
        
        
        ### Step 4.
        #print()
        #print("Step 4. Save information: running...")
        #dp_data_system.save_all_info(
        #                dir_path=dir_path,
        #                scaling_matrix=scaling_matrix)
        
        
        ### Step 5.
        #print()
        #print("Step 5. The atomic number in system:", end='\t')
        #print(dp_data_system.atomic_numbers_lst)
        
        
        ### Step 6.
        #print()
        #print("Step 6. The number of atoms in system:", end='\t')
        #print(dp_data_system.num_atoms)
        
        
        ### Step 7.
        print()
        print("\nStep 7. __getitem__():")
        sub_dp_data_system = dp_data_system[:20]
        print("\nStep 7.1. len(sub_dp_data_system):", end='\t')
        print(len(sub_dp_data_system))


class SubDeepmdDataSystemTest(unittest.TestCase):
    def test_all(self):
        movement = Movement(movement_path="/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT")
        rcut = 6.5
        
        dp_data_system = DeepmdDataSystem.from_trajectory_s(
                            trajectory_object=movement,
                            rcut=rcut)
        sub_dp_data_system = SubDeepmdDataSystem.from_indices(
                                deepmd_data_system=dp_data_system,
                                indices_lst=[*range(104)])
        ### Step 1.
        print()
        print("Step 1. Extract sub_deepmd_data_system, including {0} structures".format(len(sub_dp_data_system)))

if __name__ == "__main__":
    unittest.main()