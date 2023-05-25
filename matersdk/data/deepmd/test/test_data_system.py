import unittest

# python3 -m matersdk.data.deepmd.test.test_data_system
from ..data_system import DpLabeledSystem
from ....io.pwmat.output.movement import Movement


class DpLabeledSystemTest(unittest.TestCase):
    def test_all(self):
        movement = Movement(movement_path="/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT")
        rcut = 6.5
        scaling_matrix = [3, 3, 3]
        
        dp_labeled_system = DpLabeledSystem.from_trajectory_s(
                            trajectory_object=movement,
                            rcut=rcut)
        dir_path = "/data/home/liuhanyu/hyliu/code/test"
        
        ### Step 1.
        print()
        print("Step 1. DeepDataSystem 包含的结构数目:", end='\t')
        assert (dp_labeled_system.num_structures == len(dp_labeled_system))
        print(dp_labeled_system.num_structures, " = ", len(dp_labeled_system))
        print(dp_labeled_system)
        
        ### Step 2. save_info_dpdata()
        print()
        print("Step 2. save_info_dpdata:")
        dp_labeled_system.save_info_dpdata(dir_path=dir_path,set_size=550)

        
        ### Step 3. from_indices()
        print()
        training_set = DpLabeledSystem.from_indices(
                                    dp_labeled_system=dp_labeled_system,
                                    indices_lst=[*range(400)])
        validation_set = DpLabeledSystem.from_indices(
                                    dp_labeled_system=dp_labeled_system,
                                    indices_lst=[*range(400,550)])
        print("Step 3. from_indices():")
        print(training_set)
        print(validation_set)
        
        
        ### Step 4.
        #print()
        #print("Step 4. Save information: running...")
        #dp_data_system.save_all_info(
        #                dir_path=dir_path,
        #                scaling_matrix=scaling_matrix)
        
        


if __name__ == "__main__":
    unittest.main()