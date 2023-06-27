import unittest
import numpy as np

# python3 -m matersdk.feature.avg.test.test_msd
from ..msd import (
                Msd,
                MsdParallelFunction
)
from ....io.pwmat.output.movement import Movement



class MsdTest(unittest.TestCase):
    def test_all(self):
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
        #movement_path = "/data/home/liuhanyu/hyliu/pwmat_demo/xhm/MOVEMENT"
        movement = Movement(movement_path=movement_path)
        
        ### Step 1. calc_msd()
        print()
        print("Step 1. calc_msd():")
        msd_object = Msd(trajectory=movement)
        msd_values_lst = msd_object.calc_msd()
        print(msd_values_lst)


class ParallelFunctionTest(unittest.TestCase):
    def all(self):
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
        movement = Movement(movement_path=movement_path)
        structure_1 = movement.get_frame_structure(idx_frame=0)
        structure_2 = movement.get_frame_structure(idx_frame=100)
        
        msd_value = MsdParallelFunction.calc_msd_s(
                            structure_1=structure_1,
                            structure_2=structure_2)
        print("msd_value: ", msd_value)
        


if __name__ == "__main__":
    unittest.main()