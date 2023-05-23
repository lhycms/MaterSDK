import unittest

# python3 -m matersdk.toolkit.pwmat.test.test_bondfft
from ..bondfft import Bondfft


class BondfftTest(unittest.TestCase):
    def test_all(self):
        #movement_path = "/data/home/liuhanyu/hyliu/data_for_test/bondfft_GeTe/MOVEMENT"
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
        #element_1 = "Ge"
        #element_2 = "Te"
        element_1 = "Li"
        element_2 = "Si"
        cutoff = 3.2
        
        
        bondfft = Bondfft(
                        movement_path=movement_path,
                        element_1=element_1,
                        element_2=element_2,
                        cutoff=cutoff)
        ### Step 1.
        print( len(bondfft.frames_lst) )
        
        
        ### Step 2. 
        print(bondfft.get_frames_avg_bond())
        


if __name__ == "__main__":
    unittest.main()