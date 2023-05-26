import unittest

# python3 -m matersdk.feature.avg.test.test_avgbond
from ..avgbond import AvgBond, PairBond


class AvgBondTest(unittest.TestCase):
    def test_all(self):
        #movement_path = "/data/home/liuhanyu/hyliu/data_for_test/bondfft_GeTe/MOVEMENT"
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
        #element_1 = "Ge"
        #element_2 = "Te"
        element_1 = "Li"
        element_2 = "Si"
        rcut = 3.2
        
        
        bondfft = AvgBond(
                        movement_path=movement_path,
                        element_1=element_1,
                        element_2=element_2,
                        rcut=rcut)
        ### Step 1.
        print( len(bondfft.frames_lst) )
        
        ### Step 2. 
        print(bondfft.get_frames_avg_bond())
        
    
    def all(self):
        #movement_path = "/data/home/liuhanyu/hyliu/data_for_test/bondfft_GeTe/MOVEMENT"
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
        atom1_idx = 0
        atom2_idx = 1
        
        pairbond = PairBond(
                            movement_path=movement_path,
                            atom1_idx=atom1_idx,
                            atom2_idx=atom2_idx)
        #print( len(pairbond.frames_lst) )
        
        #print(pairbond.get_frames_pair_bond())
        


if __name__ == "__main__":
    unittest.main()