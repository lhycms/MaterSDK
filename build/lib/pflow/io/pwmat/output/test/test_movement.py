import unittest

# python3 -m pflow.io.pwmat.output.test.test_movement
from ..movement import Movement



class MovementTest(unittest.TestCase):
    def test_all(self):
        #movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo1/PWdata/data1/MOVEMENT"
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
        movement = Movement(
                                movement_path=movement_path)
        idx_frame = 100 # 帧数从 0 开始计数
        
        # 1. get_chunksize()
        print()
        print("Step 1. The chunksize of each frame:", end="\t")
        print(movement.get_chunksize())
        
        # 2. Get chunk of specific frame structure
        print()
        print("Step 2. String for structure:")
        print(movement._get_frame_str(idx_frame=idx_frame))
        
        # 3. 
        print()
        print("Step 3. Structure from {0}th frame:".format(idx_frame))
        print(movement.get_frame_structure(idx_frame=idx_frame))
        
        # 4. 
        print()
        print("Step 4. Virial tensor of {0}th frame:".format(idx_frame)) 
        frame_structure = movement.get_frame_structure(idx_frame=idx_frame)
        print(frame_structure.virial_tensor)   
        
        # 5.
        print()
        print("Step 5. Energy of {0}th")
        print(movement.get_energy(idx_frame=idx_frame))
    
        
        
if __name__ == "__main__":
    unittest.main()