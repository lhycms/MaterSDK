import unittest

# python3 -m pflow.io.pwmat.output.test.test_movement
from ..movement import MovementExtractor



class MovementExtractorTest(unittest.TestCase):
    def test_all(self):
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo1/PWdata/data1/MOVEMENT"
        movement_extractor = MovementExtractor(
                                movement_path=movement_path)
        idx_frame = 998   # 帧数从 1 开始计数
        
        # 1. get_chunksize()
        print()
        print("Step 1. The chunksize of each frame: ", end="\t")
        print(movement_extractor.get_chunksize())
        
        # 2. 
        print()
        print("Step 2. ")
        print(movement_extractor._get_frame_str(idx_frame=idx_frame))
        
        # 3. 
        print()
        print("Step 3.")
        print(movement_extractor.get_frame_structure(idx_frame=idx_frame))
        
        # 4. 
        print()
        print("Step 4.")    
    
        
        
if __name__ == "__main__":
    unittest.main()