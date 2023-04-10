import unittest

# python3 pflow.io.pwmat.output.test.test_movement
from ..movement import MovementExtractor



class MovementExtractorTest(unittest.TestCase):
    def test_all(self):
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo1/PWdata/data1/MOVEMENT"
        movement_extractor = MovementExtractor(
                                movement_path=movement_path)
        
        # 1. get_chunksize()
        print(movement_extractor.get_chunksize())