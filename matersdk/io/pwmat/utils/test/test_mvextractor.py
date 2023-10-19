import unittest

from ..mvextractor import MVExtractor


# python3 -m matersdk.io.pwmat.utils.test.test_mvextractor
class MVExtractorTest(unittest.TestCase):
    def test_all(self):
        #movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo1/PWdata/data1/MOVEMENT"
        #movement_path = "/data/home/liuhanyu/hyliu/code/mlff/PWmatMLFF_dev/test/SiC/MD/MOVEMENT"
        movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
        virial_mark = True
        magmoms_mark = False
        eatoms_mark = False
        
        
        mvextractor = MVExtractor(movement_path=movement_path, virial=virial_mark, magmoms=magmoms_mark, eatoms=eatoms_mark)
        #print(mvextractor.get_chunksizes())
        #print(mvextractor.get_chunkslices())
        #print(mvextractor.get_frame_str(fidx=0))
        box, types, coords, etot, fatoms, virial = mvextractor.get_frame_info(fidx=0)
        print(virial)


if __name__ == "__main__":
    unittest.main()
    