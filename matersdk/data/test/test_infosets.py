import unittest

# python3 -m matersdk.data.test.test_infosets
from ..infosets import InfoSets


class InfoSetsTest(unittest.TestCase):
    def test_all(self):
        dir_path = "/data/home/liuhanyu/hyliu/code/test/data"
        file_name = "MOVEMENT"
        file_format = "pwmat/movement"
        
        #info_sets = InfoSets(dir_path=dir_path, file_name=file_name, file_format=file_format)
        
        info_sets = InfoSets.from_dir(dir_path=dir_path, file_name=file_name, file_format=file_format)
        print(info_sets.get_num_frames())
        info_sets.to_dir(dir_path=dir_path)
    
if __name__ == "__main__":
    unittest.main()