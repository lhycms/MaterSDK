import unittest


# python3 -m boNEB.io.utilitys.test.test_lineLocator
from ..lineLocator import LineLocator


class LineLocatorTest(unittest.TestCase):
    def test_locate_all_lines(self):
        movement_path = "/Users/mac/我的文件/Mycode/new/new2/boNEB/test_data/movement/MOVEMENT"
        row_idxs_lst = LineLocator.locate_all_lines(file_path=movement_path,
                                                    content="Position")
        print( len(row_idxs_lst) )



if __name__ == "__main__":
    unittest.main()