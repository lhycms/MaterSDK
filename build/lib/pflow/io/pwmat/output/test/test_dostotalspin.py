import unittest

# python3 -m pflow.io.pwmat.output.test.test_dostotalspin
from ..dostotalspin import Dostotalspin

class DostotalspinTest(unittest.TestCase):
    def test_dostotalspin(self):
        dos_totalspin_path = "/data/home/liuhanyu/hyliu/pwmat_demo/dos_ispin/DOS.spinup"
        dos_totalspin = Dostotalspin(dos_totalspin_path=dos_totalspin_path)
        
        print("\n1. tdos:")
        print(dos_totalspin.get_tdos())


if __name__ == "__main__":
    unittest.main()