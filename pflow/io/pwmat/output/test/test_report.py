import unittest

# python3 -m pflow.io.pwmat.output.test.test_report
from ..report import Report


class ReportTest(unittest.TestCase):
    def test_report_all(self):
        report_path = "/data/home/liuhanyu/hyliu/pwmat_demo/band/REPORT"
        report = Report(report_path=report_path)
        
        print("\n1. 能带数:", end="\t")
        print( report.get_num_bands() )
        
        print("\n2. kpoints的数目:", end="\t")
        print(report.get_num_kpts())
        
        print("\n3. 得到所有kpoints的本征能量:")
        print(report.get_eigen_energies())
        
        
        
if __name__ == "__main__":
    unittest.main()