import unittest

# python3 -m pflow.io.pwmat.output.test.test_report
from ..report import Report


class ReportTest(unittest.TestCase):
    def test_report_all(self):
        report_path = "/data/home/liuhanyu/hyliu/pwmat_demo/band_ispin/REPORT"
        out_fermi_path = "/data/home/liuhanyu/hyliu/pwmat_demo/band_ispin/OUT.FERMI"
        report = Report(report_path=report_path)
        
        print("\n1. 能带数:", end="\t")
        print( report.get_num_bands() )
        
        print("\n2. kpoints的数目:", end="\t")
        print(report.get_num_kpts())
        
        print("\n3. 得到所有kpoints的本征能量:")
        print(report.get_eigen_energies())
        
        print("\n4. IN.ATOM: ", end="\t")
        print(report.get_in_atom())
        
        print("\n5. self._is_metal: ", end="\t")
        print(report._is_metal(out_fermi_path=out_fermi_path))
        
if __name__ == "__main__":
    unittest.main()