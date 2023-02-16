import unittest

# python3 -m pflow.io.pwmat.input.test.test_inkpt
from ..inkpt import Inkpt


class InkptTest(unittest.TestCase):
    def test_get_iflag(self):
        in_kpt_path = "/data/home/liuhanyu/hyliu/pwmat_demo/band/IN.KPT"
        atom_config_path = "/data/home/liuhanyu/hyliu/pwmat_demo/band/atom.config"
      
    
        in_kpt_path = "/data/home/liuhanyu/hyliu/pwmat_demo/scf_3d/IN.KPT"
        atom_config_path = "/data/home/liuhanyu/hyliu/pwmat_demo/scf_3d/atom.config"
        
          
        in_kpt = Inkpt(in_kpt_path=in_kpt_path)
        print("\n1. The `iflag` of IN.KPT:", end="\t")
        print(in_kpt._get_iflag())
        
        print("\n2. The `a0` of IN.KPT:", end="\t")
        print(in_kpt._get_a0())
        
        print("\n3. The number of kpoints:", end="\t")
        print(in_kpt.get_num_kpts())
        
        print("\n4. The fractional coordinations of all kpoints:", end="\n")
        print(in_kpt.get_kpt_coords_frac())

        print("\n5. The coordinations of all kpoints (unit: 埃):")
        print(in_kpt.get_kpt_coords_A(atom_config_path=atom_config_path))
        
        print("\n5. The coordinations of all kpoints (unit: Bohr):")
        print(in_kpt.get_kpt_coords_Bohr(atom_config_path=atom_config_path))

        print("\n6. The weights of kpoints:")
        print(in_kpt.get_kpt_weights())
        
        print("\n7. The High symmetry points for ")
        print(in_kpt.get_hsp())
        
        print("\n8. The distance from gamma (unit: 埃)")
        print(in_kpt.get_distance_from_gamma_A())


if __name__ == "__main__":
    unittest.main()