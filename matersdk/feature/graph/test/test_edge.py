import unittest


# python3 -m matersdk.feature.graph.test.test_edge
from ....io.publicLayer.structure import DStructure
from ..edge import AdjacentMatrix



class AdjacentMatrixTest(unittest.TestCase):
    def test_all(self):
        atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
        scaling_matrix = [3, 3, 3]
        structure = DStructure.from_file(
                            file_format="pwmat",
                            file_path=atom_config_path
                            )
        rcut = 3.2
        
        
        adjacent_matrix = AdjacentMatrix(
                                structure=structure,
                                rcut=rcut,
                                scaling_matrix=scaling_matrix
                                )
        
        ### Step 1. 
        print()
        print("Step 1. get_neigh_primitive_frac_coords:")
        #adjacent_matrix._get_neigh_primtive_cell_frac_coords()
        
        ### Step 2. 
        print()
        print("Step 2. The adjacent matrix (radius cutoff = {0})".format(rcut))
        adjacent_matrix.get_adjacent_matrix()


if __name__ == "__main__":
    unittest.main()