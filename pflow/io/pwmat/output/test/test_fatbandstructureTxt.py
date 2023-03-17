import unittest

# python3 -m pflow.io.pwmat.output.test.test_fatbandstructureTxt
from ..fatabandstructureTxt import FatbandStructure


class FatbandStructureTest(unittest.TestCase):
    def test_fatbandstructure(self):
        fatbandstructure_txt_path = "/data/home/liuhanyu/hyliu/pwmat_demo/MoS2/scf/nonscf/dos/fatbandstructure_1.txt"
        fatbandstructure = FatbandStructure(
                fatbandstructure_txt_path=fatbandstructure_txt_path,
        )        
        
        ### Step 1. 得到能带的条数
        print("\n1. 能带的总数目:", end="\t")
        print(fatbandstructure._get_num_bands())
        
        ### Step 2. 得到 kpoints 的数目
        print("\n2. KPOINTS的总数目:", end="\t")
        print(fatbandstructure._get_num_kpoints())
        
        ### Step 3. 得到 `BAND` 所在行的索引
        print("\n3. BAND 所在行的索引:")
        print(fatbandstructure._get_BAND_mark_idxs())

        ### Step 4. preproceee 读取 pd.DataFrame，注意不包括空行的`BAND行`
        print("\n4. 预处理后的 DataFrame:")
        #print(fatbandstructure._preprocess())
        
        
        ### Step 5. preproceee 读取 pd.DataFrame，注意不包括空行的`BAND行`
        print("\n5. dfs_lst:")
        print(fatbandstructure.get_dfs()[0])


if __name__ == "__main__":
    unittest.main()