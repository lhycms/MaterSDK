from matersdk.io.publicLayer.structure import DStructure
from matersdk.io.publicLayer.neigh import StructureNeighborsDescriptor

from matersdk.io.pwmat.output.movement import Movement
from matersdk.data.deepmd.data_system import DpLabeledSystem


### Step 0. 自定义参数
atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
rcut = 6.0                  # Rmax
rcut_smooth = 5.5           # Rmin
scaling_matrix = [3, 3, 3]  # 三维体系：[3, 3, 3]; 二维体系: [3, 3, 1]
reformat_mark = True        # 永远都是True
coords_are_cartesian = True # 永远都是True



### Step 1. 由 atom.config 获取结构
structure = DStructure.from_file(file_path=atom_config_path, file_format="pwmat")

### Step 2. 获取 model.forward() 需要的参数
### Step 2.1. ImageDR: 相对中心原子的距离 -- shape = (1, num_centers, max_num_nbrs, 3)
neigh_list = StructureNeighborsDescriptor.create(
                    'v1',
                    structure=structure,
                    rcut=rcut,
                    scaling_matrix=scaling_matrix,
                    reformat_mark=reformat_mark,
                    coords_are_cartesian=coords_are_cartesian
)



### Step 2.2. Ri: Rij in DeepPot-SE
movement = Movement()