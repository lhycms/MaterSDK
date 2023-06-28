import numpy as np
from matersdk.io.publicLayer.structure import DStructure
from matersdk.io.publicLayer.neigh import StructureNeighborsDescriptor

from matersdk.io.pwmat.output.movement import Movement
from matersdk.data.deepmd.data_system import DpLabeledSystem
from matersdk.feature.deepmd.preprocess import TildeRNormalizer

from matersdk.io.publicLayer.neigh import StructureNeighborsUtils

### Step 0. 自定义参数
atom_config_path = "/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config"
rcut = 6.0                  # Rmax
rcut_smooth = 5.5           # Rmin
scaling_matrix = [3, 3, 3]  # 三维体系：[3, 3, 3]; 二维体系: [3, 3, 1]
reformat_mark = True        # 永远都是True
coords_are_cartesian = True # 永远都是True
davg = None                 # 设置为 None 时，可以从 MOVEMENT 自动计算
dstd = None                 # 设置为 None 时，可以从 MOVEMENT 自动计算
movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"  # 产生训练集、验证集的 MOVEMENT
center_atomic_numbers = [3, 14] # 体系内所有元素的原子序数，从小到大排列
nbr_atomic_numbers = [3, 14] # 体系内所有元素的原子序数，从小到大排列
max_num_nbrs = [100, 100]    # 近邻原子的最大数目，与 nbr_atomic_numbers 对应




### Step 1. 由 atom.config 获取结构
movement = Movement(movement_path=movement_path)
structure = movement.get_frame_structure(idx_frame=500)

### Step 2. 获取 model.forward() 需要的参数
### Step 2.1. ImageDR: 相对中心原子的距离 -- shape = (1, num_centers, max_num_nbrs, 3)
struct_nbr = StructureNeighborsDescriptor.create(
                    'v1',
                    structure=structure,
                    rcut=rcut,
                    scaling_matrix=scaling_matrix,
                    reformat_mark=reformat_mark,
                    coords_are_cartesian=coords_are_cartesian
)



### Step 2.2. Ri: Rij in DeepPot-SE ; Ri_d: The detivative of Rij
if (not davg) or (not dstd):
    movement = Movement(movement_path=movement_path)
    dpsys = DpLabeledSystem.from_trajectory_s(trajectory_object=movement)
    tildeR_normalizer = TildeRNormalizer.from_dp_labeled_system(
                            dp_labeled_system=dpsys,
                            structure_indices=[*range(10)],
                            rcut=rcut,
                            rcut_smooth=rcut_smooth,
                            center_atomic_numbers=center_atomic_numbers,
                            nbr_atomic_numbers=nbr_atomic_numbers,
                            max_num_nbrs=max_num_nbrs,
                            scaling_matrix=scaling_matrix
    )
    # shape = (num_centers, 4)
    davg, dstd = tildeR_normalizer.davgs, tildeR_normalizer.dstds
else:
    tildeR_normalizer = TildeRNormalizer(
                            rcut=rcut,
                            rcut_smooth=rcut_smooth,
                            center_atomic_numbers=center_atomic_numbers,
                            nbr_atomic_numbers=nbr_atomic_numbers,
                            max_num_nbrs=max_num_nbrs,
                            scaling_matrix=scaling_matrix,
                            davgs=davg,
                            dstds=dstd
    )

Ri, Ri_d = tildeR_normalizer.normalize(structure=structure)
# shape = (1, num_centers, max_num_nbrs, 4)
Ri = np.expand_dims(Ri, axis=0)
# shape = (1, num_centers, max_num_nbrs, 4, 3)
Ri_d = np.expand_dims(Ri_d, axis=0)
print("Step 2.1. Ri.shape = ", end='\t')
print(Ri.shape)
print("Step 2.2. Ri_d.shape = ", end='\t')
print(Ri_d.shape)


### Step 3. neigh_list
neigh_list_ = StructureNeighborsUtils.get_nbrs_indices(
                            struct_nbr=struct_nbr,
                            center_atomic_numbers=center_atomic_numbers,
                            nbr_atomic_numbers=nbr_atomic_numbers,
                            max_num_nbrs=max_num_nbrs
)
neigh_list = np.expand_dims(neigh_list_, axis=0)
print("Step 3. neigh_list.shape = ", end='\t')
print(neigh_list.shape)


### Step 4. e.g. 72原子的Li2Si -- natoms = [72, 48, 24]
natoms = structure.get_natoms()
natoms_img = np.repeat(natoms[np.newaxis, :], 1, axis=0)
print("Step 4. natoms_img.shape = ", end='\t');
print(natoms_img.shape)
