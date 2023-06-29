import torch
import numpy as np
from matersdk.io.pwmat.output.movement import Movement
from matersdk.data.deepmd.data_system import DpLabeledSystem
from matersdk.feature.deepmd.preprocess import TildeRNormalizer
from matersdk.infer.pwmatmlff.deepmd.preprocess import InferPreprocessor



movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"  # 产生训练集、验证集的 MOVEMENT
rcut = 6.0                      # Rmax
rcut_smooth = 5.5               # Rmin
scaling_matrix = [3, 3, 3]      # 三维体系：[3, 3, 3]; 二维体系: [3, 3, 1]
davg = None                     # 设置为 None 时，可以从 MOVEMENT 自动计算
dstd = None                     # 设置为 None 时，可以从 MOVEMENT 自动计算
center_atomic_numbers = [3, 14] # 体系内所有元素的原子序数，从小到大排列
nbr_atomic_numbers = [3, 14] # 体系内所有元素的原子序数，从小到大排列
max_num_nbrs = [100, 100]    # 近邻原子的最大数目，与 nbr_atomic_numbers 对应
reformat_mark = True            # 永远都是True
coords_are_cartesian = True # 永远都是True

movement = Movement(movement_path=movement_path)
### 计算 Rij 的 davg, dstd
if (not davg) or (not dstd):
    dpsys = DpLabeledSystem.from_trajectory_s(trajectory_object=movement)
    tildeR_normalizer = TildeRNormalizer.from_dp_labeled_system(
                        dp_labeled_system=dpsys,
                        structure_indices=[*(range(10))],
                        rcut=rcut,
                        rcut_smooth=rcut_smooth,
                        center_atomic_numbers=center_atomic_numbers,
                        nbr_atomic_numbers=nbr_atomic_numbers,
                        max_num_nbrs=max_num_nbrs,
                        scaling_matrix=scaling_matrix
    )
    davg, dstd = tildeR_normalizer.davgs, tildeR_normalizer.dstds

infer_preprocessor = InferPreprocessor(
            structure=movement.get_frame_structure(idx_frame=100),
            rcut=rcut,
            rcut_smooth=rcut_smooth,
            scaling_matrix=scaling_matrix,
            davg=davg,
            dstd=dstd,
            center_atomic_numbers=center_atomic_numbers,
            nbr_atomic_numbers=nbr_atomic_numbers,
            max_num_nbrs=max_num_nbrs,
            reformat_mark=reformat_mark,
            coords_are_cartesian=coords_are_cartesian
)

device = "cuda:0"

### Step 1. 
ImageDR = infer_preprocessor.expand_rc()
print("1. ImageDR.shape = ", ImageDR.shape)
ImageDR_tensor = torch.from_numpy(ImageDR).double().to(device).requires_grad_()

### Step 2.
Ri, Ri_d = infer_preprocessor.expand_tildeR()
print("2. Ri.shape = ", Ri.shape)
print("3. Ri_d.shape =", Ri_d.shape)
Ri_tensor = torch.from_numpy(Ri).double().to(device).requires_grad_()
Ri_d_tensor = torch.from_numpy(Ri_d).double().to(device).requires_grad_()

### Step 3.
list_neigh = infer_preprocessor.expand_list_neigh()
print("4. list_neigh.shape = ", list_neigh.shape)
list_neigh_tensor = torch.from_numpy(list_neigh).double().to(device).requires_grad_()

### Step 4. 
natoms_image = infer_preprocessor.expand_natoms_img()
print("5. natoms_img = ", natoms_image.shape)
#natoms_image_tensor = torch.from_numpy(natoms_img).float()

### Step 5. 
#model = torch.load(f="/data/home/liuhanyu/hyliu/code/mlff/PWmatMLFF_dev/test/demo2/record/best.pth.tar", map_location="cpu")
model = torch.load(f="/data/home/liuhanyu/hyliu/code/mlff/PWmatMLFF_dev/test/demo2/record/checkpoint.pt")
model.to(device)
model.eval()


result = model(
            ImageDR_tensor, 
            Ri_tensor, 
            Ri_d_tensor, 
            list_neigh_tensor,
            natoms_image
)
print("\n\n")
print("Inference Result:")
print("\tStep 1. Etot = " , result[0].item())
print("\tStep 2. Ei.shape = ", result[1].shape)
print("\tStep 3. F.shape = ", result[2].shape)
print("\tStep 4. Virial.shape = \n", result[3].view((3, 3)))
