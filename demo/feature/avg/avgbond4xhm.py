import numpy as np
from scipy.fft import fft
from matersdk.feature.avg.avgbond import AvgBond
from matersdk.io.pwmat.output.movement import Movement
from matersdk.data.deepmd.data_system import DpLabeledSystem


### Part I. Custom parameters
movement_path = "/data/home/liuhanyu/hyliu/pwmat_demo/xhm/MOVEMENT"
element_1 = "Ge"
element_2 = "Te"

rcut = 3.2    # 在此范围内认为 `element_1` 与 `element_1` 成键

save_path = "./data.txt"


### Part II. 查看 MOVEMENT 信息
movement = Movement(movement_path=movement_path)
dsys = DpLabeledSystem.from_trajectory_s(movement, rcut=rcut)
print(dsys)


### Part III. 运行程序
avgbond = AvgBond(
            movement_path=movement_path,
            element_1=element_1,
            element_2=element_2,
            rcut=rcut)
print("\nCalculating the avg bond length for all frames in MOVEMENT...")
frame_avg_bonds_lst = avgbond.get_frames_avg_bond()
xs_array = np.array([*range(len(dsys))]).reshape(-1, 1)
ys_array = np.array(frame_avg_bonds_lst).reshape(-1, 1)
xys_array = np.concatenate([xs_array, ys_array], axis=1)
np.savetxt(fname=save_path, X=xys_array)
