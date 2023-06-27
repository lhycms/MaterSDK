import numpy as np
from matersdk.io.pwmat.output.movement import Movement
from matersdk.feature.avg.msd import Msd


### Step 1. 自定义参数
# 1. MOVEMENT 文件路径
movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
# 2. msd 文件保存路径
msd_file_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/msd.dat"


### Step 2. 对象初始化
movement = Movement(movement_path=movement_path)
msd_object = Msd(trajectory=movement)

### Step 3. 计算对应步数的msd
msd_values_lst = msd_object.calc_msd()
steps_lst = [*range(len(msd_values_lst))]

steps_array = np.array(steps_lst).reshape(-1, 1)
msd_values_array = np.array(msd_values_lst).reshape(-1, 1)
tot_array = np.concatenate([steps_array, msd_values_array], axis=1)

### Step 4. 保存文件
### Step 4.1. 保存 msd 文件
np.savetxt(fname=msd_file_path, X=tot_array)
