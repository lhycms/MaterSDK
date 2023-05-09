from timeit import default_timer as timer
from matersdk.io.pwmat.output.movement import Movement
from matersdk.data.deepmd.data_system import DeepmdDataSystem


### Parameters
movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo2/PWdata/data1/MOVEMENT"
rcut = 6.5
output_dir = "/data/home/liuhanyu/hyliu/code/test"

start = timer()
movement = Movement(movement_path=movement_path)
dp_data_system = DeepmdDataSystem.from_trajectory_s(
                                trajectory_object=movement,
                                rcut=rcut)
dir_path = output_dir

### Step 1.
dp_data_system.save(dir_path=dir_path)
end = timer()

print("Costing time: {0} s".format(end-start))