# Material Software Development Kit (MaterSDK)
A python library for computational materials science. You can use it to:
1. Process structure file (For PWmat, VASP, ATAT, ...)
2. Process Kohn-Sham state (For PWmat)
3. Generate feature for structure.

## Contributors

- [Liu Hanyu](https://github.com/lhycms)
- [Tian Hongzhen](https://github.com/LonxunQuantum)

# 1. Brief introduction
## 1.1. Process structure file
You can:
1. Convert file format arbitrarily as below diagram
2. Choose high symmetry for 2D/3D material
3. Choose KMesh
4. ...

<img src="./demo/pics/pic_1.png" width = "500" height = "300"  />

## 1.2. Process Kohn-Sham state
1. Support several types of ab initio software, particularlly PWmat.
2. Plot bandstructure, density of state directly.
3. Process something in reciprocal space.
4. ...

<img src="./demo/pics/图片1.png" width = "300" height = "240"  />       <img src="./demo/pics/图片2.png" width = "300" height = "240"  />


## 1.3. Pick out specific frame from `MOVEMENT`
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/movement/movement.ipynb;
2. Extract information from specific frame:
   1. Extract a `DStructure` object;
      1. `atomic_energy`;
      2. `atomic_force`;
      3. `atomic_velocity`;
      4. `magmom`;
   2. Extract `virial tensor` of specific frame;
   3. Extract `volume` of specific frame;
   4. Extract `energy` of specific frame;
3. Calculate the `max number of neighbors (max_num_nbrs_real)` and `max_num_nbrs_real_element` among all `DStructrue` in `Movement`.
4. `Movement.get_all_frame_structures_info()`:
   1. `structures_lst:List[DStructure]`
   2. `total_energys_array`
   3. `potential_energys_array`
   4. `kintic_energys_array`
   5. `virial_tensors_array`


## 1.4. Generate feature for structure. (For DeepPot-SE, GNN and so on)
Using following algorithm to handle the periodic boundary conditions (Take an 2D materials as an demo). Our new `expansion function (DStructure.make_supercell_())` requires:
1. The corresponding surfaces of structure to be parallel.
2. The element in `scaling matrix` must be odd number.

### 1.4.0. `matersdk.io.publicLayer.neigh.StructureNeighborUtils`
1. Get `max_num_nbrs_real`: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/neighs/SturtureNeighbors.ipynb
2. Get `max_num_nbrs_real_element`: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/neighs/SturtureNeighbors.ipynb
3. Mapping the site in supercell to site in primitive_cell, and get `neigh_list`: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/dpseTileRPair.ipynb

### 1.4.1. `matersdk.io.publicLayer.neigh.StructureNeighborsBase`: Get information of neighbors atom for center atom.

#### 1.4.1.1. Use `matersdk.io.publicLayer.neigh.StructureNeighborsV1` to analyse the surrouing environment for specific atom.
1. Note: If not set `max_num_nbrs`, `max_num_nbrs = max_num_nbrs_real + 1` !!!
2. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/neighs/SturtureNeighbors.ipynb

### 1.4.2. `matersdk.feature.deepmd.premise.DpFeaturePairPremiseBase`: Get pair information of neighbor atoms for center atom

#### 1.4.2.1. `matersdk.feature.deepmd.premise.DpFeaturePairPremiseV1`
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/premise.ipynb

### 1.4.3. Smooth edition $\widetilde{R}$ of Deepmd feature pair -- `matersdk.feature.deepmd.se_pair.DpseTildeRPairBase`
#### 1.4.3.1. `matersdk.feature.deepmd.se.DpseTildeRPairV1`
1. 2D system: https://github.com/lhycms/MaterSDK/blob/reduce/demo/feature/deepmd/dpseTileRPair.ipynb
2. 3D system: https://github.com/lhycms/MaterSDK/blob/reduce/demo/feature/deepmd/dpseTileRPair_3d.ipynb

#### 1.4.3.2. `matersdk.feature.deepmd.se.DpseTildeRV1`
1. 2D system: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/dpseTildeR.ipynb


### 1.4.4. `NormalizerPremise` + `TildeRPairNormalizer`: normalize the environment matrix($\tilde{R}^i$)
#### 1.4.4.1. Normalize the $\tilde{R}$ for the specific center atoms (e.g. Li-Li&Si)
1. https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/normalizer_for_single_pair.ipynb

#### 1.4.4.2. Calculate the `davgs` and `dstds` using the first 10 frames.
1. Demo: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/normalizer.ipynb

#### 1.4.4.3. `TildeRNormalizer.normalize(DStructure, ...)`
1. Demo: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/normalizer.ipynb

#### 1.4.4.4. `TildeRNormalizer.to()` & `TildeRNormalizer.from_file()`
1. Demo: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/normalizer.ipynb


## 1.5. Adjacent Matrix
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/graph/adjacent_matrix.ipynb


## 1.6. `DpLabeledSystem`
### 1.6.1. Use `DpLabeledSystem` to store some neighbor information in Image folder
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/data_system.ipynb

### 1.6.2. Use indices to get trainset and testset
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/data_system.ipynb

### 1.6.3. Almost same way of storage as Deepmd-kit
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/data_system_dp.ipynb


# 2. Installation
## 2.1. Online
```shell
$ git clone git@github.com:lhycms/MaterSDK.git
$ cd matersdk
$ pip install .
```

## 2.1. Offline
1. You can download a python interpreter containing `matersdk` from https://www.jianguoyun.com/p/DfhQFx8Q_qS-CxifgfwEIAA.matersdk.egg-info


# 3. Integrate MaterSDK with `lammps` as third-party libraries
```shell
# Step 1.
$ mkdir <your_path_lammps>/lammps/src/MATERSDK
$ cd <your_path_lammps>/lammps/src/MATERSDK

# Step 2.
$ cp <your_path_to_matersdk>/matersdk/source ./.
$ cp <your_path_to_matersdk>/matersdk/cmake ./.
$ ls .
source cmake LICENSE
$ vim <your_path_to_lammps>/src/MATERSDK/source/CMakeLiest.txt
### Modify the path of MATERSDK
add_subdirectory("<your_path_to_lammps>/src/MATERSDK/source/ext/googletest-main")
$ cd vim <your_path_to_lammps>/src/MATERSDK/source/build; cmake ..; make

# Step 3. 
$ vim <your_path_to_lammps>/cmake/CMakeLists.txt
### Add content below to CMakeLists.txt
include("<your_path_to_lammps>/src/MATERSDK/cmake/matersdk.cmake")   ### matersdk
include_directories(${MATESDK_INCLUDE_DIRS}) ### matesdk

# Link libraries of matersdk to lammps.so/lammps.a
add_library(lammps ${ALL_SOURCES})  ### matersdk
target_link_libraries(lammps PUBLIC ${MATERSDK_LIBRARIES})  ### matersdk

# Step 4. 
$ cd <your_path_to_lammps>/build
$ cmake ../cmake 
-- +++ MATERSDK_DIR : /data/home/liuhanyu/hyliu/code1/lammps_gpu/src/MATERSDK
-- +++ MATERSDK_INCLUDE_DIRS : /data/home/liuhanyu/hyliu/code1/lammps_gpu/src/MATERSDK/source/include;/data/home/liuhanyu/hyliu/code1/lammps_gpu/src/MATERSDK/source/core/include;/data/home/liuhanyu/hyliu/code1/lammps_gpu/src/MATERSDK/source/matersdk/io/publicLayer/include;/data/home/liuhanyu/hyliu/code1/lammps_gpu/src/MATERSDK/source/matersdk/feature/deepmd/include
-- +++ MATERSDK_LIBRARIES : /data/home/liuhanyu/hyliu/code1/lammps_gpu/src/MATERSDK/source/build/lib/matersdk/feature/deepmd/se4pw_op
```