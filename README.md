# Material Software Development Kit (MaterSDK)
A python library for computational materials science. You can use it to:
1. Process structure file (For PWmat, VASP, ATAT, ...)
2. Process Kohn-Sham state (For PWmat)
3. Generate feature for structure.

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
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/movement/movement.ipynb
2. Extract a `DStructure` object
   1. `atomic_energy`
   2. `atomic_force`
   3. `atomic_velocity`
   4. `magmom`
3. Extract `virial tensor` of specific frame
4. Extract `volume` of specific frame
5. Extract `energy` of specific frame


## 1.4. Generate feature for structure. (For DeepPot-SE, GNN and so on)
Using following algorithm to handle the periodic boundary conditions (Take an 2D materials as an demo). Our new `expansion function (DStructure.make_supercell_())` requires:
1. The corresponding surfaces of structure to be parallel.
2. The element in `scaling matrix` must be odd number.


### 1.4.1. Use `matersdk.io.publicLayer.neigh.StructureNeighborsV1` to analyse the surrouing environment for specific atom.
1. `matersdk.io.publicLayer.neigh.StructureNeighborsV2` is much faster than `matersdk.io.publicLayer.neigh.StructureNeighbors`.
2. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/neighs/SturtureNeighbors.ipynb

### 1.4.2. Use `matersdk.io.publicLayer.neigh.StructureNeighborsV2` to analyse the surrouing environment for specific atom.
1. `matersdk.io.publicLayer.neigh.StructureNeighborsV2` is much faster than `matersdk.io.publicLayer.neigh.StructureNeighbors`.
2. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/neighs/SturtureNeighbors.ipynb

### 1.4.3. Use `matersdk.io.publicLayer.neigh.StructureNeighborsV3` to analyse the surrouing environment for specific atom.
1. `matersdk.io.publicLayer.neigh.StructureNeighborsV3` is much faster than `matersdk.io.publicLayer.neigh.StructureNeighborsV1`. And save more space than `StrcutureNeighborsV2`
2. Demo url: 

### 1.4.4. Premise for Deepmd feature pair
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/premise.ipynb

### 1.4.5. Smooth edition $\widetilde{R}$ of Deepmd feature pair 
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/deepmdTileR_3d.ipynb

## 1.5. Adjacent Matrix
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/graph/adjacent_matrix.ipynb


## 1.6. Use `DeepmdDataSystem` to store some neighbor information in Image folder
1. Demo url: https://github.com/lhycms/MaterSDK/blob/main/demo/feature/deepmd/data_system.ipynb

# 2. Installation
## 2.1. Online
```shell
$ git clone git@github.com:lhycms/MaterSDK.git
$ cd matersdk
$ pip install .
```

## 2.1. Offline
1. You can download a python interpreter containing `matersdk` from https://www.jianguoyun.com/p/DfhQFx8Q_qS-CxifgfwEIAA.