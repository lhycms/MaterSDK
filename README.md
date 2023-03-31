# pflow

A python library for computational materials science. You can use it to:
1. Process structure file (For PWmat, VASP, ATAT, ...)
2. Process Kohn-Sham state (For PWmat, VASP, ...)
3. Process periodic boundary condition

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


## 1.3. Process perioidic boundary condination
1. Analyse the local environment for arbitrary crystalline materials, whatever the dimension of material.
2. Generate DEEPMD features conviently. Don't worry about the choice of cutoff radius.
3. You can attach any feature to sites in primitive cell.
4. ...