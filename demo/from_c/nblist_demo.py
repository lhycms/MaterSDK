from pymatgen.core import Structure
import numpy as np
import time


poscar_path = "/data/home/liuhanyu/hyliu/pwmat_demo/MoS2/scf_/POSCAR"
structure = Structure.from_file(poscar_path).make_supercell([6, 6, 1])
lattice = structure.lattice.matrix
species = np.array([tmp_specie.Z for tmp_specie in structure.species])
frac_coords = structure.frac_coords
rcut = 3.2
pbc_xyz = [True, True, False]
umax_num_neigh_atoms = 19
sort = True


import sys
sys.path.append("/data/home/liuhanyu/hyliu/code/matersdk/source/matersdk/io/publicLayer/bind/gen")
import nblist

start = time.time()
for _ in range(5000):
    #print(_)
    info = nblist.find_info4mlff(
        lattice,
        species.astype(np.int32),
        frac_coords,
        rcut,
        pbc_xyz,
        umax_num_neigh_atoms,
        sort
    )
end = time.time()
print(end-start)

print(info[3])


"""
for ii in range(12):
    for jj in range(19):
        tmp_rc = info[4][ii][jj]
        print(np.linalg.norm(tmp_rc), end=", ")
    print("\n")
"""