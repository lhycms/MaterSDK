%module nblist

%{
#include "./nblist_bind.h"
#include "../include/neighborList.h"
%}

PyObject* find_info4mlff(
        PyObject* lattice_py,
        PyObject* atomic_numbers_py,
        PyObject* frac_coords_py,
        PyObject* rcut_py,
        PyObject* pbc_xyz_py,
        PyObject* umax_num_neigh_atoms,
        PyObject* sort_py);