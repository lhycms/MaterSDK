#ifndef MATERSDK_ENVMATRIX_H
#define MATERSDK_ENVMATRIX_H
#include <cstring>
#include "./se.h"

namespace matersdk {
namespace deepPotSE {

template <typename CoordType>
class EnvMatrix {
public:
    static void find_value_deriv(
        CoordType* tilde_r,
        CoordType* tilde_r_deriv,
        int inum, 
        int* ilist,
        int* numneigh,
        CoordType* firstneigh,
        CoordType* relative_coords,
        int* types,
        int ntypes,
        int* umax_num_neigh_atoms_lst);
};  // class : EnvMatrix


/**
 * @brief Find the value and deriv wrt. R_{ij} of Environment Matrix
 * 
 * @tparam CoordType 
 * @param tilde_r : value of Environment Matrix.
 *          len = inum * sum(umax_num_neigh_atoms_lst) * 4
 * @param tilde_r_deriv : deriv wrt. R{ij} of Environment Matrix
 *          len = inum * sum(umax_num_neigh_atoms_lst) * 4 * 3
 * @param inum : number of primitive cell atoms
 * @param ilist : len = inum
 * @param numneigh : len = inum
 * @param firstneigh : len = inum * umax_num_neigh_atoms
 * @param relative_coords : len = inum * nmax_num_neigh_atoms * 3
 * @param types : len = inum, depends on `firstneigh`
 * @param ntypes : number of element kinds
 * @param umax_num_neigh_atoms_lst : len = ntypes. User specified `umax_num_neigh_atoms` for all elements.
 * 
 */
template <typename CoordType>
void EnvMatrix<CoordType>::find_value_deriv(
        CoordType* tilde_r,
        CoordType* tilde_r_deriv,
        int inum, 
        int* ilist,
        int* numneigh,
        CoordType* firstneigh,
        CoordType* relative_coords,
        int* types,
        int ntypes,
        int* umax_num_neigh_atoms_lst)
{   
    // Step 1. 
    // Step 1.1. Init tilde_r, tilde_r_deriv
    int umax_num_neigh_atoms = 0;
    for (int ii=0; ii<ntypes; ii++)
        umax_num_neigh_atoms += umax_num_neigh_atoms_lst[ii];
    memset(tilde_r, 0.0, sizeof(CoordType) * inum * umax_num_neigh_atoms * 4);
    memset(tilde_r_deriv, 0.0, sizeof(CoordType) * inum * umax_num_neigh_atoms * 4 * 3);
    
    // Step 1.2. 
    int* nstart_idxs = (int*)malloc(sizeof(int) * ntypes);
    memset(nstart_idxs, 0, sizeof(int) * ntypes);
    for (int ii=0; ii<ntypes; ii++)
        for (int jj=0; jj<ii; jj++)
            nstart_idxs[ii] += umax_num_neigh_atoms_lst[jj];
    int* nloop_idxs = (int*)malloc(sizeof(int) * ntypes);
    memset(nloop_idxs, 0, sizeof(int) * ntypes);


    // Step 2. 
    for (int ii=0; ii<inum; ii++) {
        for (int jj=0; jj<numneigh[ii]; jj++) {
            int tmp_neigh_idx = firstneigh[ii][jj];
            int tmp_neigh_type = types[tmp_neigh_idx];
            CoordType distance = relative_coords[ii][jj];
            
        }
    }


    // Step . Free memory
    free(nstart_idxs);
    free(nloop_idxs);
}


}   // namespace : deepPotSE
}   // namespace : matersdk


#endif