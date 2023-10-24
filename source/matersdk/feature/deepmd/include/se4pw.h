#ifndef MATERSDK_SE4PW_H
#define MATERSDK_SE4PW_H
#include "./se.h"

namespace matersdk {
namespace deepPotSE {


template <typename CoordType>
class Se4pw {
public:
    static void generate(
        CoordType* tilde_r,
        int inum,
        int* ilist,
        int* numneigh,
        int** firstneigh,
        CoordType** x,
        int* types, // starts from 0.
        int ntypes, // starts from 0. e.g. 2
        int* num_neigh_atoms_lst,
        CoordType rcut,
        CoordType rcut_smooth);
};


template <typename CoordType>
void Se4pw<CoordType>::generate(
        CoordType* tilde_r,
        int inum,
        int* ilist,
        int* numneigh,
        int** firstneigh,
        CoordType** x,
        int* types, // starts from 0.
        int ntypes, // starts from 0. e.g. 2
        int* num_neigh_atoms_lst,
        CoordType rcut,
        CoordType rcut_smooth)
{
    // Step 1.
    // Step 1.1. $\tilde{R}$ = (s(r_{ji}), x_{ji}, y_{ji}, z_{ji})
    //  = (tilde_s_value, tilde_x_value, tilde_y_value, tilde_z_value)
    CoordType tilde_s_value;
    CoordType tilde_x_value;
    CoordType tilde_y_value;
    CoordType tilde_z_value;
    CoordType distance_ji;
    CoordType distance_ji_recip;
    int center_atom_idx;
    int neigh_atom_idx;
    CoordType* center_cart_coords = (CoordType*)malloc(sizeof(CoordType) * 3);
    CoordType* neigh_cart_coords = (CoordType*)malloc(sizeof(CoordType) * 3);
    CoordType* diff_cart_coords = (CoordType*)malloc(sizeof(CoordType) * 3);

    int* nstart_idxs = (int*)malloc(sizeof(int) * ntypes);
    for (int ii=0; ii<ntypes; ii++)
        nstart_idxs[ii] = 0;
    for (int ii=0; ii<ntypes; ii++)
        for (int jj=0; jj<ii; jj++)
            nstart_idxs[ii] += num_neigh_atoms_lst[jj];
    int* nloop_idxs = (int*)malloc(sizeof(int) * ntypes);
    
    int tot_num_neigh_atoms = 0;
    for (int ii=0; ii<ntypes; ii++)
        tot_num_neigh_atoms += num_neigh_atoms_lst[ii];

    // Step 2. Populate `tilde_s/x/y/z`
    for (int ii=0; ii<inum; ii++) {
        center_atom_idx = ilist[ii];
        center_cart_coords = x[center_atom_idx];
        for (int ii=0; ii<ntypes; ii++)
            nloop_idxs[ii] = 0;
            
        for (int jj=0; jj<numneigh[ii]; jj++) {
            neigh_atom_idx = firstneigh[ii][jj];
            
            // Step 3.1. 计算计算 1/r (`distance_ji_recip`), s(r_ji) (`tilde_s_value`)
            neigh_cart_coords = x[neigh_atom_idx];
            for (int kk=0; kk<3; kk++)
                diff_cart_coords[kk] = neigh_cart_coords[kk] - center_cart_coords[kk];
            distance_ji = vec3Operation::norm(diff_cart_coords);
            distance_ji_recip = recip<CoordType>(distance_ji);
            tilde_s_value = smooth_func(distance_ji, rcut, rcut_smooth);

            // Step 3.2. 计算 `x_ji_s`, `y_ji_s`, `z_ji_s` 
            tilde_x_value = tilde_s_value * distance_ji_recip * diff_cart_coords[0];
            tilde_y_value = tilde_s_value * distance_ji_recip * diff_cart_coords[1];
            tilde_z_value = tilde_s_value * distance_ji_recip * diff_cart_coords[2];
            
            // Step 3.3. Assignment
            int kk = 0;
            for (kk=0; kk<ntypes; kk++) 
                if (types[neigh_atom_idx] == kk)
                    break;
            printf("[%d: %d, %d]\n", kk, nstart_idxs[kk], nloop_idxs[kk]);
            tilde_r[0 + (nstart_idxs[kk]+nloop_idxs[kk])*4 + ii*tot_num_neigh_atoms*4] = tilde_s_value;
            tilde_r[1 + (nstart_idxs[kk]+nloop_idxs[kk])*4 + ii*tot_num_neigh_atoms*4] = tilde_x_value;
            tilde_r[2 + (nstart_idxs[kk]+nloop_idxs[kk])*4 + ii*tot_num_neigh_atoms*4] = tilde_y_value;
            tilde_r[3 + (nstart_idxs[kk]+nloop_idxs[kk])*4 + ii*tot_num_neigh_atoms*4] = tilde_z_value;
            nloop_idxs[kk]++;
        }
    }


    // Step. Free Memory
    free(center_cart_coords);
}

};
};


#endif