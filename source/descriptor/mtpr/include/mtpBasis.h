#ifndef MATERSDK_MTPR_MTP_BASIS_H
#define MATERSDK_MTPR_MTP_BASIS_H

#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "./basis.h"
#include "./mtpParam.h"

namespace matersdk {
namespace mtpr {

template <typename CoordType>
class MtpBasis {
public:
    static void find_val_der(
        CoordType *mtp_val,
        CoordType (*mtp_der)[3],
        CoordType *mtp_der2coeffs,
        int chebyshev_size,
        CoordTyoe *coeffs,
        int alpha_moments_count,
        int alpha_index_basic_count,
        int (*alpha_index_basic)[4],
        int alpha_index_times_count,
        int (*alpha_index_times)[4],
        int alpha_scalar_moments,
        int *alpha_moment_mapping,
        int nmus,
        int inum,
        int *ilist,
        int *numneigh,
        int *firstneigh,
        CoordType (*relative_coords)[3],
        int *types,
        int ntypes,
        int umax_num_neigh_atoms,
        CoordType rmax,
        CoordType rmin);
};  // classs : MtpBasis


template <typename CoordType>
void MtpBasis<CoordType>::find_val_der(
    CoordType *mtp_basis_val,
    CoordType (*mtp_basis_der)[3],
    CoordType *mtp_basis_der2coeffs,
    int chebyshev_size,
    CoordTyoe *coeffs,
    int alpha_moments_count,
    int alpha_index_basic_count,
    int (*alpha_index_basic)[4],
    int alpha_index_times_count,
    int (*alpha_index_times)[4],
    int alpha_scalar_moments,
    int *alpha_moment_mapping,
    int nmus,
    int inum,
    int *ilist,
    int *numneigh,
    int *firstneigh,
    CoordType (*relative_coords)[3],
    int *types,
    int ntypes,
    int umax_num_neigh_atoms,
    CoordType rmax,
    CoordType rmin)
{
    // Step 1.
    memset(mtp_basis_val, 0, sizeof(CoordType) * alpha_scalar_moments);
    memset(mtp_basis_der, 0, sizeof(CoordType) * alpha_scalar_moments * 3);
    memset(mtp_basis_der2coeffs, 0, sizeof(CoordType) * alpha_scalar_moments * ntypes * ntypes * nmus * chebyshev_size);

    CoordType *mom_vals = (CoordType*)malloc(sizeof(CoordType) * alpha_moments_count);
    CoordType (*mom_ders)[3] = (CoordType*)malloc(sizeof(CoordType) * alpha_moments_count * 3);
    CoordType *mom_ders2coeffs = (CoordType*)malloc(sizeof(CoordType) * alpha_moments_count * ntypes * ntypes * nmus * chebyshev_size);
    memset(mom_vals, 0, sizeof(CoordType) * alpha_moments_count);
    memset(mom_ders, 0, sizeof(CoordType) * alpha_moments_count * 3);
    memset(mom_ders2coeffs, 0 , sizeof(CoordType) * alpha_moments_count * ntypes * ntypes * nmus * chebyshev_size);

    int max_alpha_index_basic = 0;
    for (int ii=0; ii<alpha_index_basic_count; ii++) {
        int now_alpha_index_basic = alpha_index_basic[ii][0] + alpha_index_basic[ii][1] + alpha_index_basic[ii][2];
        if (now_alpha_index_basic > max_alpha_index_basic) 
            max_alpha_index_basic += 1;
    }
    max_alpha_index_basic++;

    CoordType *auto_dist_powers_;
    CoordType (*auto_coords_powers_)[3];
    auto_dist_powers_ = (CoordType*)malloc(sizeof(CoordType) * max_alpha_index_basic);
    auto_coords_powers_ = (CoordType (*)[3])malloc(sizeof(CoordType) * max_alpha_index_basic * 3);
    memset(auto_dist_powers_, 0, sizeof(CoordType) * max_alpha_index_basic);
    memset(auto_coords_powers_, 0, sizeof(CoordType) * max_alpha_index_basic * 3);

    CoordType NeighbVect[3];
    CoordType distance_ij;
    int type_central;
    int type_outer;
    int num_coeffs = ntypes*ntypes*nmus*chebyshev_size;

    RQ_Chebyshev<CoordType>* p_RadialBasis = new RQ_Chebyshev<CoordType>(chebyshev_size, rmax, rmin);

    // Step 2.
    for (int ii=0; ii<inum; ii++) 
    {
        type_central = types[ilist[ii]];
        if (type_central >= ntypes)
            throw MtpException("Too few types in the MTP potential.");

        for (int jj=0; jj<numneigh[ii]; jj++) 
        {
            type_outer = types[firstneigh[ii*umax_num_neigh_atoms + jj]];
            if (type_outer >= ntypes)
                throw MtpException("Too few types in the MTP potential.");
            for (int a=0; a<3; a++)
                NeighbVect[a] = relative_coords[ii*umax_num_neigh_atoms*3 + jj*3 + a];
            distance_ij = std::sqrt( std::pow(NeighbVect[0], 2) + 
                                     std::pow(NeighbVect[1], 2) + 
                                     std::pow(NeighbVect[2], 2) );
            p_RadialBasis->build(distance_ij);
            
            auto_dist_powers_[0] = 1;
            for (int a=0; a<3; a++)
                auto_coords_powers_[0*3+a] = 1;
            for (int k=1; k<max_alpha_index_basic; k++) {
                auto_dist_powers_[k] = auto_dist_powers_[k-1] * distance_ij;
                for (int a=0; a<3; a++)
                    auto_coords_powers_[k*3 + a] = auto_coords_powers_[(k-1)*3 + a] * NeighbVect[a];
            }

            for (int i=0; i<alpha_index_basic_count; i++)
            {
                int mu = alpha_index_basic[i][0];
                int k = alpha_index_basic[i][1] + alpha_index_basic[i][2] + alpha_index_basic[i][3];
                CoordType powk = 1 / auto_dist_powers_[k];
                CoordType pow0 = auto_coords_powers_[alpha_index_basic[i][1]][0];
                CoordType pow1 = auto_coords_powers_[alpha_index_basic[i][2]][1];
                CoordType pow2 = auto_coords_powers_[alpha_index_basic[i][3]][2];
                CoordType mult0 = pow0*pow1*pow2;

                for (int xi=0; xi<chebyshev_size; xi++) {
                    int idx = (type_central*ntypes + type_outer)*nmus*chebyshev_size + mu*chebyshev_size + xi;
                    mom_vals[i] += coeffs[idx] * p_RadialBasis->vals()[xi] * powk * mult0;
                    mom_ders2coeffs[i*num_coeffs + idx] += p_RadialBasis->vals()[xi] * powk * mult0;
                    mom_ders[i][0] += NeighbVect[0] / distance_ij * 
                                    (coeffs[idx] * p_RadialBasis->ders2r()[xi] * powk * mult0 
                                    - coeffs[idx] * p_RadialBasis->vals()[xi] * k * powk / distance_ij * mult0);
                    mom_ders[i][1] += NeighbVect[1] / distance_ij * 
                                    (coeffs[idx] * p_RadialBasis->ders2r()[xi] * powk * mult0
                                    - coeffs[idx] * p_RadialBasis->vals()[xi] * k * powk / distance_ij * mult0);
                    mom_ders[i][2] += NeighbVect[2] / distance_ij *
                                    (coeffs[idx] * p_RadialBasis->ders2r()[xi] * powk * mult0
                                    - coeffs[idx] * p_RadialBasis->vals()[xi] * k * powk / distance_ij * mult0);
                    if (alpha_index_basic[i][1] != 0) {
                        mom_ders[i][0] += coeffs[idx] * p_RadialBasis->vals()[xi] * powk * alpha_index_basic[i][1]
                                        * auto_coords_powers_[alpha_index_basic[i][1] - 1][0]
                                        * pow1
                                        * pow2;
                    }
                    if (alpha_index_basic[i][2] != 0) {
                        mom_ders[i][1] += coeffs[idx] * p_RadialBasis->vals()[xi] * powk * alpha_index_basic[i][2]
                                        * pow0
                                        * auto_coords_powers_[alpha_index_basic[i][2] - 1][1]
                                        * pow2;
                    }
                    if (alpha_index_basic[i][3] != 0) {
                        mom_ders[i][2] += coeffs[idx] * p_RadialBasis->vals()[xi] * powk * alpha_index_basic[i][3]
                                        * pow0
                                        * pow1
                                        * auto_coords_powers_[alpha_index_basic[i][3] - 1][2];
                    }
                }
            }
        }

        for (int i=0; i<alpha_index_times_count; i++)
        {
            CoordType val0 = mom_vals[alpha_index_times[i][0]];
            CoordType val1 = mom_vals[alpha_index_times[i][1]];
            CoordType val2 = alpha_index_times[i][2];

            mom_vals[alpha_index_times[i][3]] += val2 * val0 * val1;
            for (int xi=0; xi<chebyshev_size; xi++) {
                int idx0 = (type_central*ntypes + type_outer)*nmus*chebyshev_size + alpha_index_times[i][0]*chebyshev_size + xi;
                int idx1 = (type_central*ntypes + type_outer)*nmus*chebyshev_size + alpha_index_times[i][1]*chebyshev_size + xi;
                mom_ders2coeffs[alpha_index_times[i][3]*num_coeffs + idx0] += val2 
                    * val1
                    * mom_ders2coeffs[alpha_index_times[i][0]*num_coeffs + idx0];
                mom_ders2coeffs[alpha_index_times[i][3]*num_coeffs + idx1] += val2 
                    * val0
                    * mom_ders2coeffs[alpha_index_times[i][0]*num_coeffs + idx1];
            }
            mom_ders[alpha_index_times[i][3]][0] += val2 * 
                    ( mom_ders[alpha_index_times[i][0]][0] * val1
                    + val0 * mom_ders[alpha_index_times[i][1]][0] );
            mom_ders[alpha_index_times[i][3]][1] += val2 * 
                    ( mom_ders[alpha_index_times[i][0]][1] * val1
                    + val0 * mom_ders[alpha_index_times[i][1]][1] );
            mom_ders[alpha_index_times[i][3]][2] += val2 * 
                    ( mom_ders[alpha_index_times[i][0]][2] * val1
                    + val0 * mom_ders[alpha_index_times[i][2]][2] );
        }

        for (int i=0; i<alpha_scalar_moments; i++) 
        {
            mtp_basis_val[i] = mom_vals[alpha_moment_mapping[i]];
            for (int a=0; a<3; a++)
                mtp_basis_der[i][a] = mom_ders[alpha_moment_mapping[i]][a];
            for (int idx=0; idx<num_coeffs; idx++)
                mtp_basis_der2coeffs[i*num_coeffs + idx] = mom_ders2coeffs[alpha_moment_mapping[i]*num_coeffs + idx];
        }
    }
    
    // Step . Free memory
    free(mom_vals);
    free(mom_ders);
    free(mom_ders2coeffs);
    free(auto_dist_powers_);
    free(auto_coords_powers_);
    delete p_RadialBasis;
}

};  // namespace : mtpr
};  // namespace : matersdk

#endif