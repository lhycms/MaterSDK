#ifndef MATERSDK_SE_H
#define MATERSDK_SE_H


#include <stdlib.h>

#include "../../../io/publicLayer/include/structure.h"
#include "../../../io/publicLayer/include/neighborList.h"
#include "../../../../core/include/vec3Operation.h"


namespace matersdk {
namespace deepPotSE{


/**
 * @brief Smooth function is DeepPot-SE
 * 
 * @tparam CoordType 
 * @param distance_ji 
 * @param rcut 
 * @param rcut_smooth 
 * @return CoordType 
 */
template <typename CoordType>
CoordType smooth_func(const CoordType& distance_ji, const CoordType& rcut, const CoordType& rcut_smooth) {
    CoordType smooth_value;     // return value
    // Step 1. calculate `r_recip`
    CoordType r_recip = 0;
    if (distance_ji == 0)
        r_recip = 0;
    else
        r_recip = 1 / distance_ji;
    
    // Step 2. uu = (r - r_s) / (r_c - r_s)
    CoordType uu = (distance_ji - rcut_smooth) / (rcut - rcut_smooth);
    
    // Step 3. Calculate `smooth_value`
    if (distance_ji < rcut_smooth)
        smooth_value = r_recip;
    else if ( (distance_ji >= rcut_smooth) && (distance_ji < rcut) )
        smooth_value = r_recip * ( std::pow(uu, 3) * (-6*std::pow(uu, 2) + 15*uu - 10) + 1);
    else
        smooth_value = 0;
    
    return smooth_value;
}


/**
 * @brief Calculate the reciprocal value.
 * 
 * @tparam CoordType 
 * @param value 
 * @return CoordType 
 */
template <typename CoordType>
CoordType recip(const CoordType& value) {
    CoordType value_recip;
    if (value == 0) {
        value_recip = 0;
    } else {
        value_recip = 1.0 / value;
    }
    return value_recip;
}


/**
 * @brief 指定 `center_atomic_number` 和 `neigh_atomic_number`，计算 DeepPot-SE 的 feature
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
class PairTildeR {
public:
    PairTildeR();

    PairTildeR(
                NeighborList<CoordType>& neighbor_list, 
                int center_atomic_number, 
                int neigh_atomic_number, 
                int num_neigh_atoms,
                CoordType rcut_smmoth);

    PairTildeR(
                NeighborList<CoordType>& neighbor_list, 
                int center_atomic_number, 
                int neigh_atomic_number,
                CoordType rcut_smooth);

    PairTildeR(
                Structure<CoordType>& structure,
                CoordType rcut,
                bool* pbc_xyz,
                bool sort,
                int center_atomic_number,
                int neigh_atomic_number,
                int num_neigh_atoms,
                CoordType rcut_smooth);
    
    PairTildeR(
                Structure<CoordType>& structure,
                CoordType rcut,
                bool* pbc_xyz,
                bool sort,
                int center_atomic_number,
                int neigh_atomic_number,
                CoordType rcut_smooth);

    void calc_num_center_atoms();

    const int get_num_center_atoms() const;     // 得到中心原子的数目

    const int get_num_neigh_atoms() const;      // 得到指定的近邻原子数目（`this->num_neigh_atoms`相当于指定了 zero-padding 的尺寸）

    void show() const;

    void show_in_value() const;

    void show_in_deriv() const;

    const int get_max_num_neigh_atoms() const;  // 得到最大近邻原子数目 （指定 `center_atomic_number`, `neigh_atomic_number`, 计算最大近邻原子数）

    CoordType*** generate() const;         // 计算 $\tilde{R^i}$ 特征 .shape = (num_center_atoms, num_neigh_atoms, 4)

    CoordType**** deriv() const;           // 计算 $\tilde{R^i}$ 特征的导数 .shape = (num_center_atoms, num_neigh_atoms, 4, 3)


private:
    NeighborList<CoordType> neighbor_list;
    int center_atomic_number = 0;
    int neigh_atomic_number = 0;
    int num_center_atoms = 0;
    int num_neigh_atoms = 0;                    // 相当于指定了 zero-padding 的尺寸
    CoordType rcut = 0;
    CoordType rcut_smooth = 0;
};  // class PairTildeR




/**
 * @brief Construct a new Pair Tilde R< Coord Type>:: Pair Tilde R object
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
PairTildeR<CoordType>::PairTildeR() {
    this->neighbor_list = NeighborList<CoordType>();
    this->center_atomic_number = 0;
    this->neigh_atomic_number = 0;
    this->num_center_atoms = 0;
    this->num_neigh_atoms = 0;
    this->rcut = 0;
    this->rcut_smooth = 0;
}


/**
 * @brief Construct a new Pair Tilde R< Coord Type>:: Pair Tilde R object
 * 
 * @tparam CoordType 
 * @param neighbor_list 
 * @param center_atomic_number 
 * @param neigh_atomic_number 
 * @param num_neigh_atoms 
 */
template <typename CoordType>
PairTildeR<CoordType>::PairTildeR(
                            NeighborList<CoordType>& neighbor_list, 
                            int center_atomic_number, 
                            int neigh_atomic_number, 
                            int num_neigh_atoms,
                            CoordType rcut_smooth) {
    this->neighbor_list = neighbor_list;
    this->center_atomic_number = center_atomic_number;
    this->neigh_atomic_number = neigh_atomic_number;
    this->calc_num_center_atoms();  // 计算 `this->num_center_atoms`
    this->num_neigh_atoms = num_neigh_atoms;
    this->rcut = this->neighbor_list.get_rcut();
    this->rcut_smooth = rcut_smooth;
}


/**
 * @brief Construct a new Pair Tilde R< Coord Type>:: Pair Tilde R object
 * 
 * @tparam CoordType 
 * @param neighbor_list 
 * @param center_atomic_number 
 * @param neigh_atomic_number 
 */
template <typename CoordType>
PairTildeR<CoordType>::PairTildeR(NeighborList<CoordType>& neighbor_list, int center_atomic_number, int neigh_atomic_number, CoordType rcut_smooth) {
    this->neighbor_list = neighbor_list;
    this->center_atomic_number = center_atomic_number;
    this->neigh_atomic_number = neigh_atomic_number;
    this->calc_num_center_atoms();  // 计算 `this->num_center_atoms`
    this->num_neigh_atoms = this->get_max_num_neigh_atoms();
    this->rcut = this->neighbor_list.get_rcut();
    this->rcut_smooth = rcut_smooth;
}


/**
 * @brief Construct a new Pair Tilde R< Coord Type>:: Pair Tilde R object
 * 
 * @tparam CoordType 
 * @param structure 
 * @param rcut 
 * @param pbc_xyz 
 * @param sort 
 * @param center_atomic_number 
 * @param neigh_atomic_number 
 * @param num_neigh_atoms 
 * @param rcut_smooth 
 */
template <typename CoordType>
PairTildeR<CoordType>::PairTildeR(
                        Structure<CoordType>& structure,
                        CoordType rcut,
                        bool* pbc_xyz,
                        bool sort,
                        int center_atomic_number,
                        int neigh_atomic_number,
                        int num_neigh_atoms,
                        CoordType rcut_smooth)
{
    this->neighbor_list = NeighborList<CoordType>(structure, rcut, pbc_xyz, sort);
    this->center_atomic_number = center_atomic_number;
    this->neigh_atomic_number = neigh_atomic_number;
    this->calc_num_center_atoms();
    this->num_neigh_atoms = num_neigh_atoms;
    this->rcut = this->neighbor_list.get_rcut();
    this->rcut_smooth = rcut_smooth;
}


/**
 * @brief Construct a new Pair Tilde R< Coord Type>:: Pair Tilde R object
 * 
 * @tparam CoordType 
 * @param structure 
 * @param rcut 
 * @param pbc_xyz 
 * @param sort 
 * @param center_atomic_number 
 * @param neigh_atomic_number 
 * @param rcut_smooth 
 */
template <typename CoordType>
PairTildeR<CoordType>::PairTildeR(
                        Structure<CoordType>& structure,
                        CoordType rcut,
                        bool* pbc_xyz,
                        bool sort,
                        int center_atomic_number,
                        int neigh_atomic_number,
                        CoordType rcut_smooth)
{
    this->neighbor_list = NeighborList<CoordType>(structure, rcut, pbc_xyz, sort);
    this->center_atomic_number = center_atomic_number;
    this->neigh_atomic_number = neigh_atomic_number;
    this->calc_num_center_atoms();
    this->num_neigh_atoms = this->get_max_num_neigh_atoms();
    this->rcut = this->neighbor_list.get_rcut();
    this->rcut_smooth = rcut_smooth;
}


/**
 * @brief 计算 `this->num_center_atoms`
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void PairTildeR<CoordType>::calc_num_center_atoms() {
    int num_center_atoms = 0;
    int prim_num_atoms = this->neighbor_list.get_binLinkedList().get_supercell().get_prim_num_atoms();
    const int* supercell_atomic_numbers = this->neighbor_list.get_binLinkedList().get_supercell().get_structure().get_atomic_numbers();

    for (int ii=0; ii<prim_num_atoms; ii++) {
        if (supercell_atomic_numbers[ii] == this->center_atomic_number)
            num_center_atoms++;
    }

    this->num_center_atoms = num_center_atoms;
}


template <typename CoordType>
const int PairTildeR<CoordType>::get_num_center_atoms() const {
    return this->num_center_atoms;
}


template <typename CoordType>
const int PairTildeR<CoordType>::get_num_neigh_atoms() const {
    return this->num_neigh_atoms;
}


template <typename CoordType>
void PairTildeR<CoordType>::show() const {
    if (this->center_atomic_number == 0)
        printf("This is a NULL PairTildeR.\n");
    else {
        int max_num_neigh_atoms = this->neighbor_list.get_max_num_neigh_atoms_ssss(this->center_atomic_number, this->neigh_atomic_number);
        printf("center_atomic_number = %d\n", this->center_atomic_number);
        printf("neigh_atomic_number = %d\n", this->neigh_atomic_number);
        printf("num_center_atoms = %d\n", this->num_center_atoms);
        printf("num_neigh_atoms = %d\n", this->num_neigh_atoms);
        printf("max_num_neigh_atoms_ss = %d\n", this->get_max_num_neigh_atoms());
        printf("rcut = %f\n", this->rcut);
        printf("rcut_smooth = %f\n", this->rcut_smooth);
    }
}


template <typename CoordType>
void PairTildeR<CoordType>::show_in_value() const {
    if (this->center_atomic_number == 0)
        printf("This is a NULL PairTildeR.\n");
    else {
        CoordType*** pair_tilde_r = this->generate();

        for (int ii=0; ii<this->num_center_atoms; ii++) {
            for (int jj=0; jj<this->num_neigh_atoms; jj++) {
                printf("[%d, %d] -- [%f, %f, %f, %f]\n", 
                        ii, jj,
                        pair_tilde_r[ii][jj][0],
                        pair_tilde_r[ii][jj][1],
                        pair_tilde_r[ii][jj][2],
                        pair_tilde_r[ii][jj][3]
                );
            }
        }
    }
}


/**
 * @brief 指定 `center_atomic_number`, `neigh_atomic_number`, 计算最大近邻原子数
 * 
 * @tparam CoordType 
 * @return const int 
 */
template <typename CoordType>
const int PairTildeR<CoordType>::get_max_num_neigh_atoms() const {
    return this->neighbor_list.get_max_num_neigh_atoms_ssss(this->center_atomic_number, this->neigh_atomic_number);
}



/**
 * @brief 
 * 
 * @tparam CoordType 
 * @return CoordType***     : shape = (num_center_atoms, num_neigh_atoms, 4)
 */
template <typename CoordType>
CoordType*** PairTildeR<CoordType>::generate() const {
    // Step 1. 
    // Step 1.1. $\tilde{R}$ = (s(r_{ji}), x_{ji}, y_{ji}, z_{ji})
    //  = (tilde_s_value, tilde_x_value, tilde_y_value, tilde_z_value)
    CoordType tilde_s_value;
    CoordType tilde_x_value;
    CoordType tilde_y_value;
    CoordType tilde_z_value;

    int prim_num_atoms = this->neighbor_list.get_binLinkedList().get_supercell().get_prim_num_atoms();
    int prim_cell_idx = this->neighbor_list.get_binLinkedList().get_supercell().get_prim_cell_idx();

    int center_atom_idx;            // 中心原子在 supercell 中的索引
    int neigh_atom_idx;             // 近邻原子在 supercell 中的索引
    const CoordType* center_cart_coord;   // 中心原子的坐标
    const CoordType* neigh_cart_coord;    // 近邻原子的坐标
    CoordType* diff_cart_coord = (CoordType*)malloc(sizeof(CoordType) * 3);     // 近邻原子的坐标 - 中心原子的坐标
    CoordType distance_ji;          // 两原子间的距离
    CoordType distance_ji_recip;

    // Step 1.2. Allocate memory for $\tilde{R}$
    CoordType*** pair_tilde_r = (CoordType***)malloc(sizeof(CoordType**) * this->num_center_atoms);
    for (int ii=0; ii<this->num_center_atoms; ii++) {
        pair_tilde_r[ii] = (CoordType**)malloc(sizeof(CoordType*) * this->neigh_atomic_number);
        for (int jj=0; jj<this->num_neigh_atoms; jj++) {
            pair_tilde_r[ii][jj] = (CoordType*)malloc(sizeof(CoordType) * 4);
        }
    }

    // Step 2. ·
    const CoordType** supercell_cart_coords = this->neighbor_list.get_binLinkedList().get_supercell().get_structure().get_cart_coords();
    const int* supercell_atomic_numbers = this->neighbor_list.get_binLinkedList().get_supercell().get_structure().get_atomic_numbers();

    // Step 3. Populate `tilde_s/x/y/z_value`
    int tmp_cidx;   // 中心原子 for loop
    int tmp_nidx;   // 近邻原子 for loop

    tmp_cidx = 0;
    for (int ii=0; ii<this->neighbor_list.get_num_center_atoms(); ii++) {
        center_atom_idx = ii + prim_cell_idx * prim_num_atoms;
        center_cart_coord = supercell_cart_coords[center_atom_idx];
        // 若中心原子 不等于 `this->center_atomic_number`，直接跳过
        if (supercell_atomic_numbers[center_atom_idx] != this->center_atomic_number)
            continue;
        
        tmp_nidx = 0;
        for (int jj=0; jj<this->neighbor_list.get_neighbor_lists()[ii].size(); jj++) {
            neigh_atom_idx = this->neighbor_list.get_neighbor_lists()[ii][jj];
            
            // 若近邻原子 不等于 `this->neigh_atomic_number`，直接跳过
            if (supercell_atomic_numbers[neigh_atom_idx] != this->neigh_atomic_number)
                continue;

            // Step 3.1. 计算 1/r (`distance_ji_recip`), s(r_ji) (`tilde_s_value`)
            neigh_cart_coord = supercell_cart_coords[neigh_atom_idx];
            for (int kk=0; kk<3; kk++)
                diff_cart_coord[kk] = neigh_cart_coord[kk] - center_cart_coord[kk];
            distance_ji = vec3Operation::norm(diff_cart_coord);
            distance_ji_recip = recip(distance_ji);
            tilde_s_value = smooth_func(distance_ji, this->rcut, this->rcut_smooth);

            // Step 3.2. 计算 `x_ji_s`, `y_ji_s`, `z_ji_s` 
            tilde_x_value = tilde_s_value * distance_ji_recip * diff_cart_coord[0];
            tilde_y_value = tilde_s_value * distance_ji_recip * diff_cart_coord[1];
            tilde_z_value = tilde_s_value * distance_ji_recip * diff_cart_coord[2];

            // Step 3.3. Assignment
            pair_tilde_r[tmp_cidx][tmp_nidx][0] = tilde_s_value;
            pair_tilde_r[tmp_cidx][tmp_nidx][1] = tilde_x_value;
            pair_tilde_r[tmp_cidx][tmp_nidx][2] = tilde_y_value;
            pair_tilde_r[tmp_cidx][tmp_nidx][3] = tilde_z_value;
            //printf("[%d, %d] -- [%f, %f, %f, %f]\n", tmp_cidx, tmp_nidx, pair_tilde_r[tmp_cidx][tmp_nidx][0], pair_tilde_r[tmp_cidx][tmp_nidx][1], pair_tilde_r[tmp_cidx][tmp_nidx][2], pair_tilde_r[tmp_cidx][tmp_nidx][3]);

            tmp_nidx++;
        }

        tmp_cidx++;
    }

    assert(tmp_cidx == this->num_center_atoms);
    assert(tmp_nidx == this->num_neigh_atoms);

    // Step . Free memory
    free(diff_cart_coord);

    return pair_tilde_r;
}


}   // namespace : deepPotSE
}   // namespace : matersdk


#endif