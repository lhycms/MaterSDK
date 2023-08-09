#ifndef MATERSDK_SE_H
#define MATERSDK_SE_H


#include <stdlib.h>
#include "../../../io/publicLayer/include/structure.h"
#include "../../../io/publicLayer/include/neighborList.h"


namespace matersdk {
namespace deepPotSE{


template <typename CoordType>
CoordType smooth_func(CoordType r_recip) {

}


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

    // 一般不会使用，速度比较慢
    PairTildeR(
                Structure<CoordType>& structure,
                CoordType rcut,
                bool* pbc_xyz,
                bool sort,
                int center_atomic_number,
                int neigh_atomic_number,
                int num_neigh_atoms,
                CoordType rcut_smooth);
    
    // 一般不会使用，速度比较慢
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

    const int get_max_num_neigh_atoms() const;  // 得到最大近邻原子数目 （指定 `center_atomic_number`, `neigh_atomic_number`, 计算最大近邻原子数）

    CoordType** generate() const;         // 产生 $\tilde{R^i}$ 特征

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
 * @brief Generate $\tilde{R}$, $\tilde{R}$ .shape = (this->num_center_atoms, this->num_neigh_atoms)
 * 
 * @tparam CoordType 
 * @return CoordType** 
 */
template <typename CoordType>
CoordType** PairTildeR<CoordType>::generate() const {
    // Step 1. Allocate memory for $\tilde{R}$
    CoordType** tilde_r = (CoordType**)malloc(sizeof(CoordType*) * this->num_center_atoms);
    for (int ii=0; ii<this->num_center_atoms; ii++) 
        tilde_r[ii] = (CoordType*)malloc(sizeof(CoordType) * this->neigh_atomic_number);

    // Step 2. 
    for (int ii=0; ii<this->num_center_atoms; ii++) {
        for (int jj=0; jj<this->num_neigh_atoms; jj++) {

        }
    }

}


}   // namespace : deepPotSE
}   // namespace : matersdk


#endif