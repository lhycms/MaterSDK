#ifndef MATERSDK_SE_H
#define MATERSDK_SE_H


#include <stdlib.h>
#include "../../../io/publicLayer/include/neighborList.h"


namespace matersdk {
namespace deepPotSE{


template <typename CoordType>
class PairTildeR {
public:
    PairTildeR();

    PairTildeR(NeighborList<CoordType>& neighbor_list, int center_atomic_number, int neigh_atomic_number, int num_neigh_atoms);

    PairTildeR(NeighborList<CoordType>& neighbor_list, int center_atomic_number, int neigh_atomic_number);

    void show() const;

    const int get_max_num_neigh_atoms() const;  // 指定 `center_atomic_number`, `neigh_atomic_number`, 计算最大近邻原子数

    const int get_num_center_atoms() const;

    const int get_num_neigh_atoms() const;      // 相当于指定了 zero-padding 的尺寸

    void generate(int num_neigh_atoms);         // 产生 $\tilde{R^i}$ 特征

private:
    NeighborList<CoordType> neighbor_list;
    int center_atomic_number = 0;
    int neigh_atomic_number = 0;
    int num_neigh_atoms = 0;                    // 相当于指定了 zero-padding 的尺寸
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
    this->num_neigh_atoms = 0;
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
PairTildeR<CoordType>::PairTildeR(NeighborList<CoordType>& neighbor_list, int center_atomic_number, int neigh_atomic_number, int num_neigh_atoms) {
    this->neighbor_list = neighbor_list;
    this->center_atomic_number = center_atomic_number;
    this->neigh_atomic_number = neigh_atomic_number;
    this->num_neigh_atoms = num_neigh_atoms;
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
PairTildeR<CoordType>::PairTildeR(NeighborList<CoordType>& neighbor_list, int center_atomic_number, int neigh_atomic_number) {
    this->neighbor_list = neighbor_list;
    this->center_atomic_number = center_atomic_number;
    this->neigh_atomic_number = neigh_atomic_number;
    this->num_neigh_atoms = this->get_max_num_neigh_atoms();
}


template <typename CoordType>
void PairTildeR<CoordType>::show() const {
    if (this->center_atomic_number == 0)
        printf("This is a NULL PairTildeR.\n");
    else {
        int max_num_neigh_atoms = this->neighbor_list.get_max_num_neigh_atoms_ssss(this->center_atomic_number, this->neigh_atomic_number);
        printf("center_atomic_number = %d\n", this->center_atomic_number);
        printf("neigh_atomic_number = %d\n", this->neigh_atomic_number);
        printf("num_center_atoms = %d\n", this->get_num_center_atoms());
        printf("num_neigh_atoms = %d\n", this->get_num_neigh_atoms());
        printf("max_num_neigh_atoms_ss = %d\n", this->get_max_num_neigh_atoms());
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
 * @brief 计算 `this->num_center_atoms`
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
const int PairTildeR<CoordType>::get_num_center_atoms() const {
    int num_center_atoms = 0;
    int prim_num_atoms = this->neighbor_list.get_binLinkedList().get_supercell().get_prim_num_atoms();
    const int* supercell_atomic_numbers = this->neighbor_list.get_binLinkedList().get_supercell().get_structure().get_atomic_numbers();

    for (int ii=0; ii<prim_num_atoms; ii++) {
        if (supercell_atomic_numbers[ii] == this->center_atomic_number)
            num_center_atoms++;
    }

    return num_center_atoms;
}


/**
 * @brief `this->num_neigh_atoms` 指定了 zero-padding 的尺寸
 * 
 * @tparam CoordType 
 * @return const int 
 */
template <typename CoordType>
const int PairTildeR<CoordType>::get_num_neigh_atoms() const {
    return this->num_neigh_atoms;
}


}   // namespace : deepPotSE
}   // namespace : matersdk


#endif