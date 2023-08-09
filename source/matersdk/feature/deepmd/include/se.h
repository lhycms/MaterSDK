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

    PairTildeR(NeighborList<CoordType>& neighbor_list, int center_atomic_number, int neigh_atomic_number);

    void show() const;

private:
    NeighborList<CoordType> neighbor_list;
    int center_atomic_number = 0;
    int neigh_atomic_number = 0;
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
}


template <typename CoordType>
void PairTildeR<CoordType>::show() const {
    if (this->center_atomic_number == 0)
        printf("This is a NULL PairTildeR.\n");
    else {
        int max_num_neigh_atoms = this->neighbor_list.get_max_num_neigh_atoms_ss(this->neigh_atomic_number);
        printf("center_atomic_number = %d\n", this->center_atomic_number);
        printf("neigh_atomic_number = %d\n", this->neigh_atomic_number);
        printf("max_num_neigh_atoms_ss = %d\n", max_num_neigh_atoms);
    }
}


}   // namespace : deepPotSE
}   // namespace : matersdk


#endif