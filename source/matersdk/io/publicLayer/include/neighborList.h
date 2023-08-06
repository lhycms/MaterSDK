#ifndef MATERSDK_NEIGHBOR_LIST_H
#define MATERSDK_NEIGHBOR_LIST_H

#include <cassert>
#include "./binLinkedList.h"




namespace matersdk {


/**
 * @brief A `neighbor list` class
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
class NeighborList {
public:
    friend class BinLinkedList<CoordType>;
    
    NeighborList(Structure<CoordType> structure, CoordType rcut, CoordType* bin_size_xyz, bool* pbc_xyz);

    ~NeighborList();

    void _build();

    void show() const;

    const BinLinkedList<CoordType>& get_binLinkedList() const;

private:
    BinLinkedList<CoordType> bin_linked_list;
    int num_atoms;                              // The number of atoms in primitive cell.
    std::vector<int>* neighbor_lists;           // `num_atoms` 个 `neighbor_lists` 构成完整的 `neighbor_list`
};  // class: NeighborList



/**
 * @brief Construct a new Neighbor List< Coord Type>:: Neighbor List object
 * 
 * @tparam CoordType 
 * @param structure 
 * @param rcut 
 * @param bin_size_xyz 
 * @param pbc_xyz 
 */
template <typename CoordType>
NeighborList<CoordType>::NeighborList(Structure<CoordType> structure, CoordType rcut, CoordType* bin_size_xyz, bool* pbc_xyz) {
    assert(structure.get_num_atoms() > 0);

    this->bin_linked_list = BinLinkedList<CoordType>(structure, rcut, bin_size_xyz, pbc_xyz);
    this->num_atoms = this->bin_linked_list.get_supercell().get_prim_num_atoms();
    this->neighbor_lists = new std::vector<int>[this->num_atoms];
    this->_build();     // Populate `this->neighbor_lists`
}


/**
 * @brief Destroy the Neighbor List< Coord Type>:: Neighbor List object
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
NeighborList<CoordType>::~NeighborList() {
    if (this->num_atoms != 0) {
        delete [] this->neighbor_lists;
        this->num_atoms = 0;
    }
}


/**
 * @brief Populate `this->neighbor_lists`
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void NeighborList<CoordType>::_build() {
    // Step 1. 得到 primitive cell 的相关信息
    // Step 1.1. primitive cell 中的原子个数、primitive cell 的 `prim_cell_idx`
    int prim_cell_idx = this->bin_linked_list.get_supercell().get_prim_cell_idx();
    int prim_num_atoms = this->bin_linked_list.get_supercell().get_prim_num_atoms();
    assert(prim_num_atoms == this->num_atoms);

    // Step 2. Populate `this->neighbor_lists`
    int center_atom_idx;    // `center_atom_index`: atom 在 supercell 中对应的 index
    int center_atomic_number;
    for (int ii=0; ii<this->num_atoms; ii++) {  // `ii`: atom 在 primitive cell 中的 index
        this->neighbor_lists[ii] = this->bin_linked_list.get_neigh_atoms(ii);
    }
}


template <typename CoordType>
void NeighborList<CoordType>::show() const {
    const int* supercell_atomic_numbers = this->bin_linked_list.get_supercell().get_structure().get_atomic_numbers();

    for (int ii=0; ii<this->num_atoms; ii++) {
        printf("%d:\t", ii);
        for (int jj: this->neighbor_lists[ii]) {
            printf("%d, ", jj);
        }
        printf("\n");
    }

    printf("\n\n");
    for (int ii=0; ii<this->num_atoms; ii++) {
        printf("%d:\t", supercell_atomic_numbers[ii]);
        for (int jj: this->neighbor_lists[ii]) {
            printf("%d, ", supercell_atomic_numbers[jj]);
        }
        printf("\n");
    }
}


template <typename CoordType>
const BinLinkedList<CoordType>& NeighborList<CoordType>::get_binLinkedList() const {
    return this->bin_linked_list;
}


}; // namespace: matersdk


#endif