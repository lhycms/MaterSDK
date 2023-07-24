#ifndef MATERSDK_BIN_LINKED_LIST_H
#define MATERSDK_BIN_LINKED_LIST_H


#include "./structure.h"


namespace matersdk {


template <typename CoordType>
class Supercell {
public:
    Supercell(Structure<CoordType>& structure, int *scaling_matrix);

    //Supercell(Structure<CoordType> Structure, int scaling_matrix[3]);
    
    //Supercell(const Structure &rhs);

    //~Supercell();
    
    //void show() const;

    void calc_prim_cell_idx_xyz();

    void calc_prim_cell_idx();   // Call this function after `this->calc_prim_cell_idx_xyz()`

    const int* get_prim_cell_idx_xyz() const;

    const int get_prim_cell_idx() const;

    const int get_prim_num_atoms() const;

    const int get_num_atoms() const;

private:
    Structure<CoordType> structure;
    int scaling_matrix[3];      // 扩包倍数；x, y, z 方向上的 primitive_cell 个数
    int prim_num_atoms;         // primitive cell 的元素数目
    int prim_cell_idx;          // primitive cell 对应的 cell index
    int prim_cell_idx_xyz[3];   // 
    int *owned_atom_idxs;       // 

}; // class: Supercell




template <typename CoordType>
matersdk::Supercell<CoordType>::Supercell(Structure<CoordType>& structure, int *scaling_matrix)
{
    this->structure = structure;
    for (int ii=0; ii<3; ii++) {
        this->scaling_matrix[ii] = scaling_matrix[ii];
    }

    this->prim_num_atoms = this->structure.get_num_atoms();
    this->calc_prim_cell_idx_xyz();  // Assign `this->prim_cell_idx_xyz`
    this->calc_prim_cell_idx();      // Assign `this->prim_cell_idx`
    
    // Step 3. make_supercell
    this->structure.make_supercell(this->scaling_matrix);
}



/**
 * @brief Calculate the `this->prim_cell_idx_xyz` and assign it.
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void Supercell<CoordType>::calc_prim_cell_idx_xyz() {
    /*
        1. scaling_factor = 奇数
            - central_idx = \frac{scaling_factor - 1}{2}
        2. scaling_factor = 偶数
            - central_idx = \frac{scaling_factor}{2} - 1
    */
    for (int ii=0; ii<3; ii++) {
        if (this->scaling_matrix[ii] % 2 != 0) {    // 奇数
            this->prim_cell_idx_xyz[ii] = (this->scaling_matrix[ii] - 1) / 2;
        } else {    // 偶数
            this->prim_cell_idx_xyz[ii] = this->scaling_matrix[ii] / 2 - 1;
        }
    }
}


/**
 * @brief Calculate the `this->prim_cell_idx`
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void Supercell<CoordType>::calc_prim_cell_idx() {
    this->prim_cell_idx = (
                this->prim_cell_idx_xyz[0] + 
                this->prim_cell_idx_xyz[1] * this->scaling_matrix[0] + 
                this->prim_cell_idx_xyz[2] * this->scaling_matrix[0] * this->scaling_matrix[1]
    );
}


template <typename CoordType>
const int* Supercell<CoordType>::get_prim_cell_idx_xyz() const {
    return (const int*)this->prim_cell_idx_xyz;   // Type conversion : `int[]` -> `int*`
}


template <typename CoordType>
const int Supercell<CoordType>::get_prim_cell_idx() const {
    return (const int)this->prim_cell_idx;
}


template <typename CoordType>
const int Supercell<CoordType>::get_prim_num_atoms() const {
    return (const int)this->prim_num_atoms;
}


template <typename CoordType>
const int Supercell<CoordType>::get_num_atoms() const {
    return (const int)(this->prim_num_atoms * this->scaling_matrix[0] * this->scaling_matrix[1] * this->scaling_matrix[2]);
}



} // namespace: matersdk

#endif