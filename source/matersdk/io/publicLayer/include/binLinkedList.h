#ifndef MATERSDK_BIN_LINKED_LIST_H
#define MATERSDK_BIN_LINKED_LIST_H


#include "./structure.h"


namespace matersdk {


template <typename CoordType>
class Supercell {
public:
    Supercell();

    Supercell(Structure<CoordType>& structure, int *scaling_matrix);

    //Supercell(Structure<CoordType> Structure, int scaling_matrix[3]);
    
    //Supercell(const Supercell &rhs);

    //Supercell& operator=(const Supercell &rhs);

    ~Supercell();

    void calc_prim_cell_idx_xyz();

    void calc_prim_cell_idx();      // Call this function after `this->calc_prim_cell_idx_xyz()`

    void calc_owned_atom_idxs();    // Call this function after `this->`

    void show() const;              // Print out information

    const int* get_prim_cell_idx_xyz() const;

    const int get_prim_cell_idx() const;

    const int get_prim_num_atoms() const;

    const int get_num_atoms() const;

    const int* get_owned_atom_idxs() const;


private:
    Structure<CoordType> structure;
    int scaling_matrix[3];      // 扩包倍数；x, y, z 方向上的 primitive_cell 个数
    int num_atoms = 0;
    int prim_num_atoms;         // primitive cell 的元素数目
    int prim_cell_idx;          // primitive cell 对应的 cell index
    int prim_cell_idx_xyz[3];   // 
    int *owned_atom_idxs;       // 

}; // class: Supercell




template <typename CoordType>
matersdk::Supercell<CoordType>::Supercell() {
    this->structure = Structure<CoordType>();

    this->scaling_matrix[0] = 1;
    this->scaling_matrix[1] = 1;
    this->scaling_matrix[2] = 1;

    this->prim_num_atoms = 0;
    this->prim_cell_idx = 0;

    this->prim_cell_idx_xyz[0] = 0;
    this->prim_cell_idx_xyz[1] = 0;
    this->prim_cell_idx_xyz[2] = 0;

    // owned_atom_idxs
}


/**
 * @brief Construct a new matersdk::Supercell<Coord Type>::Supercell object
 * 
 * @tparam CoordType 
 * @param structure 
 * @param scaling_matrix 
 */
template <typename CoordType>
matersdk::Supercell<CoordType>::Supercell(Structure<CoordType>& structure, int *scaling_matrix)
{
    this->structure = structure;
    for (int ii=0; ii<3; ii++) {
        this->scaling_matrix[ii] = scaling_matrix[ii];
    }

    this->num_atoms = this->structure.get_num_atoms() * this->scaling_matrix[0] * this->scaling_matrix[1] * this->scaling_matrix[2];
    this->prim_num_atoms = this->structure.get_num_atoms();
    this->calc_prim_cell_idx_xyz();  // Assign `this->prim_cell_idx_xyz`
    this->calc_prim_cell_idx();      // Assign `this->prim_cell_idx`
    this->owned_atom_idxs = (int*)malloc(sizeof(int) * this->prim_num_atoms);
    this->calc_owned_atom_idxs();

    // Step 3. make_supercell
    this->structure.make_supercell(this->scaling_matrix);
}


/**
 * @brief Destroy the Supercell< Coord Type>:: Supercell object
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
Supercell<CoordType>::~Supercell() {
    if (this->num_atoms != 0)
        free(this->owned_atom_idxs);
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
 * @brief Calculate the `this->prim_cell_idx` and assign it.
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


/**
 * @brief Calculate the `this->owned_atom_idxs` and assign it.
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void Supercell<CoordType>::calc_owned_atom_idxs() {
    for (int ii=0; ii<this->prim_num_atoms; ii++) {
        this->owned_atom_idxs[ii] = this->prim_cell_idx * this->prim_num_atoms + ii;
    }
}



/**
 * @brief Print out the information of `this` (class = Supercell)
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void Supercell<CoordType>::show() const {
    this->structure.show();
    printf("\n");
    printf("scaling_matrix = [%4d, %4d, %4d]\n", this->scaling_matrix[0], this->scaling_matrix[1], this->scaling_matrix[2]);
    printf("num_atoms = %15d\n", this->num_atoms);
    printf("prim_num_atoms = %15d\n", this->prim_num_atoms);
    printf("prim_cell_idx = %15d\n", this->prim_cell_idx);
    printf("prim_cell_idx_xyz = [%4d, %4d, %4d]\n", this->prim_cell_idx_xyz[0], this->prim_cell_idx_xyz[1], this->prim_cell_idx_xyz[2]);
    if (this->num_atoms != 0)
        printf("owned_atom_idxs range = %15d ~ %15d\n", this->owned_atom_idxs[0], this->owned_atom_idxs[this->prim_num_atoms-1]);
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


template <typename CoordType>
const int* Supercell<CoordType>::get_owned_atom_idxs() const {
    return (const int*)this->owned_atom_idxs;
}


} // namespace: matersdk

#endif