#ifndef MATERSDK_BIN_LINKED_LIST_H
#define MATERSDK_BIN_LINKED_LIST_H


#include "./structure.h"
#include "../../../../core/include/vec3Operation.h"


namespace matersdk {


template <typename CoordType>
class Supercell {
public:
    Supercell();

    Supercell(Structure<CoordType>& structure, int *scaling_matrix); // Note `Structure<CoordType> &structure` is a reference.

    //Supercell(Structure<CoordType> Structure, int scaling_matrix[3]);
    
    Supercell(const Supercell &rhs);

    Supercell& operator=(const Supercell &rhs);

    ~Supercell();

    void calc_prim_cell_idx_xyz();

    void calc_prim_cell_idx();      // Call this function after `this->calc_prim_cell_idx_xyz()`

    void calc_owned_atom_idxs();    // Call this function after `this->`

    void calc_prim_projected_lengths();

    void calc_prim_inter_planar_distances();

    void show() const;              // Print out information

    const int* get_prim_cell_idx_xyz() const;

    const int get_prim_cell_idx() const;

    const int get_prim_num_atoms() const;

    const int get_num_atoms() const;

    const int* get_owned_atom_idxs() const;


private:
    Structure<CoordType> structure;
    int scaling_matrix[3] = {1, 1, 1};      // 扩包倍数；x, y, z 方向上的 primitive_cell 个数
    int num_atoms = 0;
    CoordType **prim_basis_vectors;
    int prim_num_atoms = 0;         // primitive cell 的元素数目
    int prim_cell_idx = 0;          // primitive cell 对应的 cell index
    int prim_cell_idx_xyz[3] = {0, 0, 0};   // 
    CoordType prim_projected_lengths[3] = {0, 0, 0};
    CoordType prim_inter_planar_distances[3] = {0, 0, 0};
    int *owned_atom_idxs;           // 

}; // class: Supercell




/**
 * @brief Construct a new matersdk::Supercell<Coord Type>::Supercell object
 * 
 * @tparam CoordType 
 */
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

    // prim_basis_vectors
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

    this->prim_basis_vectors = (CoordType**)malloc(sizeof(CoordType*) * 3);
    for (int ii=0; ii<3; ii++) {
        this->prim_basis_vectors[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    for (int ii=0; ii<3; ii++) {
        this->prim_basis_vectors[ii][0] = this->structure.basis_vectors[ii][0];
        this->prim_basis_vectors[ii][1] = this->structure.basis_vectors[ii][1];
        this->prim_basis_vectors[ii][2] = this->structure.basis_vectors[ii][2];
    }

    this->prim_num_atoms = this->structure.get_num_atoms();
    this->calc_prim_cell_idx_xyz();             // Assign `this->prim_cell_idx_xyz`
    this->calc_prim_cell_idx();                 // Assign `this->prim_cell_idx`
    this->calc_prim_projected_lengths();        // Assign `this->prim_projected_lengths`
    this->calc_prim_inter_planar_distances();   // Assign `this->prim_inter_planar_distances`
    this->owned_atom_idxs = (int*)malloc(sizeof(int) * this->prim_num_atoms);
    this->calc_owned_atom_idxs();               // Assign `this->prim_owned_atom_idxs`

    // Step 3. make_supercell
    this->structure.make_supercell(this->scaling_matrix);
}


/**
 * @brief Construct a new Supercell< Coord Type>:: Supercell object
 * 
 * @tparam CoordType 
 * @param rhs 
 */
template <typename CoordType>
Supercell<CoordType>::Supercell(const Supercell &rhs) {
    this->structure = rhs.structure;
    for (int ii=0; ii<3; ii++) {
        this->scaling_matrix[ii] = rhs.scaling_matrix[ii];
    }
    this->num_atoms = rhs.num_atoms;

    this->prim_num_atoms = rhs.prim_num_atoms;
    this->prim_cell_idx = rhs.prim_cell_idx;
    this->prim_cell_idx_xyz[0] = rhs.prim_cell_idx_xyz[0];
    this->prim_cell_idx_xyz[1] = rhs.prim_cell_idx_xyz[1];
    this->prim_cell_idx_xyz[2] = rhs.prim_cell_idx_xyz[2];

    this->prim_projected_lengths[0] = rhs.prim_projected_lengths[0];
    this->prim_projected_lengths[1] = rhs.prim_projected_lengths[1];
    this->prim_projected_lengths[2] = rhs.prim_projected_lengths[2];

    this->prim_inter_planar_distances[0] = rhs.prim_inter_planar_distances[0];
    this->prim_inter_planar_distances[1] = rhs.prim_inter_planar_distances[1];
    this->prim_inter_planar_distances[2] = rhs.prim_inter_planar_distances[2];

    this->owned_atom_idxs = (int*)malloc(sizeof(int) * rhs.prim_num_atoms);

    if (this->num_atoms != 0) {
        this->prim_basis_vectors = (CoordType**)malloc(sizeof(CoordType*) * 3);
        for (int ii=0; ii<3; ii++) {
            this->prim_basis_vectors[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
        }
        for (int ii=0; ii<3; ii++) {
            this->prim_basis_vectors[ii][0] = rhs.prim_basis_vectors[ii][0];
            this->prim_basis_vectors[ii][1] = rhs.prim_basis_vectors[ii][1];
            this->prim_basis_vectors[ii][2] = rhs.prim_basis_vectors[ii][2];
        }

        for (int ii=0; ii<rhs.prim_num_atoms; ii++) {
            this->owned_atom_idxs[ii] = rhs.owned_atom_idxs[ii];
        }
    }

}


template <typename CoordType>
Supercell<CoordType>& Supercell<CoordType>::operator=(const Supercell<CoordType> &rhs) {
    if (this->num_atoms != 0) {
        for (int ii=0; ii<3; ii++) {
            free(this->prim_basis_vectors[ii]);
        }
        free(this->prim_basis_vectors);

        free(this->owned_atom_idxs);
    }
    
    this->structure = rhs.structure;
    for (int ii=0; ii<3; ii++) {
        this->scaling_matrix[ii] = rhs.scaling_matrix[ii];
    }
    this->num_atoms = rhs.num_atoms;

    this->prim_num_atoms = rhs.prim_num_atoms;
    this->prim_cell_idx = rhs.prim_cell_idx;
    this->prim_cell_idx_xyz[0] = rhs.prim_cell_idx_xyz[0];
    this->prim_cell_idx_xyz[1] = rhs.prim_cell_idx_xyz[1];
    this->prim_cell_idx_xyz[2] = rhs.prim_cell_idx_xyz[2];

    this->prim_projected_lengths[0] = rhs.prim_projected_lengths[0];
    this->prim_projected_lengths[1] = rhs.prim_projected_lengths[1];
    this->prim_projected_lengths[2] = rhs.prim_projected_lengths[2];

    this->prim_inter_planar_distances[0] = rhs.prim_inter_planar_distances[0];
    this->prim_inter_planar_distances[1] = rhs.prim_inter_planar_distances[1];
    this->prim_inter_planar_distances[2] = rhs.prim_inter_planar_distances[2];

    if (rhs.num_atoms != 0) {
        this->prim_basis_vectors = (CoordType**)malloc(sizeof(CoordType*) * 3);
        for (int ii=0; ii<3; ii++) {
            this->prim_basis_vectors[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
        }
        for (int ii=0; ii<3; ii++) {
            this->prim_basis_vectors[ii][0] = rhs.prim_basis_vectors[ii][0];
            this->prim_basis_vectors[ii][1] = rhs.prim_basis_vectors[ii][1];
            this->prim_basis_vectors[ii][2] = rhs.prim_basis_vectors[ii][2];
        }

        this->owned_atom_idxs = (int*)malloc(sizeof(int) * rhs.prim_num_atoms);
        for (int ii=0; ii<rhs.prim_num_atoms; ii++) {
            this->owned_atom_idxs[ii] = rhs.owned_atom_idxs[ii];
        }
    }
    
    return *this;
}

/**
 * @brief Destroy the Supercell< Coord Type>:: Supercell object
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
Supercell<CoordType>::~Supercell() {
    if (this->num_atoms != 0) {
        for (int ii=0; ii<3; ii++) {
            free(this->prim_basis_vectors[ii]);
        }
        free(this->prim_basis_vectors);

        free(this->owned_atom_idxs);
    }
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
 * @brief Calculate the `this->prim_projected_lengths` and assign it.
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void Supercell<CoordType>::calc_prim_projected_lengths() {
    CoordType* unit_vector_x = (CoordType*)malloc(sizeof(CoordType) * 3);
    unit_vector_x[0] = 1;
    unit_vector_x[1] = 0;
    unit_vector_x[2] = 0;
    CoordType* unit_vector_y = (CoordType*)malloc(sizeof(CoordType) * 3);
    unit_vector_y[0] = 0;
    unit_vector_y[1] = 1;
    unit_vector_y[2] = 0;
    CoordType* unit_vector_z = (CoordType*)malloc(sizeof(CoordType) * 3);
    unit_vector_z[0] = 0;
    unit_vector_z[1] = 0;
    unit_vector_z[2] = 1;

    this->prim_projected_lengths[0] = (
        std::abs( vec3Operation::dot(this->prim_basis_vectors[0], unit_vector_x) ) +
        std::abs( vec3Operation::dot(this->prim_basis_vectors[1], unit_vector_x) ) +
        std::abs( vec3Operation::dot(this->prim_basis_vectors[2], unit_vector_x) )
    );
    this->prim_projected_lengths[1] = (
        std::abs( vec3Operation::dot(this->prim_basis_vectors[0], unit_vector_y) ) + 
        std::abs( vec3Operation::dot(this->prim_basis_vectors[1], unit_vector_y) ) + 
        std::abs( vec3Operation::dot(this->prim_basis_vectors[2], unit_vector_y) )
    );
    this->prim_projected_lengths[2] = (
        std::abs( vec3Operation::dot(this->prim_basis_vectors[0], unit_vector_z) ) + 
        std::abs( vec3Operation::dot(this->prim_basis_vectors[1], unit_vector_z) ) +
        std::abs( vec3Operation::dot(this->prim_basis_vectors[2], unit_vector_z) )
    );
}


template <typename CoordType>
void Supercell<CoordType>::calc_prim_inter_planar_distances() {
    CoordType *vec_vertical_yz = vec3Operation::normalize(
                                    vec3Operation::cross(
                                        this->prim_basis_vectors[1],
                                        this->prim_basis_vectors[2])
                                );
    CoordType *vec_vertical_xz = vec3Operation::normalize( 
                                    vec3Operation::cross(
                                        this->prim_basis_vectors[0],
                                        this->prim_basis_vectors[2])
                                );
    CoordType *vec_vertical_xy = vec3Operation::normalize(
                                    vec3Operation::cross(
                                        this->prim_basis_vectors[0],
                                        this->prim_basis_vectors[1])
                                );

    
    this->prim_inter_planar_distances[0] = 
        std::abs( vec3Operation::dot(this->prim_basis_vectors[0], vec_vertical_yz) );
    this->prim_inter_planar_distances[1] = 
        std::abs( vec3Operation::dot(this->prim_basis_vectors[1], vec_vertical_xz) );
    this->prim_inter_planar_distances[2] = 
        std::abs( vec3Operation::dot(this->prim_basis_vectors[2], vec_vertical_xy) );

}


/**
 * @brief Print out the information of `this` (class = Supercell)
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void Supercell<CoordType>::show() const {
    this->structure.show();
    if (this->num_atoms != 0) {
        printf("prim_basis_vectors:\n");
        printf("[%15f, %15f, %15f]\n", this->prim_basis_vectors[0][0], this->prim_basis_vectors[0][1], this->prim_basis_vectors[0][2]);
        printf("[%15f, %15f, %15f]\n", this->prim_basis_vectors[1][0], this->prim_basis_vectors[1][1], this->prim_basis_vectors[1][2]);
        printf("[%15f, %15f, %15f]\n", this->prim_basis_vectors[2][0], this->prim_basis_vectors[2][1], this->prim_basis_vectors[2][2]);
    }
    printf("prim_num_atoms = %15d\n", this->prim_num_atoms);
    printf("prim_cell_idx = %15d\n", this->prim_cell_idx);
    printf("prim_cell_idx_xyz = [%4d, %4d, %4d]\n", this->prim_cell_idx_xyz[0], this->prim_cell_idx_xyz[1], this->prim_cell_idx_xyz[2]);
    if (this->num_atoms != 0)
        printf("owned_atom_idxs range = %15d ~ %15d\n", this->owned_atom_idxs[0], this->owned_atom_idxs[this->prim_num_atoms-1]);
    printf("prim_projected_lengths = [%15f, %15f, %15f]\n", this->prim_projected_lengths[0], this->prim_projected_lengths[1], this->prim_projected_lengths[2]);
    printf("prim_inter_planar_distances = [%15f, %15f, %15f]\n", this->prim_inter_planar_distances[0], this->prim_inter_planar_distances[1], this->prim_inter_planar_distances[2]);
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