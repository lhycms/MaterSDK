#include <stdlib.h>

#include "../include/structure.h"


namespace matersdk {


/**
 * @brief Construct a new Structure< Coord Type>:: Structure object
 * 
 * @tparam CoordType 
 * @param num_atoms 
 */
template <typename CoordType>
Structure<CoordType>::Structure(int num_atoms) {
    // Step 1. Allocate memory for `this->basis_vectors`
    this->basis_vectors = (CoordType**)malloc(sizeof(CoordType*) * 3);
    for (int ii=0; ii<3; ii++) {
        this->basis_vectors[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }

    // Step 2. Allocate memory for `this->atomic_numbers`
    this->atomic_numbers = (int*)malloc(sizeof(int) * num_atoms);

    // Step 3. Allocate memory for `this->cart_coords`
    this->cart_coords = (CoordType**)malloc(sizeof(CoordType*) * num_atoms);
    for (int ii=0; ii<num_atoms; ii++) {
        this->cart_coords[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
}


/**
 * @brief Construct a new Structure< Coord Type>:: Structure object
 * 
 * @tparam CoordType 
 * @param num_atoms 
 * @param basis_vectors 
 * @param atomic_numbers 
 * @param cart_coords 
 */
template <typename CoordType>
Structure<CoordType>::Structure(
        int num_atoms,
        CoordType **basis_vectors, int *atomic_numbers, CoordType **cart_coords)
{
    // Step 1. Allocate memory for `this->basis_vectors` and assign
    this->basis_vectors = (CoordType**)malloc(sizeof(CoordType*) * 3);
    for (int ii=0; ii<3; ii++) {
        this->basis_vectors[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    // Step 1.1. Assign `this->basis_vectors`
    for (int ii=0; ii<3; ii++) {
        for (int jj=0; jj<3; jj++) {
            this->basis_vectors[ii][jj] = basis_vectors[ii][jj];
        }
    }

    // Step 2. Allocate memory for `this->atomic_numbers` and assign
    this->atomic_numbers = (int*)malloc(sizeof(int) * num_atoms);
    // Step 2.1. Assign `this->atomic_numbers`
    for (int ii=0; ii<num_atoms; ii++) {
        this->atomic_numbers[ii] = atomic_numbers[ii];
    }

    // Step 3. Allocate memory for `this->cart_coords` and assign
    this->cart_coords = (CoordType**)malloc(sizeof(CoordType*) * num_atoms);
    for (int ii=0; ii<num_atoms ; ii++) {
        this->cart_coords[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    for (int ii=0; ii<num_atoms; ii++) {
        for (int jj=0; jj<3; jj++) {
            this->cart_coords[ii][jj] = cart_coords[ii][jj];
        }
    }

}



}   // namespace: matersdk