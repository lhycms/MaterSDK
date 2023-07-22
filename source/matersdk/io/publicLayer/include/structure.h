#ifndef MATERSDK_STRUCTURE_H
#define MATERSDK_STRUCTURE_H


namespace matersdk {

template <typename CoordType>
class Structure {
public:
    Structure(int num_atoms);
    
    Structure(int num_atoms, CoordType **basis_vectors, int *atomic_numbers, CoordType **cart_coords);
    
    Structure(int num_atoms, CoordType basis_vectors[3][3], int atomic_number[], CoordType cart_coords[][3]);

    Structure(const Structure &rhs);
    
    ~Structure();
    
    //void make_supercell(int *scaling_matix);
    //

private:
    int num_atoms;
    CoordType **basis_vectors;
    int *atomic_numbers;
    CoordType **cart_coords;
}; // class: Structure



} // namespace: matersdk








// Definition of Structure function
namespace matersdk {

/**
 * @brief Construct a new Structure< Coord Type>:: Structure object
 * 
 * @tparam CoordType 
 * @param num_atoms 
 */
template <typename CoordType>
Structure<CoordType>::Structure(int num_atoms) {
    this->num_atoms = num_atoms;
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
    this->num_atoms = num_atoms;
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
    // Step 3.1. Assign
    for (int ii=0; ii<num_atoms; ii++) {
        for (int jj=0; jj<3; jj++) {
            this->cart_coords[ii][jj] = cart_coords[ii][jj];
        }
    }
}


template <typename CoordType>
Structure<CoordType>::Structure(int num_atoms,
        CoordType basis_vectors[3][3], int atomic_numbers[], CoordType cart_coords[][3])
{
    this->num_atoms = num_atoms;
    // Step 1. Allocate memory for `this->basis_vectors` and assign
    this->basis_vectors = (CoordType**)malloc(sizeof(CoordType*) * 3);
    for (int ii=0; ii<3; ii++) {
        this->basis_vectors[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    for (int ii=0; ii<3; ii++) {
        for (int jj=0; jj<3; jj++) {
            this->basis_vectors[ii][jj] = basis_vectors[ii][jj];
        }
    }

    // Step 2. Allocate memory for `this->atomic_numbers`
    this->atomic_numbers = (int*)malloc(sizeof(int) * this->num_atoms);

    // Step 3. Allocate memory for `this->cart_coords`
    this->cart_coords = (CoordType**)malloc(sizeof(CoordType*) * this->num_atoms);
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->cart_coords[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    for (int ii=0; ii<this->num_atoms; ii++) {
        for (int jj=0; jj<3; jj++) {
            this->cart_coords[ii][jj] = cart_coords[ii][jj];
        }
    }
}


/**
 * @brief Construct a new Structure< Coord Type>:: Structure object
 * 
 * @tparam CoordType 
 * @param rhs 
 */
template <typename CoordType>
Structure<CoordType>::Structure(const Structure &rhs)
{  
    this->num_atoms = rhs.num_atoms;
    // Step 1. Allocate memory for `this->basis_vectors` and assign
    this->basis_vectors = (CoordType**)malloc(sizeof(CoordType*) * 3);
    for (int ii=0; ii<3; ii++) {
        this->basis_vectors[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    for (int ii=0; ii<3; ii++) {
        for (int jj=0; jj<3; jj++) {
            this->basis_vectors[ii][jj] = rhs.basis_vectors[ii][jj];
        }
    }

    // Step 2. Allocate memory for `this->atomic_numbers` and assign
    this->atomic_numbers = (int*)malloc(sizeof(int) * this->num_atoms);
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->atomic_numbers[ii] = rhs.atomic_numbers[ii];
    }

    // Step 3. Allocate memory for `this->cart_coords` and assign
    this->cart_coords = (CoordType**)malloc(sizeof(CoordType*) * this->num_atoms);
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->cart_coords[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    for (int ii=0; ii<this->num_atoms; ii++) {
        for (int jj=0; jj<3; jj++) {
            this->cart_coords[ii][jj] = rhs.cart_coords[ii][jj];
        }
    }

}


/**
 * @brief Destroy the Structure< Coord Type>:: Structure object
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
Structure<CoordType>::~Structure() {
    // Step 1. Deallocate `this->basis_vectors`
    for (int ii=0; ii<3; ii++) {
        free(this->basis_vectors[ii]);
    }
    free(this->basis_vectors);

    // Step 2. Deallocate `this->atomic_numbers`
    free(this->atomic_numbers);

    // Step 3. Deallocate `this->cart_coords`
    for (int ii=0; ii<this->num_atoms; ii++) {
        free(this->cart_coords[ii]);
    }
    free(this->cart_coords);

    // Step 4. `this->num_atoms = 0`
    this->num_atoms = 0;
}



}   // namespace: matersdk




#endif