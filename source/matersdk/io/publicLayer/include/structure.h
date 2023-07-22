#ifndef MATERSDK_STRUCTURE_H
#define MATERSDK_STRUCTURE_H


namespace matersdk {

template <typename CoordType>
class Structure {
public:
    Structure(int num_atoms);
    
    Structure(int num_atoms, CoordType **basis_vectors, int *atomic_numbers, CoordType **coords, bool is_cart_coords=true);
    
    Structure(int num_atoms, CoordType basis_vectors[3][3], int atomic_number[], CoordType coords[][3], bool is_cart_coords=true);

    Structure(const Structure &rhs);
    
    ~Structure();

    void calc_cart_coords(CoordType **frac_coords);

    void calc_cart_coords(CoordType frac_coords[][3]);
    
    //void make_supercell(int *scaling_matix);
    //

    void show();


private:
    int num_atoms;
    CoordType **basis_vectors;
    int *atomic_numbers;
    CoordType **cart_coords;
}; // class: Structure



} // namespace: matersdk








// Definition of Structure member function
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
        CoordType **basis_vectors, int *atomic_numbers, CoordType **coords,
        bool is_cart_coords)
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
    if (is_cart_coords) {
        for (int ii=0; ii<num_atoms; ii++) {
            for (int jj=0; jj<3; jj++) {
                this->cart_coords[ii][jj] = coords[ii][jj];
            }
        }
    } else {
        this->calc_cart_coords(coords);
    }
}



/**
 * @brief Construct a new Structure< Coord Type>:: Structure object
 * 
 * @tparam CoordType 
 * @param num_atoms 
 * @param basis_vectors 
 * @param atomic_numbers 
 * @param coords 
 * @param is_cart_coords 
 */
template <typename CoordType>
Structure<CoordType>::Structure(int num_atoms,
        CoordType basis_vectors[3][3], int atomic_numbers[], CoordType coords[][3],
        bool is_cart_coords)
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
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->atomic_numbers[ii] = atomic_numbers[ii];
    }

    // Step 3. Allocate memory for `this->cart_coords`
    this->cart_coords = (CoordType**)malloc(sizeof(CoordType*) * this->num_atoms);
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->cart_coords[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    if (is_cart_coords) {   // 如果 `coords` 是笛卡尔坐标
        for (int ii=0; ii<this->num_atoms; ii++) {
            for (int jj=0; jj<3; jj++) {
                this->cart_coords[ii][jj] = coords[ii][jj];
            }
        }
    } else {    // 若如果不是笛卡尔坐标
        this->calc_cart_coords(coords);
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


/**
 * @brief Convert the `fractional coordinates` to `cartesian coordinates`
 * 
 * @tparam CoordType 
 * @param frac_coords 
 */
template <typename CoordType>
void Structure<CoordType>::calc_cart_coords(CoordType **frac_coords) {
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->cart_coords[ii][0] = (
            frac_coords[ii][0] * this->basis_vectors[0][0] + 
            frac_coords[ii][1] * this->basis_vectors[1][0] + 
            frac_coords[ii][2] * this->basis_vectors[2][0]
        );
        this->cart_coords[ii][1] = (
            frac_coords[ii][0] * this->basis_vectors[0][1] + 
            frac_coords[ii][1] * this->basis_vectors[1][1] + 
            frac_coords[ii][2] * this->basis_vectors[2][1]
        );
        this->cart_coords[ii][2] = (
            frac_coords[ii][0] * this->basis_vectors[0][2] + 
            frac_coords[ii][1] * this->basis_vectors[1][2] +
            frac_coords[ii][2] * this->basis_vectors[2][2]
        );
    }
}


/**
 * @brief Convert the `fractional coordinates` to `cartesian coordinates`
 * 
 * @tparam CoordType 
 * @param frac_coords 
 */
template <typename CoordType>
void Structure<CoordType>::calc_cart_coords(CoordType frac_coords[][3]) {
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->cart_coords[ii][0] = (
            frac_coords[ii][0] * this->basis_vectors[0][0] + 
            frac_coords[ii][1] * this->basis_vectors[1][0] + 
            frac_coords[ii][2] * this->basis_vectors[2][0]
        );
        this->cart_coords[ii][1] = (
            frac_coords[ii][0] * this->basis_vectors[0][1] +
            frac_coords[ii][1] * this->basis_vectors[1][1] + 
            frac_coords[ii][2] * this->basis_vectors[2][1]
        );
        this->cart_coords[ii][2] = (
            frac_coords[ii][0] * this->basis_vectors[0][2] + 
            frac_coords[ii][1] * this->basis_vectors[1][2] +
            frac_coords[ii][2] * this->basis_vectors[2][2]
        );
    }
}


/**
 * @brief Output the information of `Sturcture`
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void Structure<CoordType>::show() {
    printf("Lattice\n");
    printf("------------------------------------------------\n");
    printf(" %-15.6f %-15.6f %-15.6f\n", this->basis_vectors[0][0], this->basis_vectors[0][1], this->basis_vectors[0][2]);
    printf(" %-15.6f %-15.6f %-15.6f\n", this->basis_vectors[1][0], this->basis_vectors[1][1], this->basis_vectors[1][2]);
    printf(" %-15.6f %-15.6f %-15.6f\n", this->basis_vectors[2][0], this->basis_vectors[2][1], this->basis_vectors[2][2]);
    printf("\nSite (Cartesian Coordinate)\n");
    printf("------------------------------------------------\n");
    for (int ii=0; ii<this->num_atoms; ii++)
        printf(" %-4d  %-15.6f %-15.6f %-15.6f\n", this->atomic_numbers[ii], this->cart_coords[ii][0], this->cart_coords[ii][1], this->cart_coords[ii][2]);
}



}   // namespace: matersdk


#endif