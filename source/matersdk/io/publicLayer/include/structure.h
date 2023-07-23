#ifndef MATERSDK_STRUCTURE_H
#define MATERSDK_STRUCTURE_H


namespace matersdk {

template <typename CoordType>
class Structure {
public:
    Structure();

    Structure(int num_atoms);
    
    Structure(int num_atoms, CoordType **basis_vectors, int *atomic_numbers, CoordType **coords, bool is_cart_coords=true);
    
    Structure(int num_atoms, CoordType basis_vectors[3][3], int atomic_number[], CoordType coords[][3], bool is_cart_coords=true);

    Structure(const Structure &rhs);

    Structure& operator=(const Structure &rhs);
    
    ~Structure();

    void calc_cart_coords(CoordType **frac_coords);

    void calc_cart_coords(CoordType frac_coords[][3]);

    // Note: `0~this->num_atoms` are owned atoms; others are ghost atoms.
    void make_supercell(const int *scaling_matix);   // Note: You can use `int[3]` to init it.

    // void make_supercell(const int scaling_matrix[3]);

    void show() const;


private:
    int num_atoms;
    CoordType **basis_vectors;
    int *atomic_numbers;
    CoordType **cart_coords;
}; // class: Structure



} // namespace: matersdk








// Definition of Structure member function
namespace matersdk {


template <typename CoordType>
Structure<CoordType>::Structure() {
    this->num_atoms = 0;
}


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
    // Step 1. Free memory
    // Step 1.1. `this->basis_vectors`
    for (int ii=0; ii<3; ii++) {
        free(this->basis_vectors[ii]);
    }
    free(this->basis_vectors);

    // Step 1.2. `this->atomic_numbers`
    free(this->atomic_numbers);

    // Step 1.3. `this->cart_coords`
    for (int ii=0; ii<this->num_atoms; ii++) {
        free(this->cart_coords[ii]);
    }
    free(this->cart_coords);
    // Step 1.4. `this->num_atoms = 0`
    this->num_atoms = 0;


    // Step 2. Reallocate and Reasign
    this->num_atoms = rhs.num_atoms;
    // Step 2.1. Allocate memory for `this->basis_vectors` and assign
    this->basis_vectors = (CoordType**)malloc(sizeof(CoordType*) * 3);
    for (int ii=0; ii<3; ii++) {
        this->basis_vectors[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    for (int ii=0; ii<3; ii++) {
        for (int jj=0; jj<3; jj++) {
            this->basis_vectors[ii][jj] = rhs.basis_vectors[ii][jj];
        }
    }

    // Step 2.2. Allocate memory for `this->atomic_numbers` and assign
    this->atomic_numbers = (int*)malloc(sizeof(int) * this->num_atoms);
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->atomic_numbers[ii] = rhs.atomic_numbers[ii];
    }

    // Step 2.3. Allocate memory for `this->cart_coords` and assign
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
 * @brief Copy assignment operator
 * 
 * @tparam CoordType 
 * @param rhs 
 * @return Structure<CoordType>& 
 */
template <typename CoordType>
Structure<CoordType>& Structure<CoordType>::operator=(const Structure &rhs) {
    // Step 1. Free memory 
    // Step 1.1. `this->basis_vectors`
    for (int ii=0; ii<3; ii++) {
        free(this->basis_vectors[ii]);
    }
    free(this->basis_vectors);

    // Step 1.2. `this->atomic_numbers`
    free(this->atomic_numbers);

    // Step 1.3. `this->cart_coords`
    for (int ii=0; ii<this->num_atoms; ii++) {
        free(this->cart_coords[ii]);
    }
    free(this->cart_coords);

    // Step 1.4. `this->num_atoms = 0`
    this->num_atoms = 0;


    // Step 2. Reallocate and reassign
    this->num_atoms = rhs.num_atoms;
    // Step 2.1. Allocate memory for `this->basis_vectors` and assign it
    this->basis_vectors = (CoordType**)malloc(sizeof(CoordType*) * 3);
    for (int ii=0; ii<3; ii++) {
        this->basis_vectors[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    for (int ii=0; ii<3; ii++) {
        this->basis_vectors[ii][0] = rhs.basis_vectors[ii][0];
        this->basis_vectors[ii][1] = rhs.basis_vectors[ii][1];
        this->basis_vectors[ii][2] = rhs.basis_vectors[ii][2];
    }

    // Step 2.2. Allocate memory for `this->atomic_numbers` and assign it
    this->atomic_numbers = (int*)malloc(sizeof(int) * this->num_atoms);
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->atomic_numbers[ii] = rhs.atomic_numbers[ii];
    }

    // Step 2.3. Allocate memory for `this->cart_coords` and assign it 
    this->cart_coords = (CoordType**)malloc(sizeof(CoordType*) * this->num_atoms);
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->cart_coords[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->cart_coords[ii][0] = rhs.cart_coords[ii][0];
        this->cart_coords[ii][1] = rhs.cart_coords[ii][1];
        this->cart_coords[ii][2] = rhs.cart_coords[ii][2];
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
 * @brief make supercell
 * 
 * @tparam CoordType 
 * @param scaling_matrix 
 */
template <typename CoordType>
void Structure<CoordType>::make_supercell(const int *scaling_matrix) {
    /*
        1. 奇数: 
            ( -\frac{num-1}{2}, \frac{num-1}{2})
        2. 偶数: 
            ( -(\frac{num}{2}+1), \frac{num}{2} )
    */
    int range[3][2];
    for (int ii=0; ii<3; ii++) {
        if (scaling_matrix[ii] % 2 == 0){   // 偶数
            range[ii][0] = -scaling_matrix[ii]/2 + 1;
            range[ii][1] = scaling_matrix[ii]/2;
        } else {    // 奇数
            range[ii][0] = -(scaling_matrix[ii]-1)/2;
            range[ii][1] = (scaling_matrix[ii]-1)/2;
        }
    }


    // Step 2. Allocate memory for `primitive_cell` and Reallocate memory for `supercell (this)`
    // Step 2.1. 利用 `num_atoms_prim`, `atomic_numbers_prim`, `cart_coords_prim` 存储原胞的信息
    int num_atoms_prim = this->num_atoms;
    // `atomic_numbers`
    int *atomic_numbers_prim = (int*)malloc(sizeof(int) * num_atoms_prim);
    for (int ii=0; ii<num_atoms_prim; ii++) {
        atomic_numbers_prim[ii] = this->atomic_numbers[ii];
    }
    // `cart_coords`
    CoordType **cart_coords_prim = (CoordType**)malloc(sizeof(CoordType*) * num_atoms_prim);
    for (int ii=0; ii<num_atoms_prim; ii++) {
        cart_coords_prim[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }
    for (int ii=0; ii<num_atoms_prim; ii++) {
        for (int jj=0; jj<3; jj++) {
            cart_coords_prim[ii][jj] = this->cart_coords[ii][jj];
        }
    }

    // Step 2.2. Reallocate memory for `this`
    this->num_atoms = num_atoms_prim * scaling_matrix[0] * scaling_matrix[1] * scaling_matrix[2];
    // `atomic_numbers`
    free(this->atomic_numbers);
    this->atomic_numbers = (int*)malloc(sizeof(int) * this->num_atoms);

    // `cart_coords`
    for (int ii=0; ii<num_atoms_prim; ii++) {
        free(this->cart_coords[ii]);
    }
    free(this->cart_coords);
    this->cart_coords = (CoordType**)malloc(sizeof(CoordType*) * this->num_atoms);
    for (int ii=0; ii<this->num_atoms; ii++) {
        this->cart_coords[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);
    }

    
    // Step 3. Reassign `basis_vectors`, `atomic_numbers`, `cart_coords`
    // Step 3.1. Calculate `atomic_numbers` and assign it to `this->atomic_numbers`
    for (int num_copies=0; num_copies<scaling_matrix[0]*scaling_matrix[1]*scaling_matrix[2]; num_copies++) {
        for (int atom_idx=0; atom_idx<num_atoms_prim; atom_idx++) {
            this->atomic_numbers[num_copies*num_atoms_prim + atom_idx] = atomic_numbers_prim[atom_idx];
        }
    }

    // Step 3.2. Calculate `cart_coords` and assign it to `this->cart_coords`
    int atom_idx = 0;
    for (int ii=range[0][0]; ii<range[0][1] + 1; ii++) {
        for (int jj=range[1][0]; jj<range[1][1] + 1; jj++) {
            for (int kk=range[2][0]; kk<range[2][1] + 1; kk++) {
                for (int prim_atom_idx=0; prim_atom_idx<num_atoms_prim; prim_atom_idx++) {
                    //std::cout << ii << ", " << jj << ", " << kk << std::endl;
                    this->cart_coords[atom_idx][0] = (
                        cart_coords_prim[prim_atom_idx][0] + 
                        this->basis_vectors[0][0] * ii + 
                        this->basis_vectors[1][0] * jj + 
                        this->basis_vectors[2][0] * kk
                    );
                    this->cart_coords[atom_idx][1] = (
                        cart_coords_prim[prim_atom_idx][1] + 
                        this->basis_vectors[0][1] * ii + 
                        this->basis_vectors[1][1] * jj +
                        this->basis_vectors[2][1] * kk
                    );
                    this->cart_coords[atom_idx][2] = (
                        cart_coords_prim[prim_atom_idx][2] + 
                        this->basis_vectors[0][2] * ii + 
                        this->basis_vectors[1][2] * jj + 
                        this->basis_vectors[2][2] * kk
                    );

                    atom_idx++;
                }
            }
        }
    }

    // Step 3.3. Calculate `basis_vectors` and assign it to `this->basis_vectors`
    for (int ii=0; ii<3; ii++) { // 三个基矢方向
        this->basis_vectors[ii][0] *= scaling_matrix[ii];
        this->basis_vectors[ii][1] *= scaling_matrix[ii];
        this->basis_vectors[ii][2] *= scaling_matrix[ii];
    }

}


/**
 * @brief Output the information of `Sturcture`
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void Structure<CoordType>::show() const {
    printf("Lattice\n");
    printf("------------------------------------------------\n");
    printf(" %-15.6f %-15.6f %-15.6f\n", this->basis_vectors[0][0], this->basis_vectors[0][1], this->basis_vectors[0][2]);
    printf(" %-15.6f %-15.6f %-15.6f\n", this->basis_vectors[1][0], this->basis_vectors[1][1], this->basis_vectors[1][2]);
    printf(" %-15.6f %-15.6f %-15.6f\n", this->basis_vectors[2][0], this->basis_vectors[2][1], this->basis_vectors[2][2]);
    printf("\nSite (Cartesian Coordinate)\n");
    printf("------------------------------------------------\n");
    for (int ii=0; ii<this->num_atoms; ii++)
        printf(" %-4d %-4d  %-15.6f %-15.6f %-15.6f\n", ii, this->atomic_numbers[ii], this->cart_coords[ii][0], this->cart_coords[ii][1], this->cart_coords[ii][2]);
}



}   // namespace: matersdk


#endif