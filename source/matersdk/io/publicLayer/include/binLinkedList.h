#ifndef MATERSDK_BIN_LINKED_LIST_H
#define MATERSDK_BIN_LINKED_LIST_H


#include "./structure.h"
#include "../../../../core/include/vec3Operation.h"


namespace matersdk {

template <typename CoordType>
class BinLinkedList;


template <typename CoordType>
class BasicStructureInfo {
public:
    BasicStructureInfo();

    BasicStructureInfo(Structure<CoordType> &structure);

    BasicStructureInfo& operator=(const BasicStructureInfo& rhs);

    BasicStructureInfo(const BasicStructureInfo& rhs);

    ~BasicStructureInfo();

    void show() const;

    friend class Supercell<CoordType>;

private:
    int num_atoms = 0;
    CoordType projected_lengths[3] = {0, 0, 0};
    CoordType interplanar_distances[3] = {0, 0, 0};
};



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

    void show() const;              // Print out information

    const int get_prim_num_atoms() const;

    const int* get_prim_cell_idx_xyz() const;

    const int get_prim_cell_idx() const;

    const int get_num_atoms() const;

    const int* get_owned_atom_idxs() const;

    friend class BinLinkedList<CoordType>;


private:
    Structure<CoordType> structure;
    BasicStructureInfo<CoordType> prim_structure_info;
    int scaling_matrix[3] = {1, 1, 1};      // 扩包倍数；x, y, z 方向上的 primitive_cell 个数
    int num_atoms = 0;
    int prim_cell_idx = 0;          // primitive cell 对应的 cell index
    int prim_cell_idx_xyz[3] = {0, 0, 0};   // 
    int *owned_atom_idxs;           // 
}; // class: Supercell



template <typename CoordType>
class BinLinkedList {
public:
    BinLinkedList();

    BinLinkedList(Structure<CoordType>& structure, CoordType rcut, CoordType* bin_size_xyz, bool* pbc_xyz);

    // BinLinkedList(const BinLinkedList& rhs);

    // BinLinkedList& operator=(const BinLinkedList& rhs);

    // ~BinLinkedList();

    int get_bin_idx(int prim_atom_idx);


private:
    Supercell<CoordType> supercell;
    CoordType rcut = 0;
    int bin_size_xyz[3];
    int num_bin_xyz[3];
    int* heads_lst;
    int* nexts_lst;

};  // class BinLinkedList








/**
 * @brief Construct a new Basic Structure Info< Coord Type>:: Basic Structure Info object
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
BasicStructureInfo<CoordType>::BasicStructureInfo() {
    this->num_atoms = 0;

    this->projected_lengths[0] = 0;
    this->projected_lengths[1] = 0;
    this->projected_lengths[2] = 0;

    this->interplanar_distances[0] = 0;
    this->interplanar_distances[1] = 0;
    this->interplanar_distances[2] = 0;
}


/**
 * @brief Construct a new Basic Structure Info< Coord Type>:: Basic Structure Info object
 * 
 * @tparam CoordType 
 * @param structure 
 */
template <typename CoordType>
BasicStructureInfo<CoordType>::BasicStructureInfo(Structure<CoordType> &structure) {
    this->num_atoms = structure.num_atoms;
    if (this->num_atoms == 0) {
        this->projected_lengths[0] = 0;
        this->projected_lengths[1] = 0;
        this->projected_lengths[2] = 0;

        this->interplanar_distances[0] = 0;
        this->interplanar_distances[1] = 0;
        this->interplanar_distances[2] = 0;
    } else {
        CoordType* projected_lengths = (CoordType*)structure.get_projected_lengths();
        CoordType* interplanar_distances_nonconst = (CoordType*)structure.get_interplanar_distances();
        
        this->projected_lengths[0] = projected_lengths[0];
        this->projected_lengths[1] = projected_lengths[1];
        this->projected_lengths[2] = projected_lengths[2];
        this->interplanar_distances[0] = interplanar_distances_nonconst[0];
        this->interplanar_distances[1] = interplanar_distances_nonconst[1];
        this->interplanar_distances[2] = interplanar_distances_nonconst[2];

        free(projected_lengths);
        free(interplanar_distances_nonconst);
    }
}


/**
 * @brief Construct a new Basic Structure Info< Coord Type>:: Basic Structure Info object
 * 
 * @tparam CoordType 
 * @param rhs 
 */
template <typename CoordType>
BasicStructureInfo<CoordType>::BasicStructureInfo(const BasicStructureInfo& rhs) {
    this->num_atoms = rhs.num_atoms;

    this->projected_lengths[0] = rhs.projected_lengths[0];
    this->projected_lengths[1] = rhs.projected_lengths[1];
    this->projected_lengths[2] = rhs.projected_lengths[2];
    this->interplanar_distances[0] = rhs.interplanar_distances[0];
    this->interplanar_distances[1] = rhs.interplanar_distances[1];
    this->interplanar_distances[2] = rhs.interplanar_distances[2];
}


/**
 * @brief Overload assignment operator.
 * 
 * @tparam CoordType 
 * @param rhs 
 * @return BasicStructureInfo<CoordType>& 
 */
template <typename CoordType>
BasicStructureInfo<CoordType>& BasicStructureInfo<CoordType>::operator=(const BasicStructureInfo& rhs) {
    this->num_atoms = rhs.num_atoms;

    this->projected_lengths[0] = rhs.projected_lengths[0];
    this->projected_lengths[1] = rhs.projected_lengths[1];
    this->projected_lengths[2] = rhs.projected_lengths[2];
    this->interplanar_distances[0] = rhs.interplanar_distances[0];
    this->interplanar_distances[1] = rhs.interplanar_distances[1];
    this->interplanar_distances[2] = rhs.interplanar_distances[2];
}


/**
 * @brief Destroy the Basic Structure Info< Coord Type>:: Basic Structure Info object
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
BasicStructureInfo<CoordType>::~BasicStructureInfo() {
    if (this->num_atoms != 0)
        this->num_atoms = 0;
}


/**
 * @brief Print out the information about `BasicStructureInfo`.
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void BasicStructureInfo<CoordType>::show() const{
    printf("num_atoms = %15d\n", this->num_atoms);
    if (this->num_atoms != 0) {
        printf("Projected Lengths : \n");
        printf("[%15f, %15f, %15f]\n", this->projected_lengths[0], this->projected_lengths[1], this->projected_lengths[2]);
        printf("Inter Planar Distances : \n");
        printf("[%15f, %15f, %15f]\n", this->interplanar_distances[0], this->interplanar_distances[1], this->interplanar_distances[2]);
    } else {
        printf("This is a null basic_structure_info.\n");
    }
}








/**
 * @brief Construct a new matersdk::Supercell<Coord Type>::Supercell object
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
matersdk::Supercell<CoordType>::Supercell() {
    this->structure = Structure<CoordType>();
    this->prim_structure_info = BasicStructureInfo<CoordType>();

    this->scaling_matrix[0] = 1;
    this->scaling_matrix[1] = 1;
    this->scaling_matrix[2] = 1;
    
    this->num_atoms = 0;
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
    this->prim_structure_info = BasicStructureInfo<CoordType>(structure);
    for (int ii=0; ii<3; ii++) {
        this->scaling_matrix[ii] = scaling_matrix[ii];
    }
    this->num_atoms = this->prim_structure_info.num_atoms * this->scaling_matrix[0] * this->scaling_matrix[1] * this->scaling_matrix[2];

    this->calc_prim_cell_idx_xyz();             // Assign `this->prim_cell_idx_xyz`
    this->calc_prim_cell_idx();                 // Assign `this->prim_cell_idx`
    this->owned_atom_idxs = (int*)malloc(sizeof(int) * this->prim_structure_info.num_atoms);
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
    this->prim_structure_info = rhs.prim_structure_info;
    for (int ii=0; ii<3; ii++) {
        this->scaling_matrix[ii] = rhs.scaling_matrix[ii];
    }
    this->num_atoms = rhs.num_atoms;

    this->prim_cell_idx = rhs.prim_cell_idx;
    this->prim_cell_idx_xyz[0] = rhs.prim_cell_idx_xyz[0];
    this->prim_cell_idx_xyz[1] = rhs.prim_cell_idx_xyz[1];
    this->prim_cell_idx_xyz[2] = rhs.prim_cell_idx_xyz[2];

    this->owned_atom_idxs = (int*)malloc(sizeof(int) * rhs.prim_structure_info.num_atoms);

    if (this->num_atoms != 0) {
        for (int ii=0; ii<rhs.prim_structure_info.num_atoms; ii++) {
            this->owned_atom_idxs[ii] = rhs.owned_atom_idxs[ii];
        }
    }

}


template <typename CoordType>
Supercell<CoordType>& Supercell<CoordType>::operator=(const Supercell<CoordType> &rhs) {
    if (this->num_atoms != 0) {
        free(this->owned_atom_idxs);
    }
    
    this->structure = rhs.structure;
    this->prim_structure_info = rhs.prim_structure_info;
    for (int ii=0; ii<3; ii++) {
        this->scaling_matrix[ii] = rhs.scaling_matrix[ii];
    }
    this->num_atoms = rhs.num_atoms;

    this->prim_cell_idx = rhs.prim_cell_idx;
    this->prim_cell_idx_xyz[0] = rhs.prim_cell_idx_xyz[0];
    this->prim_cell_idx_xyz[1] = rhs.prim_cell_idx_xyz[1];
    this->prim_cell_idx_xyz[2] = rhs.prim_cell_idx_xyz[2];

    if (rhs.num_atoms != 0) {
        this->owned_atom_idxs = (int*)malloc(sizeof(int) * rhs.prim_structure_info.num_atoms);
        for (int ii=0; ii<rhs.prim_structure_info.num_atoms; ii++) {
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
    for (int ii=0; ii<this->prim_structure_info.num_atoms; ii++) {
        this->owned_atom_idxs[ii] = this->prim_cell_idx * this->prim_structure_info.num_atoms + ii;
    }
}


/**
 * @brief Print out the information of `this` (class = Supercell)
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
void Supercell<CoordType>::show() const {
    printf("Prmitive Cell:\n");
    this->prim_structure_info.show();

    printf("\nSuper Cell:\n");
    this->structure.show();
    printf("prim_num_atoms = %15d\n", this->prim_structure_info.num_atoms);
    printf("prim_cell_idx = %15d\n", this->prim_cell_idx);
    printf("prim_cell_idx_xyz = [%4d, %4d, %4d]\n", this->prim_cell_idx_xyz[0], this->prim_cell_idx_xyz[1], this->prim_cell_idx_xyz[2]);
    if (this->num_atoms != 0)
        printf("owned_atom_idxs range = %15d ~ %15d\n", this->owned_atom_idxs[0], this->owned_atom_idxs[this->prim_structure_info.num_atoms-1]);
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
    return (const int)(this->prim_structure_info.num_atoms);
}


template <typename CoordType>
const int Supercell<CoordType>::get_num_atoms() const {
    return (const int)(this->structure.num_atoms);
}


template <typename CoordType>
const int* Supercell<CoordType>::get_owned_atom_idxs() const {
    return (const int*)this->owned_atom_idxs;
}


/**
 * @brief Construct a new Bin Linked List< Coord Type>:: Bin Linked List object
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
BinLinkedList<CoordType>::BinLinkedList() {

}


/**
 * @brief Construct a new Bin Linked List< Coord Type>:: Bin Linked List object
 * 
 * @tparam CoordType 
 * @param structure 
 * @param rcut 
 * @param bin_size_xyz 
 * @param pbc_xyz 
 */
template <typename CoordType>
BinLinkedList<CoordType>::BinLinkedList(Structure<CoordType>& structure, CoordType rcut, CoordType* bin_size_xyz, bool* pbc_xyz) {
    // Step 1. 计算 `scaling_matrix` -- 根据 `rcut` 和 `interplanar_distances`
    this->rcut = rcut;
    this->bin_size_xyz[0] = bin_size_xyz[0];
    this->bin_size_xyz[1] = bin_size_xyz[1];
    this->bin_size_xyz[2] = bin_size_xyz[2];
    CoordType* prim_interplanar_distances = (CoordType*)structure.get_interplanar_distances();
    int* scaling_matrix = (int*)malloc(sizeof(int) * 3);
    int* extending_matrix = (int*)malloc(sizeof(int) * 3);
    for (int ii=0; ii<3; ii++) {
        extending_matrix[ii] = std::ceil(rcut / prim_interplanar_distances[ii]);
        scaling_matrix[ii] = extending_matrix[ii] * 2 + 1;
        if (pbc_xyz[ii] == false) {
            scaling_matrix[ii] = 1;
        }
    }
    printf("[%f, %f, %f]\n", prim_interplanar_distances[0], prim_interplanar_distances[1], prim_interplanar_distances[2]);
    printf("[%d, %d, %d]\n", scaling_matrix[0], scaling_matrix[1], scaling_matrix[2]);
    free(extending_matrix);
    free(prim_interplanar_distances);

    // Step 2. 初始化 supercell
    Supercell<CoordType> supercell(structure, scaling_matrix);
    this->supercell = supercell;
    CoordType* projected_lengths = (CoordType*)supercell.structure.get_projected_lengths();
    for (int ii=0; ii<3; ii++) {
        printf("%f, %f\n", projected_lengths[ii], bin_size_xyz[ii]);
        this->num_bin_xyz[ii] = std::ceil( projected_lengths[ii] / bin_size_xyz[ii] );
    }
    printf("[%d, %d, %d]\n", this->num_bin_xyz[0], this->num_bin_xyz[1], this->num_bin_xyz[2]);
    free(projected_lengths);
}


/*
template <typename CoordType>
BinLinkedList<CoordType>::~BinLinkedList() {    
    if (this->supercell.num_atoms != 0) {
        free(this->heads_lst);
        free(this->nexts_lst);
    }
    this->supercell.~Supercell();
}
*/


/**
 * @brief 
 *          S1. 找到 prim_atom_idx 在 `center_cell` 中对应的 `atom_idx`
 *          S2. 利用 `atom_cart_coord` 计算 `atom` 所在的 `bin_idx`，并返回
 * 
 * @tparam CoordType 
 * @param prim_atom_idx 
 * @return int 
 */
template <typename CoordType>
int BinLinkedList<CoordType>::get_bin_idx(int prim_atom_idx) {
    // Step 1. 获取 `prim_atom_idx` 在 supercell 中对应的 `atom_idx`，并获取其坐标 `atom_cart_coord`    
    int atom_idx = prim_atom_idx + (this->supercell.prim_cell_idx * this->supercell.get_prim_num_atoms());
    CoordType* atom_cart_coord = (CoordType*)this->supercell.structure.get_cart_coords()[atom_idx];
    printf("%d: [%f, %f, %f]\n", atom_idx, atom_cart_coord[0], atom_cart_coord[1], atom_cart_coord[2]);
    
    // Step 2. 根据 `atom_cart_coord` 计算原子所属的 bin_idx
    int bin_idx_xyz[3];
    for (int ii=0; ii<3; ii++) {
        bin_idx_xyz[ii] = std::floor( [ii] / this->bin_size_xyz[ii] );
    }
    free(atom_cart_coord);

    return (
        bin_idx_xyz[0] +
        bin_idx_xyz[1] * this->num_bin_xyz[0] + 
        bin_idx_xyz[2] * this->num_bin_xyz[0] * this->num_bin_xyz[1]
    );
}

} // namespace: matersdk

#endif