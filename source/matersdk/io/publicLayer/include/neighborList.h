#ifndef MATERSDK_NEIGHBOR_LIST_H
#define MATERSDK_NEIGHBOR_LIST_H

#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include "./binLinkedList.h"


namespace matersdk {


template <typename CoordType>
class SingleNeighborListSortBasis {
public:
    SingleNeighborListSortBasis(CoordType* distances) : distances(distances)
    {}

    bool operator()(int index_i, int index_j) const {
        return (this->distances[index_i] < this->distances[index_j]);
    }

private:
    CoordType* distances;
};  // class: SingleNeighborListSortBasis


template <typename CoordType>
class SingleNeighborListArrangement {
public:
    SingleNeighborListArrangement(std::vector<int>& single_neighbor_list, int *new_indices)
    : single_neighbor_list(single_neighbor_list), new_indices(new_indices)
    {}


    std::vector<int> arrange() const {
        std::vector<int> new_single_neighbor_list(
                    this->single_neighbor_list.size(), -1);

        for (int ii=0; ii<this->single_neighbor_list.size(); ii++) {
            new_single_neighbor_list[ii] = this->single_neighbor_list[this->new_indices[ii]];
        }

        return new_single_neighbor_list;
    }


private:
    std::vector<int> single_neighbor_list;
    int* new_indices;
};  // class: SingleNeighborListArrangement


/**
 * @brief A `neighbor list` class
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
class NeighborList {
public:
    friend class BinLinkedList<CoordType>;
    
    NeighborList();

    NeighborList(Structure<CoordType>& structure, CoordType rcut, CoordType* bin_size_xyz, bool* pbc_xyz, bool sort=false);

    NeighborList(Structure<CoordType>& structure, CoordType rcut, bool* pbc_xyz, bool sort=false);

    NeighborList(const NeighborList& rhs);

    NeighborList& operator=(const NeighborList& rhs);

    ~NeighborList();

    void _build(bool sort=false);       // Populate `this->neighbor_list` (`std::vector<int>* this->neighbor_list = new std::vector<int>[this->num_atoms];`)

    const std::vector<int>* get_neighbor_lists() const;

    const int get_num_center_atoms() const;

    const CoordType get_rcut() const;

    void show_in_index() const;         // 展示 supercell 中的 atom_index

    void show_in_prim_index() const;    // 展示 primitive cell 中的 atom_index

    void show_in_an() const;            // 展示 atomic_number

    void show_in_distances() const;     // 展示 距中心原子的距离

    const BinLinkedList<CoordType>& get_binLinkedList() const;

    const int get_max_num_neigh_atoms() const;

    const int get_max_num_neigh_atoms_ss(int neigh_atomic_number) const; // `ss`: specified specie

    const int get_max_num_neigh_atoms_ssss(int center_atomic_number, int neigh_atomic_number) const;

    const int get_num_neigh_atoms(int prim_atom_idx) const;         // 由 `this->neighbor_lists[ii]` 得到 prim_atom_idx 的近邻原子的数目

    int* get_neigh_atomic_numbers(int prim_atom_idx) const;         // 由 `this->neighbor_lists[ii]` 得到 prim_atom_idx 的近邻原子的原子序数

    CoordType* get_neigh_distances(int prim_atom_idx) const;        // 由 `this->neighbor_lists[ii]` 得到 prim_atom_idx 的近邻原子的距中心原子的距离

    CoordType** get_neigh_relative_cart_coords(int prim_atom_idx) const;     // 由 `this->neighbor_lists[ii]` 得到 prim_atom_idx 的近邻原子的距中心原子的坐标 (r_j - r_i)

private:
    BinLinkedList<CoordType> bin_linked_list;
    int num_atoms = 0;                              // The number of atoms in primitive cell.
    std::vector<int>* neighbor_lists = nullptr;           // `num_atoms` 个 `neighbor_lists` 构成完整的 `neighbor_list`
    CoordType rcut = 0;
};  // class: NeighborList


/**
 * @brief Construct a new Neighbor List< Coord Type>:: Neighbor List object
 * 
 * @tparam CoordType 
 */
template <typename CoordType>
NeighborList<CoordType>::NeighborList() {
    this->bin_linked_list = BinLinkedList<CoordType>();
    this->num_atoms = 0;
    this->neighbor_lists = nullptr;
    this->rcut = 0;
}


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
NeighborList<CoordType>::NeighborList(Structure<CoordType>& structure, CoordType rcut, CoordType* bin_size_xyz, bool* pbc_xyz, bool sort) {
    assert(structure.get_num_atoms() > 0);
    
    this->bin_linked_list = BinLinkedList<CoordType>(structure, rcut, bin_size_xyz, pbc_xyz);
    this->num_atoms = this->bin_linked_list.get_supercell().get_prim_num_atoms();
    this->neighbor_lists = new std::vector<int>[this->num_atoms];
    this->rcut = rcut;
    this->_build(sort);     // Populate `this->neighbor_lists`
}


/**
 * @brief Construct a new Neighbor List< Coord Type>:: Neighbor List object
 * 
 * @tparam CoordType 
 * @param structure 
 * @param rcut 
 * @param pbc_xyz 
 * @param sort 
 */
template <typename CoordType>
NeighborList<CoordType>::NeighborList(Structure<CoordType>& structure, CoordType rcut, bool* pbc_xyz, bool sort) {
    assert(structure.get_num_atoms() > 0);

    CoordType bin_size_xyz[3];
    // bin_size_xyz[0] = bin_size_xyz[1] = bin_size_xyz[2] = rcut / 2;
    // this->bin_linked_list = BinLinkedList<CoordType>(structure, rcut, bin_size_xyz, pbc_xyz);
    this->bin_linked_list = BinLinkedList<CoordType>(structure, rcut, pbc_xyz);
    this->num_atoms = this->bin_linked_list.get_supercell().get_prim_num_atoms();
    this->neighbor_lists = new std::vector<int>[this->num_atoms];
    this->rcut = rcut;
    this->_build(sort);
}


/**
 * @brief Copy constructor
 * 
 * @tparam CoordType 
 * @param rhs 
 */
template <typename CoordType>
NeighborList<CoordType>::NeighborList(const NeighborList<CoordType>& rhs) {
    this->bin_linked_list = rhs.bin_linked_list;
    this->num_atoms = rhs.num_atoms;
    this->rcut = rhs.rcut;
    
    if (this->num_atoms != 0) {
        this->neighbor_lists = new std::vector<int>[this->num_atoms];
        for (int ii=0; ii<this->num_atoms; ii++) {  // 遍历中心原子
            this->neighbor_lists[ii].resize( rhs.neighbor_lists[ii].size(), -1 );
            for (int jj=0; jj<rhs.neighbor_lists[ii].size(); jj++) { // 遍历近邻原子
                this->neighbor_lists[ii][jj] = rhs.neighbor_lists[ii][jj];
            }
        }
    }
}


/**
 * @brief Overloading the assignment operator
 * 
 * @tparam CoordType 
 * @param rhs 
 * @return NeighborList<CoordType>& 
 */
template <typename CoordType>
NeighborList<CoordType>& NeighborList<CoordType>::operator=(const NeighborList<CoordType>& rhs) {
    this->bin_linked_list = rhs.bin_linked_list;
    this->num_atoms = rhs.num_atoms;
    this->rcut = rhs.rcut;

    if (this->num_atoms != 0) {
        this->neighbor_lists = new std::vector<int>[this->num_atoms];
        for (int ii=0; ii<this->num_atoms; ii++) {
            this->neighbor_lists[ii].resize( rhs.neighbor_lists[ii].size() , -1 );
            for (int jj=0; jj<rhs.neighbor_lists[ii].size(); jj++) {
                this->neighbor_lists[ii][jj] = rhs.neighbor_lists[ii][jj];
            }
        }
    }

    return *this;
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
 * @brief Populate `this->neighbor_lists` && Sort the `single_neighbor_list` according to `distances`
 * 
 * @tparam CoordType 
 * @param sort 是否按照距中心原子的距离排序
 */
template <typename CoordType>
void NeighborList<CoordType>::_build(bool sort) {
    // Step 1. 得到 primitive cell 的相关信息
    // Step 1.1. primitive cell 中的原子个数、primitive cell 的 `prim_cell_idx`
    int prim_cell_idx = this->bin_linked_list.get_supercell().get_prim_cell_idx();
    int prim_num_atoms = this->bin_linked_list.get_supercell().get_prim_num_atoms();
    assert(prim_num_atoms == this->num_atoms);

    // Step 2. Populate `this->neighbor_lists`
    for (int ii=0; ii<this->num_atoms; ii++) {  // `ii`: atom 在 primitive cell 中的 index
        this->neighbor_lists[ii] = this->bin_linked_list.get_neigh_atoms(ii);
    }
    
    // Step 3. Sort the `single_neighbor_list` according to distances
    // Step 3.1. 得到 supercell 中所有sites的坐标
    const CoordType** supercell_cart_coords = this->bin_linked_list.get_supercell().get_structure().get_cart_coords();

    // Step 3.2. Sort the `this->neighbor_lists[ii]` according to `distances`
    if (sort) {
        int tmp_num_neigh_atoms;        // 某中心原子近邻的原子数目
        int tmp_center_atom_idx;        // `center_atom_idx`: 中心 atom 在 supercell 中对应的 index
        int tmp_neigh_atom_idx;
        const CoordType* tmp_center_cart_coord; // 中心原子的笛卡尔坐标 
        const CoordType* tmp_neigh_cart_coord;  // 近邻原子的笛卡尔坐标
        CoordType* distances;   // neigh_atom 距 center_atom 的距离
        int* indices;           // neigh_atom 的 index

        for (int ii=0; ii<this->num_atoms; ii++) {
            // Step 3.2.1. 得到 ii 原子在 supercell 中对应的 `center_atom_index`
            tmp_center_atom_idx = ii + prim_cell_idx * prim_num_atoms;
            tmp_center_cart_coord = supercell_cart_coords[tmp_center_atom_idx];
            
            tmp_num_neigh_atoms = this->neighbor_lists[ii].size();
            distances = (CoordType*)malloc(sizeof(CoordType) * tmp_num_neigh_atoms);
            
            indices = (int*)malloc(sizeof(int) * tmp_num_neigh_atoms);
            for (int jj=0; jj<tmp_num_neigh_atoms; jj++) {
                indices[jj] = jj;
            }
            
            // Step 3.2.2. Calculate the distances
            for (int jj=0; jj<tmp_num_neigh_atoms; jj++) {
                tmp_neigh_atom_idx = this->neighbor_lists[ii][jj];
                tmp_neigh_cart_coord = supercell_cart_coords[tmp_neigh_atom_idx];

                distances[jj] = std::sqrt(
                    std::pow(tmp_neigh_cart_coord[0] - tmp_center_cart_coord[0], 2) + 
                    std::pow(tmp_neigh_cart_coord[1] - tmp_center_cart_coord[1], 2) + 
                    std::pow(tmp_neigh_cart_coord[2] - tmp_center_cart_coord[2], 2)
                );
            }

            // Step 3.2.3. Sort 
            std::sort(indices, indices + tmp_num_neigh_atoms, SingleNeighborListSortBasis<CoordType>(distances));
            SingleNeighborListArrangement<CoordType> snl_arrangement(this->neighbor_lists[ii], indices);
            this->neighbor_lists[ii] = snl_arrangement.arrange();

            // Step 3.2.4. Free `distances`
            free(distances);
            free(indices);
        }

    // Step . Free memory
    }
}


template <typename CoordType>
const std::vector<int>* NeighborList<CoordType>::get_neighbor_lists() const {
    return this->neighbor_lists;
}


template <typename CoordType>
const int NeighborList<CoordType>::get_num_center_atoms() const {
    return this->num_atoms;
}


template <typename CoordType>
const CoordType NeighborList<CoordType>::get_rcut() const {
    return this->rcut;
}


template <typename CoordType>
void NeighborList<CoordType>::show_in_index() const {
    if (this->num_atoms == 0) {
        printf("This is NULL NeighborList.\n");
    }
    else {
        for (int ii=0; ii<this->num_atoms; ii++) {
            printf("%4d:\t", ii);
            for (int jj: this->neighbor_lists[ii]) {
                printf("%4d, ", jj);
            }
            printf("\n");
        }
        printf("R_cutoff = %f\n", this->rcut);
    }
}


template <typename CoordType>
void NeighborList<CoordType>::show_in_prim_index() const {
    if (this->num_atoms == 0) {
        printf("This is NULL NeighborList.\n");
    } else {
        for (int ii=0; ii<this->num_atoms; ii++) {
            printf("%4d:\t", ii);
            for (int jj: this->neighbor_lists[ii]) {
                printf("%4d, ", jj % this->num_atoms);
            }
            printf("\n");
        }
        printf("R_cutoff = %f\n", this->rcut);
    }
}


template <typename CoordType>
void NeighborList<CoordType>::show_in_an() const {
    if (this->num_atoms == 0) {
        printf("This is NULL NeighborList.\n");
    } else {
        const int* supercell_atomic_numbers = this->bin_linked_list.get_supercell().get_structure().get_atomic_numbers();

        for (int ii=0; ii<this->num_atoms; ii++) {
            printf("%4d:\t", supercell_atomic_numbers[ii]);
            for (int jj: this->neighbor_lists[ii]) {
                printf("%4d, ", supercell_atomic_numbers[jj]);
            }
            printf("\n");
        }
    }
    printf("R_cutoff = %f\n", this->rcut);
}


template <typename CoordType>
void NeighborList<CoordType>::show_in_distances() const {
    if (this->num_atoms == 0) {
        printf("This is NULL NeighborList.\n");
    } else {
        const CoordType** supercell_cart_coords = this->bin_linked_list.get_supercell().get_structure().get_cart_coords();
        const int* supercell_atomic_numbers = this->bin_linked_list.get_supercell().get_structure().get_atomic_numbers();
        CoordType* center_cart_coord = (CoordType*)malloc(sizeof(CoordType) * 3);   // 中心原子的坐标
        CoordType* neigh_cart_coord = (CoordType*)malloc(sizeof(CoordType) * 3);    // 近邻原子的坐标
        const int prim_cell_idx = this->bin_linked_list.get_supercell().get_prim_cell_idx();
        CoordType tmp_distance;

        for (int ii=0; ii<this->num_atoms; ii++) {
            printf("%d:\t", supercell_atomic_numbers[ii]);
            center_cart_coord[0] = supercell_cart_coords[ii + this->num_atoms * prim_cell_idx][0];
            center_cart_coord[1] = supercell_cart_coords[ii + this->num_atoms * prim_cell_idx][1];
            center_cart_coord[2] = supercell_cart_coords[ii + this->num_atoms * prim_cell_idx][2];

            for (int jj: this->neighbor_lists[ii]) {
                neigh_cart_coord[0] = supercell_cart_coords[jj][0];
                neigh_cart_coord[1] = supercell_cart_coords[jj][1];
                neigh_cart_coord[2] = supercell_cart_coords[jj][2];
                
                tmp_distance = std::sqrt(
                    std::pow(neigh_cart_coord[0] - center_cart_coord[0], 2) + 
                    std::pow(neigh_cart_coord[1] - center_cart_coord[1], 2) + 
                    std::pow(neigh_cart_coord[2] - center_cart_coord[2], 2)
                );
                printf("%6.6f, ", tmp_distance);
            }
            printf("\n");
        }

        printf("R_cutoff = %f\n", this->rcut);
        free(center_cart_coord);
        free(neigh_cart_coord);
    }
}


template <typename CoordType>
const BinLinkedList<CoordType>& NeighborList<CoordType>::get_binLinkedList() const {
    return this->bin_linked_list;
}


/**
 * @brief 计算最大近邻原子数 （包含所有元素种类）
 * 
 * @tparam CoordType 
 * @return const int 
 */
template <typename CoordType>
const int NeighborList<CoordType>::get_max_num_neigh_atoms() const {
    int max_num_neigh_atoms = 0;
    for (int ii=0; ii<this->num_atoms; ii++) {
        if ( this->neighbor_lists[ii].size() > max_num_neigh_atoms ) {
            max_num_neigh_atoms = this->neighbor_lists[ii].size();
        }
    }
    return max_num_neigh_atoms;
}


/**
 * @brief 计算 元素为`neigh_atomic_number` 的最大近邻原子数 （包含某种元素 -- `neigh_atomioc_number`）
 * 
 * @tparam CoordType 
 * @param neigh_atomic_number 
 * @return const int 
 */
template <typename CoordType>
const int NeighborList<CoordType>::get_max_num_neigh_atoms_ss(int neigh_atomic_number) const {
    // Step 1. Initialize the value 
    // Step 1.1. 
    int max_num_neigh_atoms_ss = 0;             // return value
    int tmp_max_num_neigh_atoms_ss = 0;         // 在循环中，暂时存储某个中心原子的 `tmp_max_num_neigh_atoms_ss`
    std::vector<int> tmp_neigh_atomic_numbers;  // 将 `this->neighbor_lists[ii]` 的index换成 `atomic_number` 存储在 `std::vector<int>` 中
    // Step 1.2. 
    const int* supercell_atomic_numbers = this->bin_linked_list.get_supercell().get_structure().get_atomic_numbers();

    // Step 2. Update the `max_num_neigh_atoms_ss`
    for (int ii=0; ii<this->num_atoms; ii++) {
        tmp_neigh_atomic_numbers.clear();   // Clears all elements from the vector, making it empty
        // Step 2.1. Populate `tmp_neigh_atomic_numbers`.
        for (int supercell_atom_index: this->neighbor_lists[ii]) {
            tmp_neigh_atomic_numbers.push_back(supercell_atomic_numbers[supercell_atom_index]);
        }

        // Step 2.2. 
        tmp_max_num_neigh_atoms_ss = std::count(
                    tmp_neigh_atomic_numbers.begin(), 
                    tmp_neigh_atomic_numbers.end(), 
                    neigh_atomic_number
        );

        // Step 2.3. 
        if (tmp_max_num_neigh_atoms_ss > max_num_neigh_atoms_ss) 
            max_num_neigh_atoms_ss = tmp_max_num_neigh_atoms_ss;
    }
    
    return max_num_neigh_atoms_ss;
}


/**
 * @brief 指定 `atomic_number`, `neigh_atomic_number`，计算最大近邻原子数
 * 
 * @tparam CoordType 
 * @param center_atomic_number 
 * @param neigh_atomic_number 
 * @return const int 
 */
template <typename CoordType>
const int NeighborList<CoordType>::get_max_num_neigh_atoms_ssss(int center_atomic_number, int neigh_atomic_number) const {
    // Step 1. Initialize the value 
    // Step 1.1.
    int max_num_neigh_atoms_ssss = 0;
    int tmp_max_num_neigh_atoms_ssss = 0;
    std::vector<int> tmp_neigh_atomic_numbers;
    // Step 1.2. 
    const int* supercell_atomic_numbers = this->bin_linked_list.get_supercell().get_structure().get_atomic_numbers();

    // Step 2. Update `max_num_nrigh_atoms_ssss`
    for (int ii=0; ii<this->num_atoms; ii++) {
        tmp_neigh_atomic_numbers.resize(this->neighbor_lists[ii].size());
        if (supercell_atomic_numbers[ii] == center_atomic_number) {
            // Step 2.1. Populate `tmp_neigh_atomic_numbers`
            for (int jj=0; jj<this->neighbor_lists[ii].size(); jj++)
                tmp_neigh_atomic_numbers[jj] = supercell_atomic_numbers[this->neighbor_lists[ii][jj]];

            // Step 2.2. 
            tmp_max_num_neigh_atoms_ssss = std::count(
                        tmp_neigh_atomic_numbers.begin(),
                        tmp_neigh_atomic_numbers.end(),
                        neigh_atomic_number
            );

            // Step 2.3. 
            if (tmp_max_num_neigh_atoms_ssss > max_num_neigh_atoms_ssss)
                max_num_neigh_atoms_ssss = tmp_max_num_neigh_atoms_ssss;
        }
    }

    return max_num_neigh_atoms_ssss;
}


/**
 * @brief Get the number of neigh atoms given `prim_atom_idx`
 * 
 * @tparam CoordType 
 * @param prim_atom_idx 
 * @return const int 
 */
template <typename CoordType>
const int NeighborList<CoordType>::get_num_neigh_atoms(int prim_atom_idx) const {
    return this->neighbor_lists[prim_atom_idx].size();
}


/**
 * @brief Given `prim_atom_idx`，返回近邻的元素种类
 * 
 * @tparam CoordType 
 * @param prim_atom_idx 
 * @return int* 
 */
template <typename CoordType>
int* NeighborList<CoordType>::get_neigh_atomic_numbers(int prim_atom_idx) const {
    // Step 1. 
    // Step 1.1. prim_atom_idx 的近邻原子数
    int num_neigh_atoms = this->neighbor_lists[prim_atom_idx].size();
    // Step 1.2. 得到 supercell 的原子序数
    const int* supercell_atomic_numbers = this->bin_linked_list.get_supercell().get_structure().get_atomic_numbers();
    
    // Step 2. Populate `neigh_atomic_numbers`
    int* neigh_atomic_numbers = (int*)malloc(sizeof(int) * num_neigh_atoms);
    for (int ii=0; ii<num_neigh_atoms; ii++) {
        neigh_atomic_numbers[ii] = supercell_atomic_numbers[this->neighbor_lists[prim_atom_idx][ii]];
    }

    return neigh_atomic_numbers;
}



/**
 * @brief Get distances given `prim_atom_idx`
 * 
 * @tparam CoordType 
 * @param prim_atom_index 
 * @return CoordType* 
 */
template <typename CoordType>
CoordType* NeighborList<CoordType>::get_neigh_distances(int prim_atom_index) const {
    // Step 1. 
    // Step 1.1. prim_atom_idx 的近邻原子数
    int num_neigh_atoms = this->neighbor_lists[prim_atom_index].size();
    // Step 1.2. 得到 supercell 中所有原子的坐标
    const CoordType** supercell_cart_coords = this->bin_linked_list.get_supercell().get_structure().get_cart_coords();
    
    // Step 2. Populate `neigh_atomic_distances`
    CoordType* distances = (CoordType*)malloc(sizeof(CoordType) * num_neigh_atoms);
    int prim_cell_idx = this->bin_linked_list.get_supercell().get_prim_cell_idx();
    int prim_num_atoms = this->bin_linked_list.get_supercell().get_prim_num_atoms();
    int center_atom_idx = prim_atom_index + prim_cell_idx * prim_num_atoms;
    int neigh_atom_idx;
    const CoordType* center_cart_coord = supercell_cart_coords[center_atom_idx];
    const CoordType* neigh_cart_coord;

    for (int ii=0; ii<num_neigh_atoms; ii++) {  // 遍历 prim_atom_idx 的所有邻居原子
        neigh_atom_idx = this->neighbor_lists[prim_atom_index][ii];
        neigh_cart_coord = supercell_cart_coords[neigh_atom_idx];
        distances[ii] = std::sqrt(
                std::pow(neigh_cart_coord[0] - center_cart_coord[0], 2) + 
                std::pow(neigh_cart_coord[1] - center_cart_coord[1], 2) + 
                std::pow(neigh_cart_coord[2] - center_cart_coord[2], 2)
        );
    }

    // Step . Free memory

    return distances;
}


/**
 * @brief Get relative cart coords given `prim_atom_idx`.
 * 
 * @tparam CoordType 
 * @param prim_atom_idx 
 * @return CoordType** 
 */
template <typename CoordType>
CoordType** NeighborList<CoordType>::get_neigh_relative_cart_coords(int prim_atom_idx) const {
    // Step 1. 
    // Step 1.1. Get the number of neigh atoms
    int num_neigh_atoms = this->neighbor_lists[prim_atom_idx].size();
    // Step 1.2. Get the cart coords of atoms in supercell
    const CoordType** supercell_cart_coords = this->bin_linked_list.get_supercell().get_structure().get_cart_coords();

    // Step 2. Populate `relative_cart_coords`
    // Step 2.1. Allocate memory for `relative_cart_coords`
    CoordType** relative_cart_coords = (CoordType**)malloc(sizeof(CoordType*) * num_neigh_atoms);
    for (int ii=0; ii<num_neigh_atoms; ii++) 
        relative_cart_coords[ii] = (CoordType*)malloc(sizeof(CoordType) * 3);

    // Step 2.2. 
    int prim_cell_idx = this->bin_linked_list.get_supercell().get_prim_cell_idx();
    int prim_num_atoms = this->bin_linked_list.get_supercell().get_prim_num_atoms();
    int center_atom_idx = prim_atom_idx + prim_cell_idx * prim_num_atoms;
    int neigh_atom_idx;
    const CoordType* center_cart_coord = supercell_cart_coords[center_atom_idx];
    const CoordType* neigh_cart_coord;

    // Step 2.3. Populate `relative_cart_coords`
    for (int ii=0; ii<num_neigh_atoms; ii++) {
        neigh_atom_idx = this->neighbor_lists[prim_atom_idx][ii];
        neigh_cart_coord = supercell_cart_coords[neigh_atom_idx];

        relative_cart_coords[ii][0] = ( neigh_cart_coord[0] - center_cart_coord[0] );
        relative_cart_coords[ii][1] = ( neigh_cart_coord[1] - center_cart_coord[1] );
        relative_cart_coords[ii][2] = ( neigh_cart_coord[2] - center_cart_coord[2] );
    }

    // Step . Free memory

    return relative_cart_coords;
}


}; // namespace: matersdk


#endif