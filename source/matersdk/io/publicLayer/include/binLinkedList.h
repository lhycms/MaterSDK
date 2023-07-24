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

private:
    Structure<CoordType> structure;
    int scaling_matrix[3];      // 扩包倍数；x, y, z 方向上的 primitive_cell 个数
    int prim_cell_idx;          // primitive cell 对应的 cell index
    int prim_num_atoms;         // primitive cell 的元素数目
    int *owned_atom_idxs;       // 

}; // class: Supercell



template <typename CoordType>
matersdk::Supercell<CoordType>::Supercell(Structure<CoordType>& structure, int *scaling_matrix)
{
    this->structure = structure;
    for (int ii=0; ii<3; ii++) {
        this->scaling_matrix[ii] = scaling_matrix[ii];
    }

    this->structure.make_supercell(this->scaling_matrix);
    this->structure.show();
}



} // namespace: matersdk

#endif