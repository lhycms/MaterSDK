#ifndef MATERSDK_COORD_NUM_H
#define MATERSDK_COORD_NUM_H

#include <stdlib.h>
#include <vector>

#include "../../../io/publicLayer/include/structure.h"
#include "../../../io/publicLayer/include/neighborList.h"
#include "../../../../core/include/arrayUtils.h"
#include "../../../../core/include/vec3Operation.h"


namespace matersdk {
namespace sro {


template <typename CoordType>
class SublatticeCoordNum {
public:
    SublatticeCoordNum();

    SublatticeCoordNum(NeighborList<CoordType>& neighbor_list, int ntypes, int* sublattice);

    SublatticeCoordNum(const SublatticeCoordNum& rhs);

    SublatticeCoordNum& operator=(const SublatticeCoordNum& rhs);

    ~SublatticeCoordNum();

    void generate();

private:
    NeighborList<CoordType> neighbor_list;
    CoordType rcut;
    int ntypes;
    int* sublattice;    // sublattice = {Re, Nb} or {S, Se} 
};



template <typename CoordType>
SublatticeCoordNum<CoordType>::SublatticeCoordNum() {
    this->neighbor_list = NULL;
    this->rcut = 0;
    this->ntype = 0;
    this->sublattice = nullptr;
}



template <typename CoordType>
SublatticeCoordNum<CoordType>::SublatticeCoordNum(
                    NeighborList<CoordType>& neighbor_list, 
                    int ntypes,
                    int* sublattice)
{
    this->neighbor_list = neighbor_list;
    this->rcut = this->neighbor_list->get_rcut();
    this->ntypes = ntypes;
    this->sublattice = (int*)malloc(sizeof(int) * this->ntypes);
    for (int ii=0; ii<this->ntypes; ii++)
        this->sublattice[ii] = sublattice[ii];
}


template <typename CoordType>
SublatticeCoordNum<CoordType>::SublatticeCoordNum(const SublatticeCoordNum& rhs)
{
    this->neighbor_list = rhs.neighbor_list;
    this->rcut = rhs.rcut;
    this->ntypes = rhs.ntypes;
    if (this->ntypes != 0) {
        this->sublattice = (int*)malloc(sizeof(int) * this->ntypes);
        for (int ii=0; ii<this->ntypes; ii++)
            this->sublattice[ii] = rhs.sublattice[ii];
    } else 
        this->sublattice = nullptr;
}


template <typename CoordType>
SublatticeCoordNum<CoordType>& SublatticeCoordNum<CoordType>::operator=(const SublatticeCoordNum& rhs) 
{
    if (this->ntypes != 0)
        free(this->sublattice);
    
    this->neighbor_list = rhs.neighbor_list;
    this->rcut = rhs.rcut;
    this->ntypes = rhs.ntypes;
    if (this->ntypes != 0) {
        this->sublattice = (int*)malloc(sizeof(int) * this->ntypes);
        for (int ii=0; ii<this->ntypes; ii++) {
            this->sublattice[ii] = rhs.sublattice[ii];
        }
    }

    return *this;
}


template <typename CoordType>
SublatticeCoordNum<CoordType>::~SublatticeCoordNum() {
    if (this->ntypes != 0)
        free(this->sublattice);
}


}   // namespace: sro
}   // namespace: maetrsdk

#endif