#ifndef MATERSDK_MTP_H
#define MATERSDK_MTP_H

#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include "../../../io/publicLayer/include/structure.h"
#include "../../../io/publicLayer/include/neighborList.h"
#include "../../../../core/include/vec3Operation.h"
#include "../../../../core/include/arrayUtils.h"


namespace matersdk {
namespace mtp {


template <typename CoordType>
class RadialBasisChebyshev {
public:
    RadialBasisChebyshev();

    RadialBasisChebyshev(
                CoordType rcut,
                CoordType rcut_smooth);

    void show_in_value(const CoordType r) const;

    void show_in_deriv(const CoordType r) const;


private:
    CoordType rcut = 0;
    CoordType rcut_smooth = 0;
};  // RadialBasisChebyshev




template <typename CoordType>
RadialBasisChebyshev<CoordType>::RadialBasisChebyshev() {
    this->rcut = 0;
    this->rcut_smooth = 0;
}


template <typename CoordType>
RadialBasisChebyshev<CoordType>::RadialBasisChebyshev(
                CoordType rcut,
                CoordType rcut_smooth) 
{
    this->rcut = rcut;
    this->rcut_smooth = rcut_smooth;
    assert(this->rcut > this->rcut_smooth);
}


template <typename CoordType>
void RadialBasisChebyshev<CoordType>::show_in_value(const CoordType r) const {
    printf("rcut = %5f\n", this->rcut);
    printf("rcut_smooth = %5f\n", this->rcut_smooth);
}

template <typename CoordType>
void RadialBasisChebyshev<CoordType>::show_in_deriv(const CoordType r) const {
    printf("rcut = %5f\n", this->rcut);
    printf("rcut_smooth = %5f\n", this->rcut_smooth);
}

}   // namespace : mtp
}   // namespace : matersdk

#endif