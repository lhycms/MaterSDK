#ifndef MATERSDK_MTP_H
#define MATERSDK_MTP_H

#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <cmath>
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
                CoordType rcut_smooth,
                int hmju,
                CoordType rji);  // Highest mju (切比雪夫多项式的最高次数)

    void calc_rb_vals();

    CoordType deriv(int mju, const CoordType rji) const;

    void show_in_value() const;

    void show_in_deriv() const;


private:
    CoordType rcut = 0;
    CoordType rcut_smooth = 0;
    int hmju = 0;
    CoordType rji = 0;
    CoordType* rb_vals = nullptr;   // Store the rb_vals[0], rb_vals[1], speed up the calculation
};  // RadialBasisChebyshev




template <typename CoordType>
RadialBasisChebyshev<CoordType>::RadialBasisChebyshev() {
    this->rcut = 0;
    this->rcut_smooth = 0;
    this->hmju = 0;
    this->rji = 0;
    this->rb_vals = nullptr;
}


template <typename CoordType>
RadialBasisChebyshev<CoordType>::RadialBasisChebyshev(
                CoordType rcut,
                CoordType rcut_smooth,
                int hmju,
                CoordType rji) 
{
    this->rcut = rcut;
    this->rcut_smooth = rcut_smooth;
    assert(this->rcut > this->rcut_smooth);
    this->hmju = hmju;
    this->rji = rji;

    this->rb_vals = (CoordType*)malloc( sizeof(CoordType) * (this->hmju+1) );
    this->calc_rb_vals();
}


template <typename CoordType>
void RadialBasisChebyshev<CoordType>::calc_rb_vals() 
{
    CoordType ksi = (2*this->rji - (this->rcut+this->rcut_smooth)) / (this->rcut-this->rcut_smooth);
    printf("+++ %f\n", ksi);
    if (this->hmju == 0) {          // hmju = 0
        this->rb_vals[0] = 1;
        return ;
    }
    else if (this->hmju == 1) {     // hmju = 1
        this->rb_vals[0] = 1;
        this->rb_vals[1] = ksi;
        return ;
    }
    else {                          // hmju > 1
        this->rb_vals[0] = 1;
        this->rb_vals[1] = ksi;
        for (int ii=2; ii<(this->hmju+1); ii++) 
            this->rb_vals[ii] = 2 * ksi * this->rb_vals[ii-1] - this->rb_vals[ii-2];
        return ;
    }
}


template <typename CoordType>
CoordType RadialBasisChebyshev<CoordType>::deriv(int mju, const CoordType rji) const 
{

}


template <typename CoordType>
void RadialBasisChebyshev<CoordType>::show_in_value() const {
    printf("rcut = %5f\n", this->rcut);
    printf("rcut_smooth = %5f\n", this->rcut_smooth);
    printf("highest mju = %3d\n", this->hmju);
    printf("rji = %5f\n", this->rji);
    printf("ksi = %5f\n", (2*this->rji - this->rcut - this->rcut_smooth) / (this->rcut - this->rcut_smooth));
    
    printf("rb_vals = [");
    for (int ii=0; ii<(this->hmju+1); ii++)
        printf("%8f, ", this->rb_vals[ii]);
    printf("]\n");
}

template <typename CoordType>
void RadialBasisChebyshev<CoordType>::show_in_deriv() const {
    printf("rcut = %5f\n", this->rcut);
    printf("rcut_smooth = %5f\n", this->rcut_smooth);
    printf("highest mju = %3d\n", this->hmju);
}

}   // namespace : mtp
}   // namespace : matersdk

#endif