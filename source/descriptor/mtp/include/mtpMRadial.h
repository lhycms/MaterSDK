#ifndef MATERSDK_MTPM_RADIAL_H
#define MATERSDK_MTPM_RADIAL_H
#include <cmath>

namespace matersdk {
namespace mtp {

template <typename CoordType>
class MtpSwitchFunc1 {
public:
    MtpSwitchFunc1(CoordType rcut, CoordType rcut_smooth);

    CoordType get_result(CoordType distance_ij) const;

    CoordType get_deriv2r() const;

private:
    CoordType rcut = 0;
    CoordType rcut_smooth = 0;
};  // class SwitchFunction1


template <typename CoordType>
class MtpSwitchFunc2 {
public:
    MtpSwitchFunc2(CoordType rcut, CoordType rcut_smooth);

    CoordType get_result(CoordType distance_ij) const;

    CoordType get_deriv2r(CoordType distance_ij) const;

private:
    CoordType rcut = 0;
    CoordType rcut_smooth = 0;
}; // class MtpSwitchFunc2



/**
 * @brief Construct a new Mtp Switch Func 1< Coord Type>:: Mtp Switch Func 1 object
 * $$\xi = \frac{2*r_{ij} - (R_{cut} + R_{cut_smooth})}{R_{cut} - R_{cut_smooth}}$$
 * Scaling to [-1, 1]
 * @tparam CoordType 
 * @param rcut 
 * @param rcut_smooth 
 */
template <typename CoordType>
MtpSwitchFunc1<CoordType>::MtpSwitchFunc1(CoordType rcut, CoordType rcut_smooth)
{
    this->rcut = rcut;
    this->rcut_smooth = rcut_smooth;
}

template <typename CoordType>
CoordType MtpSwitchFunc1<CoordType>::get_result(CoordType distance_ij) const
{
    return (2 * distance_ij - (this->rcut + this->rcut_smooth)) / (this->rcut - this->rcut_smooth);
}

template <typename CoordType>
CoordType MtpSwitchFunc1<CoordType>::get_deriv2r() const
{
    return 2 / (this->rcut - this->rcut_smooth);
}

/**
 * @brief Construct a new Mtp Switch Func 2< Coord Type>:: Mtp Switch Func 2 object
 *      $$uu = /frac{r_ij - rcut}{rcut - rcut_smooth}$$
 *      vv = 
 *          1, r_ij < rcut_smooth
 *          uu^3(-6uu^2 + 15uu -10) + 1, rcut_smooth <= r_ij < rcut
 *          0, r_ij >= rcut
 * @tparam CoordType 
 * @param rcut 
 * @param rcut_smooth 
 */
template <typename CoordType>
MtpSwitchFunc2<CoordType>::MtpSwitchFunc2(CoordType rcut, CoordType rcut_smooth)
{
    this->rcut = rcut;
    this->rcut_smooth = rcut_smooth;
}

template <typename CoordType>
CoordType MtpSwitchFunc2<CoordType>::get_result(CoordType distance_ij) const
{
    CoordType result;
    CoordType uu = (distance_ij - this->rcut_smooth) / (this->rcut - this->rcut_smooth);

    if (distance_ij < this->rcut_smooth)
        result = 1;
    else if ((distance_ij>=this->rcut_smooth) && (distance_ij<this->rcut))
        result = std::pow(uu, 3) * (-6*std::pow(uu, 2) + 15*uu - 10) + 1;
    else
        result = 0;
    return result;
}

template <typename CoordType>
CoordType MtpSwitchFunc2<CoordType>::get_deriv2r(CoordType distance_ij) const
{
    CoordType deriv2r;
    CoordType uu = (distance_ij - this->rcut_smooth) / (this->rcut - this->rcut_smooth);

    if (distance_ij < this->rcut_smooth)
        deriv2r = 0;
    else if ((distance_ij>=this->rcut_smooth) && (distance_ij<this->rcut))
        deriv2r = 1 / (this->rcut - this->rcut_smooth) * ( -30*std::pow(uu, 4) + 60*std::pow(uu, 3) - 30*std::pow(uu, 2));
    else
        deriv2r = 0;
    return deriv2r;
}

};  // namespace : mtp
};  // namespace : matersdk

#endif
