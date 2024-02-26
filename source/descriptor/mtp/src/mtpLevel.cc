#include "../include/mtpLevel.h"


namespace matersdk {
namespace mtp {

MtpMCoeffPair::MtpMCoeffPair()
{
    this->coeff_pair = std::pair<int, int>(0, 0);
    this->_level = 0;
}

MtpMCoeffPair::MtpMCoeffPair(int mu, int nu)
{
    this->coeff_pair.first = mu;
    this->coeff_pair.second = nu;
    this->_calc_level();
}

MtpMCoeffPair::MtpMCoeffPair(std::pair<int, int> coeff_pair)
{
    this->coeff_pair = coeff_pair;
    this->_calc_level();
}

MtpMCoeffPair::MtpMCoeffPair(const MtpMCoeffPair& rhs)
{
    this->coeff_pair = rhs.coeff_pair;
    this->_calc_level();
}

MtpMCoeffPair::~MtpMCoeffPair()
{
    // Cleanup code here
    // In this case, there's nothing to clean up, so it's empty.
}

MtpMCoeffPair& MtpMCoeffPair::operator=(const MtpMCoeffPair& rhs)
{
    this->coeff_pair = rhs.coeff_pair;
    this->_calc_level();
}

void MtpMCoeffPair::_calc_level() 
{
    this->_level = 2 + 4*this->coeff_pair.first + this->coeff_pair.second;
}

const int MtpMCoeffPair::level() const
{
    return this->_level;
}

};  // namespace : mtp
};  // namespace : matersdk