#ifndef MATERSDK_MTP_LEVEL_H
#define MATERSDK_MTP_LEVEL_H
#include <utility>


namespace matersdk {
namespace mtp {

class MtpMCoeffPair {
public:
    MtpMCoeffPair();

    MtpMCoeffPair(int mu, int nu);

    MtpMCoeffPair(std::pair<int, int> coeff_pair);

    MtpMCoeffPair(const MtpMCoeffPair& rhs);

    ~MtpMCoeffPair();

    MtpMCoeffPair& operator=(const MtpMCoeffPair& rhs);

    void _calc_level();

    const int level() const;

private:
    std::pair<int, int> coeff_pair = std::pair<int, int>(0, 0);
    int _level = 0;
};  // class : MtpMCoeffPair

};  // namespace : mtp
};  // namespace : matersdk

#endif