#include <iterator>
#include <cassert>
#include "../include/mtpLevel.h"


namespace matersdk {
namespace mtp {

MtpMCoeffPair::MtpMCoeffPair()
{
    this->_coeff_pair = std::pair<int, int>(0, 0);
    this->_level = 0;
}

MtpMCoeffPair::MtpMCoeffPair(int mu, int nu)
{
    this->_coeff_pair.first = mu;
    this->_coeff_pair.second = nu;
    this->_calc_level();
}

MtpMCoeffPair::MtpMCoeffPair(std::pair<int, int> coeff_pair)
{
    this->_coeff_pair = coeff_pair;
    this->_calc_level();
}

MtpMCoeffPair::MtpMCoeffPair(const MtpMCoeffPair& rhs)
{
    this->_coeff_pair = rhs._coeff_pair;
    this->_calc_level();
}

MtpMCoeffPair::~MtpMCoeffPair()
{
    // Cleanup code here
    // In this case, there's nothing to clean up, so it's empty.
}

MtpMCoeffPair& MtpMCoeffPair::operator=(const MtpMCoeffPair& rhs)
{
    this->_coeff_pair = rhs._coeff_pair;
    this->_calc_level();
}

void MtpMCoeffPair::_calc_level() 
{
    this->_level = 2 + 4*this->_coeff_pair.first + this->_coeff_pair.second;
}

const int MtpMCoeffPair::level() const
{
    return this->_level;
}

std::pair<int, int> MtpMCoeffPair::coeff_pair() const
{
    return this->_coeff_pair;
}


MtpMCoeffPairCombs::MtpMCoeffPairCombs() {
    this->_max_level = 0;
    this->_coeff_pair_combs = std::vector<std::vector<MtpMCoeffPair>>();
}

MtpMCoeffPairCombs::MtpMCoeffPairCombs(int max_level) 
{
    this->_max_level = max_level;
    this->_build();
}

MtpMCoeffPairCombs::MtpMCoeffPairCombs(const MtpMCoeffPairCombs& rhs)
{
    this->_max_level = rhs._max_level;
    this->_coeff_pair_combs = rhs._coeff_pair_combs;
}

MtpMCoeffPairCombs::~MtpMCoeffPairCombs()
{
    // There is nothing to clean up.
}

MtpMCoeffPairCombs& MtpMCoeffPairCombs::operator=(const MtpMCoeffPairCombs& rhs)
{
    this->_max_level = rhs._max_level;
    this->_coeff_pair_combs = rhs._coeff_pair_combs;
}

/**
 * @brief Given a level of MtpM, Return all corresponding coeff_pair combinations
 * 
 * @param aim_level 
 * @param start_idx_mu 
 * @param start_idx_nu 
 * @return std::vector<std::vector<MtpMCoeffPair>> 
 */
std::vector<std::vector<MtpMCoeffPair>> MtpMCoeffPairCombs::get_all_schemes4lev(
    int aim_level,
    int start_idx_mu,
    int start_idx_nu)
{
    std::vector<std::vector<MtpMCoeffPair>> coeff_pair_combs_now;
    coeff_pair_combs_now.clear();
    
    if (aim_level == 0) {
        coeff_pair_combs_now.resize(1);
        coeff_pair_combs_now[0] = std::vector<MtpMCoeffPair>();
    } else {
        for (int tmp_mu=start_idx_mu; tmp_mu<aim_level; tmp_mu++) {
            for (int tmp_nu=start_idx_nu; tmp_nu<aim_level; tmp_nu++) {
                if ((aim_level - (2+4*tmp_mu+tmp_nu)) >= 0) {
                    std::vector<std::vector<MtpMCoeffPair>> coeff_pair_combs_before = 
                        MtpMCoeffPairCombs::get_all_schemes4lev(
                            aim_level - (2+4*tmp_mu+tmp_nu),
                            tmp_mu,
                            tmp_nu);
                    for (auto& tmp_coeff_pair_comb : coeff_pair_combs_before) {
                        tmp_coeff_pair_comb.push_back(MtpMCoeffPair(tmp_mu, tmp_nu));
                        coeff_pair_combs_now.push_back(tmp_coeff_pair_comb);
                    }
                }
            }
        }
    }
    return coeff_pair_combs_now;
}

std::vector<std::vector<MtpMCoeffPair>> MtpMCoeffPairCombs::get_contracted_combs(
    const std::vector<std::vector<MtpMCoeffPair>>& coeff_pair_combs)
{
    std::vector<std::vector<MtpMCoeffPair>> coeff_pair_combs_contracted;
    coeff_pair_combs_contracted.clear();
    for (int ii=0; ii<coeff_pair_combs.size(); ii++) {
        std::vector<MtpMCoeffPair> tmp_coeff_pair_comb = coeff_pair_combs[ii];
        if (tmp_coeff_pair_comb.size() == 0)
            continue;
        else if ((tmp_coeff_pair_comb.size() == 1) and (tmp_coeff_pair_comb[0].coeff_pair().second != 0))
            continue;
        else if ((tmp_coeff_pair_comb.size() == 1) and (tmp_coeff_pair_comb[0].coeff_pair().second == 0))
            coeff_pair_combs_contracted.push_back(tmp_coeff_pair_comb);
        else {
            int nu_no0 = tmp_coeff_pair_comb[0].coeff_pair().second;
            int tot_nu_but_no0 = 0;
            for (int jj=1; jj<tmp_coeff_pair_comb.size(); jj++) {
                tot_nu_but_no0 += tmp_coeff_pair_comb[jj].coeff_pair().second;
            }
            if (tot_nu_but_no0 == nu_no0) 
                coeff_pair_combs_contracted.push_back(tmp_coeff_pair_comb);
        }
    }
    return coeff_pair_combs_contracted;
}

void MtpMCoeffPairCombs::_build() 
{
    this->_coeff_pair_combs.clear();
    for (int ii=0; ii<=this->_max_level; ii++) {
        std::vector<std::vector<MtpMCoeffPair>> tmp_combs = get_all_schemes4lev(ii, 0, 0);
        std::vector<std::vector<MtpMCoeffPair>> tmp_contracted_combs = get_contracted_combs(tmp_combs);
        for (std::vector<MtpMCoeffPair> tmp_contracted_comb : tmp_contracted_combs)
            this->_coeff_pair_combs.push_back(tmp_contracted_comb);
    }
}

const std::vector<std::vector<MtpMCoeffPair>>& MtpMCoeffPairCombs::coeff_pair_combs() const
{
    return this->_coeff_pair_combs;
}

void MtpMCoeffPairCombs::show() const {
    int count = 0;
    for (auto& tmp_coeff_pair_comb : this->_coeff_pair_combs) {
        printf("Comb#%5d:\n\t", count);
        for (auto& tmp_coeff_pair : tmp_coeff_pair_comb) {
            printf("[%3d, %3d], ", tmp_coeff_pair.coeff_pair().first, tmp_coeff_pair.coeff_pair().second);
        }
        printf("\n");
        count++;
    }
    printf("Max_Level = %3d\n", this->_max_level);
}

};  // namespace : mtp
};  // namespace : matersdk
