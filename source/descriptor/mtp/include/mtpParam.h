#ifndef MATERSDK_MTP_PARAM_H
#define MATERSDK_MTP_PARAM_H
#include <vector>

namespace matersdk {
namespace mtp {
class MtpParam {
public:
    static void find_param(
        int mtp_level,
        int& alpha_moments_count,
        int& alpha_index_basic_count,
        std::vector<std::vector<int>>& alpha_index_basic,
        int& alpha_index_times_count,
        std::vector<std::vector<int>>& alpha_index_times,
        int& alpha_scalar_moments,
        std::vector<int>& alpha_moment_mapping);
};  // class : MtpParam

};  // namespace : mtp
};  // namespace : matersdk
#endif