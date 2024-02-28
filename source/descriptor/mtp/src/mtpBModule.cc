#include <torch/torch.h>
#include <vector>
#include "../include/mtpLevel.h"
#include "../include/mtpMModule.h"
#include "../include/mtpBModule.h"


namespace matersdk {
namespace mtp {

MtpBModuleImpl::MtpBModuleImpl(
    int64_t max_level,
    //int64_t nmus,
    int64_t ntypes,
    int64_t size,
    at::Tensor& rcuts_tensor)
{   
    this->_max_level = max_level;
    this->_coeff_pair_combs = MtpMCoeffPairCombs(max_level);
    this->_nmus = this->_coeff_pair_combs.nmus();
    this->_ntypes = ntypes;
    this->_size = size;
    this->_rcuts_tensor = rcuts_tensor;
    this->_mtp_m_list->push_back(
        MtpMModule(
            this->_nmus,
            this->_ntypes,
            this->_size,
            this->_rcuts_tensor));
    this->register_module("mtp_m_list", this->_mtp_m_list);
}

at::Tensor MtpBModuleImpl::forward(
    int64_t iidx,
    at::Tensor& ifirstneigh_tensor,
    at::Tensor& types,
    at::Tensor& ircs_tensor)
{
    /*
    c10::TensorOptions options = c10::TensorOptions()
        .dtype(irc_tensor.scalar_type())
        .dtype(irc_tensor.device());
    int64_t nbs = this->_coeff_pair_combs.size();
    at::Tensor imtp_b_tensor = at::zeros({nbs}, options);
    for (int ii=0; ii<nbs; ii++) {
        std::vector<MtpMCoeffPair> tmp_coeff_pair_comb = this->_coeff_pair_combs[ii];
        at::Tensor imtp_b_tensor_ii = torch::tensor(1, options);
        for (int jj=0; jj<tmp_coeff_pair_comb.size(); jj++) {

        }
    }
    return imtp_b_tensor;
    */
    return at::Tensor();
}

};  // namespace : mtp
};  // namespace : matersdk