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
    at::Tensor& types_tensor,
    at::Tensor& ircs_tensor)
{
    c10::TensorOptions options = c10::TensorOptions()
        .dtype(ircs_tensor.scalar_type())
        .device(ircs_tensor.device());
    int64_t nbs = this->_coeff_pair_combs.size();
    at::Tensor imtp_b_tensor = at::zeros({nbs}, options);
this->_coeff_pair_combs.show();
    for (int ii=0; ii<nbs; ii++) {
        std::vector<MtpMCoeffPair> tmp_coeff_pair_comb = this->_coeff_pair_combs.coeff_pair_combs()[ii];
        at::Tensor imtp_b_tensor_ii = torch::tensor(1, options);    // Be used to loop specific comb.
        for (int jj=0; jj<tmp_coeff_pair_comb.size(); jj++) {
            int64_t tmp_mu = tmp_coeff_pair_comb[jj].coeff_pair().first;
            int64_t tmp_nu = tmp_coeff_pair_comb[jj].coeff_pair().second;
            at::Tensor tmp_mtp_m_tensor = this->_mtp_m_list[0]->as<MtpMModule>()->forward(
                tmp_mu, 
                tmp_nu, 
                iidx, 
                ifirstneigh_tensor, 
                types_tensor, 
                ircs_tensor);

            int64_t dim1 = imtp_b_tensor_ii.dim();
            int64_t dim2 = tmp_mtp_m_tensor.dim();
printf("%3d, %3d: %3ld, %3ld:\n\t", ii, jj, tmp_mu, tmp_nu);
std::cout << imtp_b_tensor_ii.sizes() << ", " << tmp_mtp_m_tensor.sizes() << std::endl;
            if (dim1 >= dim2) {
                std::vector<int64_t> dims2;
                dims2.resize(dim2, 0);
                std::vector<int64_t> dims1;
                dims1.resize(dim1, 0);
                for (int kk=0; kk<dim1; kk++) {
                    dims2[kk] = kk;
                    dims1[kk] = dim1 - dim2 + kk;
                }
                imtp_b_tensor_ii = at::tensordot(
                    imtp_b_tensor_ii, tmp_mtp_m_tensor,
                    dims1, dims2);
            } 
            else {
                assert(dim1 == 0);
                imtp_b_tensor_ii = imtp_b_tensor * tmp_mtp_m_tensor;
            }
            //printf("[%3d, %3d]:\t", ii, jj);
            //std::cout << imtp_b_tensor_ii.sizes() << std::endl;
        }
    }
    //return imtp_b_tensor;
    return at::Tensor();
}

};  // namespace : mtp
};  // namespace : matersdk